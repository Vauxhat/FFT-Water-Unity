using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Vauxhat;

using Complex = System.Numerics.Complex;

public class FourierWaterGPU : MonoBehaviour
{
    // Render textures, used by compute shaders.
    public RenderTexture[] _fourierBuffer;
    public RenderTexture[] _timeSpectrum;
    public RenderTexture _waveDispersion;
    public RenderTexture _displacement;
    public RenderTexture _surfaceNormal;
    public RenderTexture _stationarySpectrum;

    // Normal textures, used by the CPU.
    public Texture2D _butterflyTexture;
    public Texture2D _gaussianNoise;

    // Compute shaders.
    public ComputeShader _stationarySpectrumShader;
    public ComputeShader _fourierTransformShader;
    public ComputeShader _timeSpectrumShader;
    public ComputeShader _dispersionShader;
    public ComputeShader _textureMergeShader;

    // Pre-computed wave variables.
    private uint[] _reversedIndex;

    // Public variables for handling wave simulation.
    [Range(0.0f, 360.0f)] public float _windAngle = 0.0f;
    public float _windSpeed = 3.0f;
    public float _gravity = 9.81f;
    [Min(0.0f)] public int _patchSize = 10;
    [Min(0.0f)] public float _repeatTime = 200.0f;
    [Min(0.0f)] public float _steepness = 0.0f;
    [Min(0.0f)] public float _minWavelength = 0.001f;
    [Min(0.0f)] public float _phillipsConstant = 0.0002f;

    // Private wave variables, used only in code.
    private Vector2 _windDirection;
    private int _textureSize = 256;

    // Previous variable states, used to detect changes.
    private float _prevWindAngle = 0.0f;
    private float _prevWindSpeed = 0.0f;
    private float _prevGravity = 9.81f;
    private int _prevPatchSize = 10;
    private float _prevRepeatTime = 200.0f;
    private float _prevMinWavelength = 0.001f;
    private float _prevPhillipsConstant = 0.0002f;

    private Material _material;

    // Start is called before the first frame update
    void Start()
    {
        // Initialise gaussian noise.
        InitialiseGaussianNoise();

        // Initialse array of reversed indices.
        _reversedIndex = new uint[_textureSize];

        // Calculate the number of bits which it takes to represent the numbers.
        int bits = (int)Mathf.Log(_textureSize, 2);

        // Calculate the reversed position for each index in the array.
        for (uint i = 0; i < _reversedIndex.Length; i++)
        {
            _reversedIndex[i] = MathsExt.BitReverse(i, bits);
        }

        // Calculate wind direction based on angle. Angles should be easier to work with in editor.
        _windDirection = new Vector2(Mathf.Sin(Mathf.Deg2Rad * _windAngle), Mathf.Cos(Mathf.Deg2Rad * _windAngle));

        // Store local reference to material.
        _material = this.gameObject.GetComponent<Renderer>().sharedMaterial;

        // Initialise butterfly texture.
        _butterflyTexture = new Texture2D(bits, _textureSize, TextureFormat.RGBAFloat, false);
        _butterflyTexture.filterMode = FilterMode.Point;
        UpdateButterfly();

        // Initialise dispacement texture.
        _displacement = new RenderTexture(_textureSize, _textureSize, 0, RenderTextureFormat.ARGBFloat);
        _displacement.enableRandomWrite = true;
        _displacement.wrapMode = TextureWrapMode.Repeat;

        // Initialise surface normal texture.
        _surfaceNormal = new RenderTexture(_displacement);

        // Initialise compute shaders.
        InitialiseStationarySpectrumShader();
        InitialiseDispersionShader();
        InitialiseTimeSpectrumShader();
        InitialiseTextureMergeShader();
        InitialiseFourierShader();
    }

    private void InitialiseGaussianNoise()
    {
        // Initialise gaussian noise texture.
        _gaussianNoise = new Texture2D(_textureSize, _textureSize, TextureFormat.RGBAFloat, false);

        // Fill array with independent gaussian values.
        for (int i = 0; i < _textureSize * _textureSize; i++)
        {
            int x = i % _textureSize;
            int y = i / _textureSize;

            float r, g, b, a;

            // Generate two pairs of gaussian noise.
            MathsExt.GaussianRandom(0.0f, 1.0f, out r, out g);
            MathsExt.GaussianRandom(0.0f, 1.0f, out b, out a);

            // Store gaussian variables in texture.
            _gaussianNoise.SetPixel(x, y, new Color(r, g, b, a));
        }

        // Apply changes to texture.
        _gaussianNoise.Apply();
    }

    private void InitialiseStationarySpectrumShader()
    {
        // Find kernel index for compute function.
        int kernelIndex = _stationarySpectrumShader.FindKernel("ComputeStationary");

        // Inirialise stationary spectrum texture.
        _stationarySpectrum = new RenderTexture(_textureSize, _textureSize, 0, RenderTextureFormat.ARGBFloat);
        _stationarySpectrum.enableRandomWrite = true;

        // Pass textures to compute shader.
        _stationarySpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_stationarySpectrum"), _stationarySpectrum);
        //_stationarySpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_gaussianNoise"), _gaussianNoise);

        // Calculate initial spectrum.
        UpdateStationarySpectrum();
    }

    private void InitialiseDispersionShader()
    {
        // Find kernel index for compute function.
        int kernelIndex = _dispersionShader.FindKernel("ComputeDispersion");

        // Initialise wave dispersion texture.
        _waveDispersion = new RenderTexture(_textureSize, _textureSize, 0, RenderTextureFormat.RFloat);
        _waveDispersion.enableRandomWrite = true;

        // Pass texture to compute shader.
        _dispersionShader.SetTexture(kernelIndex, Shader.PropertyToID("_waveDispersion"), _waveDispersion);

        // Initialise global shader variables.
        _dispersionShader.SetFloat(Shader.PropertyToID("_gravity"), _gravity);
        _dispersionShader.SetFloat(Shader.PropertyToID("_repeatTime"), _repeatTime);
        _dispersionShader.SetInt(Shader.PropertyToID("_textureSize"), _textureSize);
        _dispersionShader.SetInt(Shader.PropertyToID("_patchSize"), _patchSize);

        // Calculate initial dispersion.
        UpdateWaveDispersion();
    }

    private void InitialiseTimeSpectrumShader()
    {
        // Find kernel index for compute function.
        int kernelIndex = _timeSpectrumShader.FindKernel("ComputeTimeSpectrum");

        // Pass input textures to compute shader.
        _timeSpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_stationarySpectrum"), _stationarySpectrum);
        _timeSpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_waveDispersion"), _waveDispersion);

        // Initialise textures for time spectrum.
        _timeSpectrum = new RenderTexture[5];

        _timeSpectrum[0] = new RenderTexture(_textureSize, _textureSize, 0, RenderTextureFormat.RGFloat);
        _timeSpectrum[0].enableRandomWrite = true;

        _timeSpectrum[1] = new RenderTexture(_timeSpectrum[0]);
        _timeSpectrum[2] = new RenderTexture(_timeSpectrum[0]);
        _timeSpectrum[3] = new RenderTexture(_timeSpectrum[0]);
        _timeSpectrum[4] = new RenderTexture(_timeSpectrum[0]);

        // Pass output textures to compute shader.
        _timeSpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_timeSpectrumX"), _timeSpectrum[0]);
        _timeSpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_timeSpectrumY"), _timeSpectrum[1]);
        _timeSpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_timeSpectrumZ"), _timeSpectrum[2]);

        _timeSpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_slopeX"), _timeSpectrum[3]);
        _timeSpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_slopeZ"), _timeSpectrum[4]);

        // Initialise global shader variables.
        _timeSpectrumShader.SetInt(Shader.PropertyToID("_textureSize"), _textureSize);
        _timeSpectrumShader.SetInt(Shader.PropertyToID("_patchSize"), _patchSize);
        _timeSpectrumShader.SetFloat(Shader.PropertyToID("_time"), Time.time);

        // Calculate initial time spectrum.
        UpdateTimeSpectrum();
    }

    private void InitialiseTextureMergeShader()
    {
        // Find kernel index for compute functions.
        int displacementKernal = _textureMergeShader.FindKernel("ComputeDisplacement");
        int normalKernel = _textureMergeShader.FindKernel("ComputeSurfaceNormal");

        // Send input textures to displacement kernel.
        _textureMergeShader.SetTexture(displacementKernal, Shader.PropertyToID("_displacementX"), _timeSpectrum[0]);
        _textureMergeShader.SetTexture(displacementKernal, Shader.PropertyToID("_displacementY"), _timeSpectrum[1]);
        _textureMergeShader.SetTexture(displacementKernal, Shader.PropertyToID("_displacementZ"), _timeSpectrum[2]);

        // Send input textures to normal kernel.
        _textureMergeShader.SetTexture(normalKernel, Shader.PropertyToID("_normalX"), _timeSpectrum[3]);
        _textureMergeShader.SetTexture(normalKernel, Shader.PropertyToID("_normalZ"), _timeSpectrum[4]);

        // Pass output textures to kernel.
        _textureMergeShader.SetTexture(displacementKernal, Shader.PropertyToID("_displacement"), _displacement);
        _textureMergeShader.SetTexture(normalKernel, Shader.PropertyToID("_surfaceNormal"), _surfaceNormal);
    }

    private void InitialiseFourierShader()
    {
        // Find kernel index for compute function.
        int kernelIndex = _fourierTransformShader.FindKernel("ComputeFourierTransform");

        // Pass butterfly texture to compute shader.
        _fourierTransformShader.SetTexture(kernelIndex, Shader.PropertyToID("_butterflyTexture"), _butterflyTexture);

        // Initialise pingpong buffer for fourier shader.
        _fourierBuffer = new RenderTexture[2];

        _fourierBuffer[0] = new RenderTexture(_textureSize, _textureSize, 0, RenderTextureFormat.RGFloat);
        _fourierBuffer[0].enableRandomWrite = true;

        _fourierBuffer[1] = new RenderTexture(_fourierBuffer[0]);

        // Pass textures to compute shader.
        _fourierTransformShader.SetTexture(kernelIndex, Shader.PropertyToID("_buffer0"), _fourierBuffer[0]);
        _fourierTransformShader.SetTexture(kernelIndex, Shader.PropertyToID("_buffer1"), _fourierBuffer[1]);

        // Calculate initial heightfield.
        UpdateHeightfield();
    }

    // Update is called once per frame
    void Update()
    {
        // Check if gravity has been changed since last update.
        if (!Mathf.Approximately(_prevGravity, _gravity) || !Mathf.Approximately(_prevPatchSize, _patchSize) || !Mathf.Approximately(_prevRepeatTime, _repeatTime))
        {
            // Re-calculate wave dispersion.
            UpdateWaveDispersion();

            // Update previous variables to current values.
            _prevGravity = _gravity;
            _prevRepeatTime = _repeatTime;
        }

        // Check if global variables have changed since last update.
        if (!Mathf.Approximately(_prevWindAngle, _windAngle) || !Mathf.Approximately(_prevWindSpeed, _windSpeed) || !Mathf.Approximately(_prevPatchSize, _patchSize) || !Mathf.Approximately(_prevMinWavelength, _minWavelength) || !Mathf.Approximately(_prevPhillipsConstant, _phillipsConstant))
        {
            // Calculate wind direction based on angle. Angles should be easier to work with in editor.
            _windDirection = new Vector2(Mathf.Sin(Mathf.Deg2Rad * _windAngle), Mathf.Cos(Mathf.Deg2Rad * _windAngle));

            // Update stationary spectrum.
            UpdateStationarySpectrum();

            // Update previous variables to current values.
            _prevWindAngle = _windAngle;
            _prevWindSpeed = _windSpeed;
            _prevMinWavelength = _minWavelength;
            _prevPhillipsConstant = _phillipsConstant;
        }

        // Update previous patch size.
        _prevPatchSize = _patchSize;

        // Update time-dependent spectrum.
        UpdateTimeSpectrum();

        // Update heightfield.
        UpdateHeightfield();

        // Check if material exists.
        if (_material)
        {
            // Pass heightfield texture to material shader.
            _material.SetTexture("_MainTex", _displacement);
        }
    }

    private void FixedUpdate()
    {
        
    }

    /*public Vector3 GetDisplacement(Vector3 worldPosition, float minWavelength)
    {
        // Calculate pixel coordinates of sample.
        int u = Mathf.FloorToInt(_textureSize * (worldPosition.x % _patchSize) / _patchSize);
        int v = Mathf.FloorToInt(_textureSize * (worldPosition.y % _patchSize) / _patchSize);

        // Sample stationary spectrum.
        Color stationarySample;



        // Calculate time dependent value.

        Complex dx, dy, dz;

        // Calculate wave vector k.
        float kx = Mathf.PI * (2.0f * u - _textureSize) / _patchSize;
        float ky = Mathf.PI * (2.0f * v - _textureSize) / _patchSize;

        // Calculate wavelength of k vector.
        float wavelength = Mathf.Sqrt(kx * kx + ky * ky);

        // Calculate time factor.
        float exponent = _waveDispersion[id] * Time.time;

        // Calculate ~h0 and conjugate using stationary spectrum, seed based on time.
        Complex tildeh0 = new Complex(stationarySample.r, stationarySample.g) * new Complex(Mathf.Cos(exponent), Mathf.Sin(exponent));
        Complex tildeh0i = new Complex(stationarySample.b, stationarySample.a) * new Complex(Mathf.Cos(exponent), Mathf.Sin(-exponent));

        // Calculate ~h from composite variables.
        Complex tildeh = tildeh0 + tildeh0i;

        // Set height component of time spectrum.
        dy = tildeh;

        // Make sure wavelength is a non-zero number, avoid division error.
        if (wavelength < 0.000001f)
        {
            // Set displacement to zero.
            dx = new Complex(0.0f, 0.0f);
            dz = new Complex(0.0f, 0.0f);
        }
        else
        {
            // Set x and y components of time spectrum.
            dx = tildeh * new Complex(0.0f, -kx / wavelength);
            dz = tildeh * new Complex(0.0f, -ky / wavelength);
        }

        // Calculate discrete fourier transform.
    }*/

    private float CalculatePhillips(Vector2 k)
    {
        // Initialise phillips spectrum.
        float phillipsSpectrum = 0.0f;

        // Calculate the length of the vector k.
        float waveLength = Mathf.Sqrt(k.x * k.x + k.y * k.y);

        // Set minimum wavelength (default 1mm).
        float l = _minWavelength;

        // Check if wavelength is greater than zero, avoid division error.
        if (waveLength > l)
        {
            // Calculate L factor (Squared wind speed divided by gravity).
            float L = (_windSpeed * _windSpeed) / _gravity;

            // Calculate normalized value for k, store as local variable.
            Vector2 waveDirection = k / waveLength;

            // Calculate wavelength to the power of four.
            float k2 = waveLength * waveLength;
            float k4 = k2 * k2;

            float kL2 = k2 * L * L;

            // Calculate the dot product of the wave direction and wind direction.
            float dotProduct = Vector2.Dot(waveDirection, _windDirection);
            dotProduct = dotProduct * dotProduct;

            // Define phillips constant.
            //float A = 0.0081f;
            float A = 0.0002f;

            // Calculate suppression. This should suppress waves at lower wavelenghts.
            float suppression = Mathf.Exp(-1.0f * k2 * l * l);

            // Calculate factor from phillips spectrum.
            phillipsSpectrum = A * (Mathf.Exp(-1.0f / kL2) / k4) * dotProduct * suppression;
        }

        return phillipsSpectrum;
    }

    private void UpdateStationarySpectrum()
    {
        uint x, y, z;

        // Find kernel index for compute function.
        int kernelIndex = _stationarySpectrumShader.FindKernel("ComputeStationary");

        // Pass textures to compute shader.
        _stationarySpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_gaussianNoise"), _gaussianNoise);

        // Send input variables to shader.
        _stationarySpectrumShader.SetFloat(Shader.PropertyToID("_windX"), _windDirection.x);
        _stationarySpectrumShader.SetFloat(Shader.PropertyToID("_windY"), _windDirection.y);
        _stationarySpectrumShader.SetFloat(Shader.PropertyToID("_windSpeed"), _windSpeed);
        _stationarySpectrumShader.SetFloat(Shader.PropertyToID("_gravity"), _gravity);
        _stationarySpectrumShader.SetFloat(Shader.PropertyToID("_minWavelength"), _minWavelength);
        _stationarySpectrumShader.SetFloat(Shader.PropertyToID("_phillipsConstant"), _phillipsConstant);

        _stationarySpectrumShader.SetInt(Shader.PropertyToID("_textureSize"), _textureSize);
        _stationarySpectrumShader.SetInt(Shader.PropertyToID("_patchSize"), _patchSize);

        // Extract thread group size.
        _stationarySpectrumShader.GetKernelThreadGroupSizes(kernelIndex, out x, out y, out z);

        // Execute compute function.
        _stationarySpectrumShader.Dispatch(kernelIndex, (int)(_textureSize / x), (int)(_textureSize / y), 1);
    }

    private void UpdateWaveDispersion()
    {
        uint x, y, z;

        // Update shader variables.
        _dispersionShader.SetFloat(Shader.PropertyToID("_gravity"), _gravity);
        _dispersionShader.SetFloat(Shader.PropertyToID("_repeatTime"), _repeatTime);
        _dispersionShader.SetFloat(Shader.PropertyToID("_patchSize"), _patchSize);
        _dispersionShader.SetFloat(Shader.PropertyToID("_textureSize"), _textureSize);

        // Find kernel index for compute function.
        int kernelIndex = _dispersionShader.FindKernel("ComputeDispersion");

        // Extract thread group size.
        _dispersionShader.GetKernelThreadGroupSizes(kernelIndex, out x, out y, out z);

        // Execute compute function.
        _dispersionShader.Dispatch(kernelIndex, (int)(_textureSize / x), (int)(_textureSize / y), 1);
    }

    // Update butterfly texture. Used to perform FFT calculations, only update if texture size changes.
    private void UpdateButterfly()
    {
        // Loop for height of butterfly texture (N).
        for (int y = 0; y < _butterflyTexture.height; y++)
        {
            // Loop for width of butterfly texture (Log N).
            for (int x = 0; x < _butterflyTexture.width; x++)
            {
                // Calculate 2 to the power of x + 1, use bit shifting for simpler operation.
                int pow = 2 << x;

                // Calculate complex root of unity using k.
                int k = y * (_textureSize / pow) % _textureSize;
                float twiddleFactor = (2.0f * Mathf.PI * (float)k) / (float)_textureSize;

                // Calculate butterfly span using 2 to the power of x, use bit shifting for simpler operation.
                int butterflySpan = 1 << x;

                if (x == 0)
                {
                    if ((y % pow) < butterflySpan)
                    {
                        _butterflyTexture.SetPixel(x, y, new Color(Mathf.Cos(twiddleFactor), Mathf.Sin(twiddleFactor), _reversedIndex[y], _reversedIndex[y + 1]));
                    }
                    else
                    {
                        _butterflyTexture.SetPixel(x, y, new Color(Mathf.Cos(twiddleFactor), Mathf.Sin(twiddleFactor), _reversedIndex[y - 1], _reversedIndex[y]));
                    }
                }
                else
                {
                    // Calculate bottom index based on top index.
                    if ((y % pow) < butterflySpan)
                    {
                        _butterflyTexture.SetPixel(x, y, new Color(Mathf.Cos(twiddleFactor), Mathf.Sin(twiddleFactor), y, y + butterflySpan));
                    }
                    else
                    {
                        _butterflyTexture.SetPixel(x, y, new Color(Mathf.Cos(twiddleFactor), Mathf.Sin(twiddleFactor), y - butterflySpan, y));
                    }
                }
            }
        }

        // Apply changes to texture.
        _butterflyTexture.Apply();
    }

    private void UpdateTimeSpectrum()
    {
        uint x, y, z;

        // Update shader variables.
        _timeSpectrumShader.SetFloat(Shader.PropertyToID("_time"), Time.time);
        _timeSpectrumShader.SetFloat(Shader.PropertyToID("_textureSize"), _textureSize);
        _timeSpectrumShader.SetFloat(Shader.PropertyToID("_patchSize"), _patchSize);

        // Find kernel index for compute function.
        int kernelIndex = _timeSpectrumShader.FindKernel("ComputeTimeSpectrum");

        // Pass textures to compute shader.
        _timeSpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_stationarySpectrum"), _stationarySpectrum);
        _timeSpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_waveDispersion"), _waveDispersion);

        // Extract thread group size.
        _timeSpectrumShader.GetKernelThreadGroupSizes(kernelIndex, out x, out y, out z);

        // Execute compute function.
        _timeSpectrumShader.Dispatch(kernelIndex, (int)(_textureSize / x), (int)(_textureSize / y), 1);
    }

    private void UpdateHeightfield()
    {
        // Compute FFT for x, y, and z displacement.
        ComputeFFT(_timeSpectrum[0]);
        ComputeFFT(_timeSpectrum[1]);
        ComputeFFT(_timeSpectrum[2]);

        // Compute FFT for surface normal.
        ComputeFFT(_timeSpectrum[3]);
        ComputeFFT(_timeSpectrum[4]);

        // Merge real values for displacement into a single texture.
        MergeTextures();

        // Merge real values for surface normal into a single texture.
    }

    private void ComputeFFT(RenderTexture texture)
    {
        // Set initial value for pingpong.
        int pingpong = 0;

        // Find kernel index for compute function.
        int kernelIndex = _fourierTransformShader.FindKernel("ComputeFourierTransform");

        // Copy data from source texture into ping pong buffer.
        Graphics.Blit(texture, _fourierBuffer[0]);

        // Set direction to horizontal.
        _fourierTransformShader.SetBool(Shader.PropertyToID("_direction"), true);

        uint x, y, z;
        _fourierTransformShader.GetKernelThreadGroupSizes(kernelIndex, out x, out y, out z);

        // Perform horizontal pass.
        for (int i = 0; i < _butterflyTexture.width; i++)
        {
            // Set pingpong value and index stage.
            _fourierTransformShader.SetBool(Shader.PropertyToID("_direction"), true);
            _fourierTransformShader.SetInt(Shader.PropertyToID("_pingpong"), pingpong);
            _fourierTransformShader.SetInt(Shader.PropertyToID("_stage"), i);

            // Pass texture to compute shader.
            _fourierTransformShader.SetTexture(kernelIndex, Shader.PropertyToID("_butterflyTexture"), _butterflyTexture);

            // Execute compute function.
            _fourierTransformShader.Dispatch(kernelIndex, (int)(_textureSize / x), (int)(_textureSize / y), 1);

            // Switch flag on pingpong variable.
            pingpong = (pingpong + 1) % 2;
        }

        // Set direction to vertical.
        _fourierTransformShader.SetBool(Shader.PropertyToID("_direction"), false);

        // Perform vertical pass.
        for (int i = 0; i < _butterflyTexture.width; i++)
        {
            // Set pingpong value and index stage.
            _fourierTransformShader.SetBool(Shader.PropertyToID("_direction"), false);
            _fourierTransformShader.SetInt(Shader.PropertyToID("_pingpong"), pingpong);
            _fourierTransformShader.SetInt(Shader.PropertyToID("_stage"), i);

            // Pass texture to compute shader.
            _fourierTransformShader.SetTexture(kernelIndex, Shader.PropertyToID("_butterflyTexture"), _butterflyTexture);

            // Execute compute function.
            _fourierTransformShader.Dispatch(kernelIndex, (int)(_textureSize / x), (int)(_textureSize / y), 1);

            // Switch flag on pingpong variable.
            pingpong = (pingpong + 1) % 2;
        }

        // Copy final output from buffer back into main texture.
        Graphics.Blit(_fourierBuffer[pingpong], texture);
    }

    private void MergeTextures()
    {
        uint x, y, z;

        // Find kernel index for compute functions.
        int displacementIndex = _textureMergeShader.FindKernel("ComputeDisplacement");
        int normalIndex = _textureMergeShader.FindKernel("ComputeSurfaceNormal");

        // Set lambda.
        _textureMergeShader.SetFloat(Shader.PropertyToID("_lambda"), _steepness);

        // Extract thread group size.
        _textureMergeShader.GetKernelThreadGroupSizes(displacementIndex, out x, out y, out z);

        // Execute compute function.
        _textureMergeShader.Dispatch(displacementIndex, (int)(_textureSize / x), (int)(_textureSize / y), 1);

        // Extract thread group size.
        _textureMergeShader.GetKernelThreadGroupSizes(normalIndex, out x, out y, out z);

        // Execute compute function.
        _textureMergeShader.Dispatch(normalIndex, (int)(_textureSize / x), (int)(_textureSize / y), 1);
    }
}