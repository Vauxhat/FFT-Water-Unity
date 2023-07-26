using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Vauxhat;

using Complex = System.Numerics.Complex;

public class CascadeOcean : MonoBehaviour
{
    [System.Serializable]
    public class WaveData
    {
        [Min(0.0f)] public float gravity = 9.81f;
        [Min(0.0f)] public int patchSize = 128;
        [Min(0.0f)] public float repeatTime = 200.0f;
        [Min(0.0f)] public float minWavelength = 0.001f;
        [Min(0.0f)] public float phillipsConstant = 0.0002f;
    }

    // Render textures, used by compute shaders.
    public RenderTexture[] _fourierBuffer;
    public RenderTexture[] _timeSpectrum;
    public RenderTexture[] _waveDispersion;
    public RenderTexture[] _displacement;
    public RenderTexture[] _surfaceNormal;
    public RenderTexture[] _stationarySpectrum;

    // Normal textures, used by the CPU.
    public Texture2D _butterflyTexture;
    public Texture2D[] _gaussianNoise;

    // Compute shaders.
    public ComputeShader _stationarySpectrumShader;
    public ComputeShader _fourierTransformShader;
    public ComputeShader _timeSpectrumShader;
    public ComputeShader _dispersionShader;
    public ComputeShader _textureMergeShader;

    // Pre-computed wave variables.
    private int[] _reversedIndex;

    // Public variables for handling wave simulation.
    [Range(0.0f, 360.0f)] public float _windAngle = 0.0f;
    [Min(0.0f)] public float _windSpeed = 6.0f;
    [Min(0.0f)] public float _steepness = 0.0f;

    // Previous variable states, used to detect changes.
    private float _prevWindAngle = 0.0f;
    private float _prevWindSpeed = 6.0f;

    // Private wave variables, used only in code.
    private Vector2 _windDirection;
    private int _textureSize = 256;

    public WaveData[] _waveData = new WaveData[4];
    public WaveData[] _prevWaveData = new WaveData[4];

    private Material _material;

    // Start is called before the first frame update
    void Start()
    {
        // Initialse array of reversed indices.
        _reversedIndex = new int[_textureSize];

        // Calculate the number of bits which it takes to represent the numbers.
        int bits = (int)Mathf.Log(_textureSize, 2);

        // Calculate the reversed position for each index in the array.
        for (int i = 0; i < _reversedIndex.Length; i++)
        {
            _reversedIndex[i] = MathsExt.BitReverse(i, bits);
        }

        // Calculate wind direction based on angle. Angles should be easier to work with in editor.
        _windDirection = new Vector2(Mathf.Sin(Mathf.Deg2Rad * _windAngle), Mathf.Cos(Mathf.Deg2Rad * _windAngle));

        // Initialise butterfly texture.
        _butterflyTexture = new Texture2D(bits, _textureSize, TextureFormat.RGBAFloat, false);
        _butterflyTexture.filterMode = FilterMode.Point;
        UpdateButterfly();

        // Initialise textures.
        InitialiseTextures();

        // Initialise compute shaders.
        InitialiseTimeSpectrumShader();
        InitialiseTextureMergeShader();
        InitialiseFourierShader();

        for (int i = 0; i < 4; i++)
        {
            // Initialise wave dispersion.
            UpdateWaveDispersion(i);

            // Initialise stationary spectrum.
            UpdateStationarySpectrum(i);

            // Initialise time-dependent spectrum.
            UpdateTimeSpectrum(i);

            // Initialise heightfield.
            UpdateHeightfield(i);
        }

        // Store local reference to material.
        _material = this.gameObject.GetComponent<Renderer>().sharedMaterial;
    }

    private void InitialiseTextures()
    {
        int numCascades = 4;

        // Create new array for noise textures.
        _gaussianNoise = new Texture2D[numCascades];

        for (int i = 0; i < _gaussianNoise.Length; i++)
        {
            // Initialise gaussian noise texture.
            _gaussianNoise[i] = new Texture2D(_textureSize, _textureSize, TextureFormat.RGBAFloat, false);

            // Fill array with independent gaussian values.
            for (int y = 0; y < _textureSize; y++)
            {
                for (int x = 0; x < _textureSize; x++)
                {
                    float r, g, b, a;

                    // Generate two pairs of gaussian noise.
                    MathsExt.GaussianRandom(out r, out g);
                    MathsExt.GaussianRandom(out b, out a);

                    // Store gaussian variables in texture.
                    _gaussianNoise[i].SetPixel(x, y, new Color(r, g, b, a));
                }
            }

            // Apply changes to texture.
            _gaussianNoise[i].Apply();
        }

        // Create new array for cascade textures.
        _displacement = new RenderTexture[numCascades];
        _surfaceNormal = new RenderTexture[numCascades];
        _waveDispersion = new RenderTexture[numCascades];
        _stationarySpectrum = new RenderTexture[numCascades];

        for (int i = 0; i < numCascades; i++)
        {
            // Create dispacement texture.
            _displacement[i] = new RenderTexture(_textureSize, _textureSize, 0, RenderTextureFormat.ARGBFloat);
            _displacement[i].enableRandomWrite = true;
            _displacement[i].wrapMode = TextureWrapMode.Repeat;

            // Create surface normal texture.
            _surfaceNormal[i] = new RenderTexture(_textureSize, _textureSize, 0, RenderTextureFormat.ARGBFloat);
            _surfaceNormal[i].enableRandomWrite = true;
            _surfaceNormal[i].wrapMode = TextureWrapMode.Repeat;

            // Create wave dispersion texture.
            _waveDispersion[i] = new RenderTexture(_textureSize, _textureSize, 0, RenderTextureFormat.RFloat);
            _waveDispersion[i].enableRandomWrite = true;

            // Create stationary spectrum texture.
            _stationarySpectrum[i] = new RenderTexture(_textureSize, _textureSize, 0, RenderTextureFormat.ARGBFloat);
            _stationarySpectrum[i].enableRandomWrite = true;
        }

        // Initialise textures for time spectrum.
        _timeSpectrum = new RenderTexture[5];

        for (int i = 0; i < _timeSpectrum.Length; i++)
        {
            // Create empty texture, allow random read-write.
            _timeSpectrum[i] = new RenderTexture(_textureSize, _textureSize, 0, RenderTextureFormat.RGFloat);
            _timeSpectrum[i].enableRandomWrite = true;
        }

        // Initialise pingpong buffer for fourier shader.
        _fourierBuffer = new RenderTexture[2];

        for (int i = 0; i < _fourierBuffer.Length; i++)
        {
            // Create empty texture, allow random read-write.
            _fourierBuffer[i] = new RenderTexture(_textureSize, _textureSize, 0, RenderTextureFormat.RGFloat);
            _fourierBuffer[i].enableRandomWrite = true;
        }
    }

    private void InitialiseTimeSpectrumShader()
    {
        // Find kernel index for compute function.
        int kernelIndex = _timeSpectrumShader.FindKernel("ComputeTimeSpectrum");

        // Pass output textures to compute shader.
        _timeSpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_timeSpectrumX"), _timeSpectrum[0]);
        _timeSpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_timeSpectrumY"), _timeSpectrum[1]);
        _timeSpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_timeSpectrumZ"), _timeSpectrum[2]);

        _timeSpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_slopeX"), _timeSpectrum[3]);
        _timeSpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_slopeZ"), _timeSpectrum[4]);
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
    }

    private void InitialiseFourierShader()
    {
        // Find kernel index for compute function.
        int kernelIndex = _fourierTransformShader.FindKernel("ComputeFourierTransform");

        // Pass textures to compute shader.
        _fourierTransformShader.SetTexture(kernelIndex, Shader.PropertyToID("_buffer0"), _fourierBuffer[0]);
        _fourierTransformShader.SetTexture(kernelIndex, Shader.PropertyToID("_buffer1"), _fourierBuffer[1]);
    }

    // Update is called once per frame
    void Update()
    {
        bool windChange = false;

        // Check if shared wind speed or angle has been changed.
        if (!Mathf.Approximately(_prevWindSpeed, _windSpeed) || !Mathf.Approximately(_prevWindAngle, _windAngle))
        {
            // Calculate wind direction based on angle. Angles should be easier to work with in editor.
            _windDirection = new Vector2(Mathf.Sin(Mathf.Deg2Rad * _windAngle), Mathf.Cos(Mathf.Deg2Rad * _windAngle));

            // Update previous variables to current values.
            _prevWindAngle = _windAngle;
            _prevWindSpeed = _windSpeed;

            windChange = true;
        }

        // Loop through each cascade.
        for (int i = 0; i < _waveData.Length; i++)
        {
            // Check if gravity has been changed since last update.
            if (!Mathf.Approximately(_prevWaveData[i].gravity, _waveData[i].gravity) || !Mathf.Approximately(_prevWaveData[i].patchSize, _waveData[i].patchSize) || !Mathf.Approximately(_prevWaveData[i].repeatTime, _waveData[i].repeatTime))
            {
                // Re-calculate wave dispersion.
                UpdateWaveDispersion(i);

                // Update previous variables to current values.
                _prevWaveData[i].gravity = _waveData[i].gravity;
                _prevWaveData[i].repeatTime = _waveData[i].repeatTime;
            }

            // Check if global variables have changed since last update.
            if (windChange || !Mathf.Approximately(_prevWaveData[i].patchSize, _waveData[i].patchSize) || !Mathf.Approximately(_prevWaveData[i].minWavelength, _waveData[i].minWavelength) || !Mathf.Approximately(_prevWaveData[i].phillipsConstant, _waveData[i].phillipsConstant))
            {
                // Update stationary spectrum.
                UpdateStationarySpectrum(i);

                // Update previous variables to current values.
                _prevWaveData[i].minWavelength = _waveData[i].minWavelength;
                _prevWaveData[i].phillipsConstant = _waveData[i].phillipsConstant;
            }

            // Update previous patch size.
            _prevWaveData[i].patchSize = _waveData[i].patchSize;

            // Update time-dependent spectrum.
            UpdateTimeSpectrum(i);

            // Update heightfield.
            UpdateHeightfield(i);
        }

        // Check if material exists.
        if (_material)
        {
            // Pass heightfield texture to material shader.
            _material.SetTexture("_LOD0", _displacement[0]);
            _material.SetTexture("_LOD1", _displacement[1]);
            _material.SetTexture("_LOD2", _displacement[2]);
        }
    }

    private void UpdateStationarySpectrum(int index = 0)
    {
        uint x, y, z;

        // Find kernel index for compute function.
        int kernelIndex = _stationarySpectrumShader.FindKernel("ComputeStationary");

        // Pass textures to compute shader.
        _stationarySpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_gaussianNoise"), _gaussianNoise[index]);

        // Send input variables to shader.
        _stationarySpectrumShader.SetFloat(Shader.PropertyToID("_windX"), _windDirection.x);
        _stationarySpectrumShader.SetFloat(Shader.PropertyToID("_windY"), _windDirection.y);
        _stationarySpectrumShader.SetFloat(Shader.PropertyToID("_windSpeed"), _windSpeed);
        _stationarySpectrumShader.SetFloat(Shader.PropertyToID("_gravity"), _waveData[index].gravity);
        _stationarySpectrumShader.SetFloat(Shader.PropertyToID("_minWavelength"), _waveData[index].minWavelength);
        _stationarySpectrumShader.SetFloat(Shader.PropertyToID("_phillipsConstant"), _waveData[index].phillipsConstant);

        _stationarySpectrumShader.SetInt(Shader.PropertyToID("_textureSize"), _textureSize);
        _stationarySpectrumShader.SetInt(Shader.PropertyToID("_patchSize"), _waveData[index].patchSize);

        // Set output texture.
        _stationarySpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_stationarySpectrum"), _stationarySpectrum[index]);

        // Extract thread group size.
        _stationarySpectrumShader.GetKernelThreadGroupSizes(kernelIndex, out x, out y, out z);

        // Execute compute function.
        _stationarySpectrumShader.Dispatch(kernelIndex, (int)(_textureSize / x), (int)(_textureSize / y), 1);
    }

    private void UpdateWaveDispersion(int index = 0)
    {
        uint x, y, z;

        // Update shader variables.
        _dispersionShader.SetFloat(Shader.PropertyToID("_gravity"), _waveData[index].gravity);
        _dispersionShader.SetFloat(Shader.PropertyToID("_repeatTime"), _waveData[index].repeatTime);
        _dispersionShader.SetFloat(Shader.PropertyToID("_patchSize"), _waveData[index].patchSize);
        _dispersionShader.SetFloat(Shader.PropertyToID("_textureSize"), _textureSize);

        // Find kernel index for compute function.
        int kernelIndex = _dispersionShader.FindKernel("ComputeDispersion");

        // Set output texture.
        _dispersionShader.SetTexture(kernelIndex, Shader.PropertyToID("_waveDispersion"), _waveDispersion[index]);

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

    private void UpdateTimeSpectrum(int index = 0)
    {
        uint x, y, z;

        // Update shader variables.
        _timeSpectrumShader.SetFloat(Shader.PropertyToID("_time"), Time.time);
        _timeSpectrumShader.SetFloat(Shader.PropertyToID("_textureSize"), _textureSize);
        _timeSpectrumShader.SetFloat(Shader.PropertyToID("_patchSize"), _waveData[index].patchSize);

        // Find kernel index for compute function.
        int kernelIndex = _timeSpectrumShader.FindKernel("ComputeTimeSpectrum");

        // Pass textures to compute shader.
        _timeSpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_stationarySpectrum"), _stationarySpectrum[index]);
        _timeSpectrumShader.SetTexture(kernelIndex, Shader.PropertyToID("_waveDispersion"), _waveDispersion[index]);

        // Extract thread group size.
        _timeSpectrumShader.GetKernelThreadGroupSizes(kernelIndex, out x, out y, out z);

        // Execute compute function.
        _timeSpectrumShader.Dispatch(kernelIndex, (int)(_textureSize / x), (int)(_textureSize / y), 1);
    }

    private void UpdateHeightfield(int index = 0)
    {
        // Compute FFT for x, y, and z displacement.
        ComputeFFT(_timeSpectrum[0]);
        ComputeFFT(_timeSpectrum[1]);
        ComputeFFT(_timeSpectrum[2]);

        // Compute FFT for surface normal.
        ComputeFFT(_timeSpectrum[3]);
        ComputeFFT(_timeSpectrum[4]);

        // Merge real values for displacement into a single texture.
        MergeTextures(index);
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

    private void MergeTextures(int index = 0)
    {
        uint x, y, z;

        // Find kernel index for compute functions.
        int displacementIndex = _textureMergeShader.FindKernel("ComputeDisplacement");
        int normalIndex = _textureMergeShader.FindKernel("ComputeSurfaceNormal");

        // Set lambda.
        _textureMergeShader.SetFloat(Shader.PropertyToID("_lambda"), _steepness);

        // Set output textures.
        _textureMergeShader.SetTexture(displacementIndex, Shader.PropertyToID("_displacement"), _displacement[index]);
        _textureMergeShader.SetTexture(normalIndex, Shader.PropertyToID("_surfaceNormal"), _surfaceNormal[index]);

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