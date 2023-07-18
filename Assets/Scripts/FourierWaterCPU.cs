using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Vauxhat;

using Complex = System.Numerics.Complex;

public class FourierWaterCPU : MonoBehaviour
{
    // Texture variables.
    public Texture2D _fourierTexture;
    public Texture2D _heightfield;
    public Texture2D _butterflyTexture;

    // Pre-computed wave variables. Only re-generate if necessary.
    private Complex[] _gaussianNoise;
    private Complex[] _frequencySpectrum;
    private float[] _waveDispersion;
    private Complex[,] _twiddleFactor;
    //private Complex[] _twiddleFactor;
    private uint[] _reversedIndex;

    // Public variables for handling wave simulation.
    [Range(0.0f, 360.0f)] public float _windAngle = 0.0f;
    public float _windSpeed = 0.0f;
    public float _gravity = 9.81f;
    public int _patchSize = 10;
    public float _repeatTime = 200.0f;
    [Range(-16.0f, 16.0f)] public float _steepness = -1.0f;
    public float _minWavelength = 0.001f;

    // Private wave variables, used only in code.
    private Vector2 _windDirection;
    private int _textureSize = 128;
    private int _fourierStages;

    // Previous variable states, used to detect changes.
    private float _prevWindAngle = 0.0f;
    private float _prevWindSpeed = 0.0f;
    private float _prevGravity = 9.81f;
    private int _prevPatchSize = 10;
    private float _prevRepeatTime = 200.0f;
    private float _prevMinWavelength = 0.001f;

    private Material _material;

    // Start is called before the first frame update
    void Start()
    {
        // Initialise array of gaussian values, generating two per pixel.
        _gaussianNoise = new Complex[_textureSize * _textureSize * 2];

        // Fill array with independent gaussian values.
        for (int i = 0; i < _textureSize * _textureSize; i++)
        {
            _gaussianNoise[2 * i] = new Complex(MathsExt.BoxMullerTransform(), MathsExt.BoxMullerTransform());
            _gaussianNoise[2 * i + 1] = new Complex(MathsExt.BoxMullerTransform(), MathsExt.BoxMullerTransform());
        }

        // Initialise wave disperion, pre-compute for faster computation.
        _waveDispersion = new float[_textureSize * _textureSize];
        UpdateWaveDispersion();

        // Calculate wind direction based on angle. Angles should be easier to work with in editor.
        _windDirection = new Vector2(Mathf.Sin(Mathf.Deg2Rad * _windAngle), Mathf.Cos(Mathf.Deg2Rad * _windAngle));

        // Initialise frequency spectrum. Interweave ~h and ~h* for better locality.
        _frequencySpectrum = new Complex[_textureSize * _textureSize * 2];

        // Initialise fourier texture.
        _fourierTexture = new Texture2D(_textureSize, _textureSize, TextureFormat.RGBAHalf, false);
        UpdateFreqeuncySpectrum();

        // Initialse array of reversed indices.
        _reversedIndex = new uint[_textureSize];

        // Calculate the number of bits which it takes to represent the numbers.
        int bits = (int)Mathf.Log(_reversedIndex.Length, 2);

        // Calculate the reversed position for each index in the array.
        for (uint i = 0; i < _reversedIndex.Length; i++)
        {
            _reversedIndex[i] = MathsExt.BitReverse(i, bits);
        }

        // Initialise heightfield texture.
        _heightfield = new Texture2D(_textureSize, _textureSize, TextureFormat.RGBAHalf, false);

        // Store local reference to material.
        _material = this.gameObject.GetComponent<Renderer>().sharedMaterial;

        _fourierStages = (int)Mathf.Log(_textureSize, 2);

        // Initialise butterfly texture.
        _butterflyTexture = new Texture2D((int)Mathf.Log(_textureSize, 2), _textureSize, TextureFormat.RGBAHalf, false);
        _butterflyTexture.filterMode = FilterMode.Point;
        //_twiddleFactor = new Complex[_fourierStages * _textureSize];
        UpdateButterfly();

        // Initialise twiddle factors.
        _twiddleFactor = new Complex[_fourierStages, _textureSize];

        int pow2 = 1;

        for (int i = 0; i < _fourierStages; i++)
        {
            for (int j = 0; j < pow2; j++)
            {
                float exponent = (2.0f * Mathf.PI * j) / (pow2 * 2.0f);
                _twiddleFactor[i, j] = new Complex(Mathf.Cos(exponent), Mathf.Sin(exponent));
            }
            pow2 *= 2;
        }

        // Send butterfly texture to vertical and horizontal kernel.
        //_fourierShader.SetTexture(_fourierShader.FindKernel("VerticalPass"), Shader.PropertyToID("_butterflyTexture"), _butterflyTexture);
        //_fourierShader.SetTexture(_fourierShader.FindKernel("HorizontalPass"), Shader.PropertyToID("_butterflyTexture"), _butterflyTexture);
        //_fourierShader.SetTexture(_fourierShader.FindKernel("ComputeFFT"), Shader.PropertyToID("_butterflyTexture"), _butterflyTexture);
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
        if (!Mathf.Approximately(_prevWindAngle, _windAngle) || !Mathf.Approximately(_prevWindSpeed, _windSpeed) || !Mathf.Approximately(_prevPatchSize, _patchSize) || !Mathf.Approximately(_prevMinWavelength, _minWavelength))
        {
            // Calculate wind direction based on angle. Angles should be easier to work with in editor.
            _windDirection = new Vector2(Mathf.Sin(Mathf.Deg2Rad * _windAngle), Mathf.Cos(Mathf.Deg2Rad * _windAngle));
            //_windDirection = new Vector2(0.0f, 1.0f);

            // Update fourier texture.
            UpdateFreqeuncySpectrum();

            // Update previous variables to current values.
            _prevWindAngle = _windAngle;
            _prevWindSpeed = _windSpeed;
            _prevMinWavelength = _minWavelength;
        }

        // Update previous patch size.
        _prevPatchSize = _patchSize;

        // Update heightfield.
        UpdateHeightfield();

        // Check if material exists.
        if (_material)
        {
            // Pass heightfield texture to material shader.
            _material.SetTexture("_MainTex", _heightfield);
        }
    }

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

            // Define phillips constant.
            //float A = 0.0081f;
            float A = 0.0005f;

            // Calculate suppression. This should suppress waves at lower wavelenghts.
            float suppression = Mathf.Exp(-1.0f * k2 * l * l);

            // Calculate factor from phillips spectrum.
            phillipsSpectrum = A * (Mathf.Exp(-1.0f / kL2) / k4) * dotProduct * dotProduct * dotProduct * dotProduct * dotProduct * dotProduct * suppression;
        }

        return phillipsSpectrum;
    }

    private void UpdateFreqeuncySpectrum()
    {
        // The wave vector.
        Vector2 k;

        // Loop through each pixel in heightfield texture.
        for (int y = 0; y < _textureSize; y++)
        {
            // Calculate y componenet of vector k.
            //k.y = (2.0f * Mathf.PI * (y - _textureSize * 0.5f)) / _patchSize;
            //k.y = Mathf.PI * (2.0f * y - _textureSize) / _patchSize;
            k.y = (2.0f * Mathf.PI * y - Mathf.PI * _textureSize) / _patchSize;

            for (int x = 0; x < _textureSize; x++)
            {
                // Calculate x component of vector k.
                //k.x = (2.0f * Mathf.PI * (x - _textureSize * 0.5f)) / _patchSize;
                //k.x = Mathf.PI * (2.0f * x - _textureSize) / _patchSize;
                k.x = (2.0f * Mathf.PI * x - Mathf.PI * _textureSize) / _patchSize;

                // Calculate index in one dimensional array.
                int index = 2 * (x + y * _textureSize);

                // Calculate final value, sqrt(PhillipsSpectrum) and 1/sqrt(2) have been simplified into a single calculation.
                _frequencySpectrum[index] = Mathf.Sqrt(CalculatePhillips(k) * 0.5f) * _gaussianNoise[index];
                _frequencySpectrum[index + 1] = Complex.Conjugate(Mathf.Sqrt(CalculatePhillips(-k) * 0.5f) * _gaussianNoise[index + 1]);

                // Calculate final output, apply to noise texture.
                _fourierTexture.SetPixel(x, y, new Color((float)_frequencySpectrum[index].Real, (float)_frequencySpectrum[index].Imaginary, 0.0f, 1.0f));
                //_frequencyTexture.
            }
        }

        // Apply changes to texture.
        _fourierTexture.Apply();
    }

    void UpdateWaveDispersion()
    {
        // Loop through each pixel in the texture.
        for (int y = 0; y < _textureSize; y++)
        {
            // Get the y component of the wave vector.
            //float ky = (2.0f * Mathf.PI * (y - _textureSize * 0.5f)) / _patchSize;
            //float ky = Mathf.PI * (2.0f * y - _textureSize) / _patchSize;
            float ky = (2.0f * Mathf.PI * y - Mathf.PI * _textureSize) / _patchSize;

            for (int x = 0; x < _textureSize; x++)
            {
                // Get the x component of the wave vector.
                //float kx = (2.0f * Mathf.PI * (x - _textureSize * 0.5f)) / _patchSize;
                //float kx = Mathf.PI * (2.0f * x - _textureSize) / _patchSize;
                float kx = (2.0f * Mathf.PI * x - Mathf.PI * _textureSize) / _patchSize;

                // Calculate the length of the wave vector.
                float k_length = Mathf.Sqrt(kx * kx + ky * ky);

                // Calculate wave dispersion based on length and gravity.
                float dispersion = _gravity * k_length;

                // Set scale value to 1cm.
                float L = 0.01f;

                // Check if wave length is less than 1cm.
                if (k_length < L)
                {
                    // Modify dispersion relation, accomodate for surface tension.
                    dispersion *= (1 + k_length * k_length * L * L);
                }

                // Calculate final dispersion using squared dispersion.
                dispersion = Mathf.Sqrt(dispersion);

                // Calculate the basic wave frequency.
                float w_0 = (2.0f * Mathf.PI) / _repeatTime;

                // Calculate final dispersion value, ensuring wave frequency is a multiple of the base frequency (w_0).
                _waveDispersion[x + y * _textureSize] = Mathf.Floor(dispersion / w_0) * w_0;
            }
        }
    }

    private void UpdateHeightfield()
    {
        // Check if heightfield exists.
        if (_heightfield)
        {
            Complex[] height = new Complex[_textureSize * _textureSize];
            Complex[] slope_x = new Complex[_textureSize * _textureSize];
            Complex[] slope_y = new Complex[_textureSize * _textureSize];

            // Loop through each pixel in the heightfield.
            for (int y = 0; y < _textureSize; y++)
            {
                //float ky = (2.0f * Mathf.PI * (y - _textureSize * 0.5f)) / _patchSize;
                //float ky = Mathf.PI * (2.0f * (float)y - (float)_textureSize) / (float)_patchSize;
                float ky = (2.0f * Mathf.PI * y - Mathf.PI * _textureSize) / _patchSize;

                for (int x = 0; x < _textureSize; x++)
                {
                    //float kx = (2.0f * Mathf.PI * (x - _textureSize * 0.5f)) / _patchSize;
                    //float kx = Mathf.PI * (2.0f * (float)x - (float)_textureSize) / (float)_patchSize;
                    float kx = (2.0f * Mathf.PI * x - Mathf.PI * _textureSize) / _patchSize;

                    float wavelength = Mathf.Sqrt(kx * kx + ky * ky);

                    // Calculate index in one dimensional array.
                    int index = x + y * _textureSize;

                    // Calculate time factor.
                    float time_factor = _waveDispersion[index] * Time.time;

                    float cosine = Mathf.Cos(time_factor);
                    float sine = Mathf.Sin(time_factor);

                    // Calculate complex variables for ~h and the time exponent.
                    Complex tilde_h_0 = _frequencySpectrum[2 * index] * new Complex(Mathf.Cos(time_factor), Mathf.Sin(time_factor));
                    Complex tilde_h_0_i = _frequencySpectrum[2 * index + 1] * new Complex(Mathf.Cos(time_factor), Mathf.Sin(-time_factor));

                    Complex tilde_h = tilde_h_0 + tilde_h_0_i;

                    height[index] = tilde_h;

                    if (wavelength < 0.000001f)
                    {
                        slope_x[index] = new Complex(0.0f, 0.0f);
                        slope_y[index] = new Complex(0.0f, 0.0f);
                    }
                    else
                    {
                        slope_x[index] = tilde_h * new Complex(0.0f, -1.0f * kx / wavelength);
                        slope_y[index] = tilde_h * new Complex(0.0f, -1.0f * ky / wavelength);
                    }

                    // Set current pixel in heightfield texture.
                    _fourierTexture.SetPixel(x, y, new Color((float)tilde_h.Real, (float)tilde_h.Imaginary, 0.0f, 1.0f));
                }
            }

            // Apply changes to texture.
            _fourierTexture.Apply();

            // Horizontal pass.
            for (int y = 0; y < _textureSize; y++)
            {
                FastFourierTransform4(height, 1, y * _textureSize);
                FastFourierTransform4(slope_x, 1, y * _textureSize);
                FastFourierTransform4(slope_y, 1, y * _textureSize);
            }

            // Vertical pass.
            for (int x = 0; x < _textureSize; x++)
            {
                FastFourierTransform4(height, _textureSize, x);
                FastFourierTransform4(slope_x, _textureSize, x);
                FastFourierTransform4(slope_y, _textureSize, x);
            }

            float[] signs = { 1.0f, -1.0f };

            for (int y = 0; y < _textureSize; y++)
            {
                for (int x = 0; x < _textureSize; x++)
                {
                    int index = (x + y) % 2;
                    float sign = signs[index];

                    float lambda = 2.0f;

                    // Calculate final colour using sign and FFT result, divide by number of waves.
                    //float h = sign * ((float)height[x + y * _textureSize].Real / height.Length);
                    height[x + y * _textureSize] *= sign;
                    float h = (float)height[x + y * _textureSize].Real;

                    
                    float dx = sign * (float)slope_x[x + y * _textureSize].Real * _steepness;
                    float dy = sign * (float)slope_y[x + y * _textureSize].Real * _steepness;


                    // Set current pixel in heightfield texture.
                    _heightfield.SetPixel(x, y, new Color(h, dx, dy, 1.0f));
                }
            }

            // Apply changes to texture.
            _heightfield.Apply();
        }
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

        /*for (int y = 0; y < _textureSize; y++)
        {
            for (int x = 0; x < _fourierStages; x++)
            {
                // Calculate 2 to the power of x + 1, use bit shifting for simpler operation.
                int pow = 2 << x;

                // Calculate complex root of unity using k.
                int k = y * (_textureSize / pow) % _textureSize;
                float twiddleFactor = (2.0f * Mathf.PI * k) / _textureSize;

                _twiddleFactor[x + y * _fourierStages] = new Complex(Mathf.Cos(twiddleFactor), Mathf.Sin(twiddleFactor));
            }
        }*/
    }

    private void FastFourierTransform4(Complex[] input, int stride, int offset)
    {
        // Create a pingpong buffer of size N.
        Complex[,] buffer = new Complex[2, _textureSize];

        // Copy current row of input into buffer using bit-reversal permutation.
        for (int i = 0; i < _textureSize; i++)
        {
            buffer[0, i] = input[_reversedIndex[i] * stride + offset];
        }

        int pingpong = 0;

        // Perform FFT for each stage (zero to Log N).
        /*for (int s = 0; s < _butterfly.width; s++)
        {
            // Calculate 2 to the power of s.
            int m = 2 << s;

            // Loop through each index in the line, in groups of size m.
            for (int k = 0; k < _textureSize; k += m)
            {
                // loop through the first half of indices in the group.
                for (int j = 0; j < (m / 2); j++)
                {
                    // Get top and bottom index from sample.
                    int topIndex = k + j;
                    int bottomIndex = k + j + (m / 2);

                    // Get complex variables from previous stage.
                    Complex a = buffer[pingpong, topIndex];
                    Complex b = buffer[pingpong, bottomIndex] * _twiddleFactor[s + (k + j) * _fourierStages];

                    // Calculate new values for top and bottom index.
                    buffer[(pingpong + 1) % 2, topIndex] = a + b;
                    buffer[(pingpong + 1) % 2, bottomIndex] = a - b;
                }
            }

            // Switch to next array in buffer.
            pingpong = (pingpong + 1) % 2;
        }*/

        int loops = _textureSize >> 1;
        int size = 1 << 1;
        int size_over_2 = 1;
        int w_ = 0;

        for (int i = 1; i <= _fourierStages; i++)
        {
            pingpong ^= 1;
            for (int j = 0; j < loops; j++)
            {
                for (int k = 0; k < size_over_2; k++)
                {
                    buffer[pingpong, size * j + k] = buffer[pingpong ^ 1, size * j + k] +
                                  buffer[pingpong ^ 1, size * j + size_over_2 + k] * _twiddleFactor[w_, k];
                }

                for (int k = size_over_2; k < size; k++)
                {
                    buffer[pingpong, size * j + k] = buffer[pingpong ^ 1, size * j - size_over_2 + k] -
                                  buffer[pingpong ^ 1, size * j + k] * _twiddleFactor[w_, k - size_over_2];
                }
            }
            loops >>= 1;
            size <<= 1;
            size_over_2 <<= 1;
            w_++;
        }

        // Copy final FFT for the current sequence into the output buffer.
        for (int i = 0; i < _textureSize; i++)
        {
            input[i * stride + offset] = buffer[pingpong, i];
        }
    }
}