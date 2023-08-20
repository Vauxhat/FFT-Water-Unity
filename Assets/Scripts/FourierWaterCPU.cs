using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Jobs;
using Unity.Collections;
using Vauxhat;

using Complex = System.Numerics.Complex;

[System.Serializable]
public class FourierWaterCPU
{
    // Texture variables.
    public Texture2D _displacementTexture;
    private Vector3[] _displacement;

    // Pre-computed wave variables. Only re-generate if necessary.
    private Complex[] _gaussianNoise;
    private Complex[] _stationarySpectrum;
    private float[] _waveDispersion;
    private int[] _reversedIndex;

    // Public variables for handling wave simulation.
    [Range(0.0f, 360.0f)] public float _windAngle = 0.0f;
    [Min(0.0f)] public float _windSpeed = 6.0f;
    [Min(0.0f)] public float _gravity = 9.81f;
    [Min(0.0f)] public int _patchSize = 16;
    [Min(0.0f)] public float _repeatTime = 200.0f;
    [Min(0.0f)] public float _steepness = 0.0f;
    [Min(0.0f)] public float _minWavelength = 0.001f;
    [Min(0.0f)] public float _phillipsConstant = 0.0002f;

    // Private wave variables, used only in code.
    private Vector2 _windDirection;
    public int _textureSize = 128;

    // Previous variable states, used to detect changes.
    private float _prevWindAngle = 0.0f;
    private float _prevWindSpeed = 6.0f;
    private float _prevGravity = 9.81f;
    private int _prevPatchSize = 16;
    private float _prevRepeatTime = 200.0f;
    private float _prevMinWavelength = 0.001f;
    private float _prevPhillipsConstant = 0.0002f;

    //private Material _material;

    // Start is called before the first frame update
    public void Initialise()
    {
        // Initialise gaussian noise.
        //InitialiseGaussianNoise();

        // Pre-compute reversed index.
        _reversedIndex = FastFourierTransform.PrecomputeReversedIndex(_textureSize);

        // Initialise wave disperion, pre-compute for faster computation.
        _waveDispersion = new float[_textureSize * _textureSize];
        UpdateWaveDispersion();

        // Calculate wind direction based on angle. Angles should be easier to work with in editor.
        _windDirection = new Vector2(Mathf.Sin(Mathf.Deg2Rad * _windAngle), Mathf.Cos(Mathf.Deg2Rad * _windAngle));

        // Initialise frequency spectrum. Interweave ~h and ~h* for better locality.
        _stationarySpectrum = new Complex[_textureSize * _textureSize * 2];
        UpdateFrequencySpectrum();

        // Initialise displacement texture.
        _displacementTexture = new Texture2D(_textureSize, _textureSize, TextureFormat.RGBAHalf, false);
        _displacement = new Vector3[_textureSize * _textureSize];
        UpdateDisplacement(Time.time);
    }

    // Update is called once per frame
    public void Update(float time)
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
            UpdateFrequencySpectrum();

            // Update previous variables to current values.
            _prevWindAngle = _windAngle;
            _prevWindSpeed = _windSpeed;
            _prevMinWavelength = _minWavelength;
            _prevPhillipsConstant = _phillipsConstant;
        }

        // Update previous patch size.
        _prevPatchSize = _patchSize;

        // Update heightfield.
        UpdateDisplacement(time);
    }

    private void InitialiseGaussianNoise()
    {
        // Initialise array of gaussian values, generating two per pixel.
        _gaussianNoise = new Complex[_textureSize * _textureSize * 2];

        // Fill array with independent gaussian values.
        for (int i = 0; i < _textureSize * _textureSize; i++)
        {
            float r, g, b, a;

            // Generate two pairs of gaussian noise.
            MathsExt.GaussianRandom(out r, out g);
            MathsExt.GaussianRandom(out b, out a);

            // Store gaussian values in array.
            _gaussianNoise[2 * i] = new Complex(r, g);
            _gaussianNoise[2 * i + 1] = new Complex(b, a);
        }
    }

    private float CalculatePhillips(Vector2 k)
    {
        // Initialise phillips spectrum.
        float phillipsSpectrum = 0.0f;

        // Calculate the length of the vector k.
        float waveLength = k.magnitude;

        // Check if wavelength is greater than zero, avoid division error.
        if (waveLength > _minWavelength)
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

            // Calculate squared dot product, increasing the power increases directional dependency of waves.
            dotProduct = dotProduct * dotProduct * dotProduct * dotProduct * dotProduct * dotProduct;

            // Calculate suppression. This should suppress waves at lower wavelenghts.
            float suppression = Mathf.Exp(-1.0f * k2 * _minWavelength * _minWavelength);

            // Calculate factor from phillips spectrum.
            phillipsSpectrum = _phillipsConstant * (Mathf.Exp(-1.0f / kL2) / k4) * dotProduct * suppression;
        }

        // Return final value.
        return phillipsSpectrum;
    }

    private void UpdateFrequencySpectrum()
    {
        // The wave vector.
        Vector2 k;

        // Loop through each pixel in heightfield texture.
        for (int y = 0; y < _textureSize; y++)
        {
            // Calculate y componenet of wave vector.
            k.y = (2.0f * Mathf.PI * (y - _textureSize * 0.5f)) / _patchSize;

            for (int x = 0; x < _textureSize; x++)
            {
                // Calculate x component of wave vector.
                k.x = (2.0f * Mathf.PI * (x - _textureSize * 0.5f)) / _patchSize;

                // Calculate index in one dimensional array.
                int index = 2 * (x + y * _textureSize);

                // Calculate final value, sqrt(PhillipsSpectrum) and 1/sqrt(2) have been simplified into a single calculation.
                _stationarySpectrum[index] = Mathf.Sqrt(CalculatePhillips(k) * 0.5f) * _gaussianNoise[index];
                _stationarySpectrum[index + 1] = Complex.Conjugate(Mathf.Sqrt(CalculatePhillips(-k) * 0.5f) * _gaussianNoise[index + 1]);
            }
        }
    }

    private void UpdateWaveDispersion()
    {
        // Loop through each pixel in the texture.
        for (int y = 0; y < _textureSize; y++)
        {
            // Get the y component of the wave vector.
            float ky = (2.0f * Mathf.PI * (y - _textureSize * 0.5f)) / _patchSize;

            for (int x = 0; x < _textureSize; x++)
            {
                // Get the x component of the wave vector.
                float kx = (2.0f * Mathf.PI * (x - _textureSize * 0.5f)) / _patchSize;

                // Calculate the length of the wave vector.
                float kLength = Mathf.Sqrt(kx * kx + ky * ky);

                // Calculate wave dispersion based on length and gravity.
                float dispersion = _gravity * kLength;

                // Set scale value to 1cm.
                float scale = 0.01f;

                // Check if wave length is less than scale.
                if (kLength < scale)
                {
                    // Modify dispersion relation, accomodate for surface tension.
                    dispersion *= (1.0f + kLength * kLength * scale * scale);
                }

                // Extract dispersion (w) from squared value.
                dispersion = Mathf.Sqrt(dispersion);

                // Calculate the basic wave frequency.
                float baseFrequency = (2.0f * Mathf.PI) / _repeatTime;

                // Calculate final dispersion value, ensuring wave frequency is a multiple of the base frequency.
                _waveDispersion[x + y * _textureSize] = Mathf.Floor(dispersion / baseFrequency) * baseFrequency;
            }
        }
    }

    private struct TimeSpectrumJob : IJobParallelFor
    {
        // Input buffers.
        [ReadOnly] public NativeArray<float> _waveDispersion;
        [ReadOnly] public NativeArray<Complex> _stationarySpectrum;

        // Output buffers.
        public NativeArray<Complex> _dx;
        public NativeArray<Complex> _dy;
        public NativeArray<Complex> _dz;

        // Read only input variables.
        [ReadOnly] public float _time;
        [ReadOnly] public int _textureSize;
        [ReadOnly] public int _patchSize;

        public void Execute(int index)
        {
            // Calculate 2D coordinates based on 1D index.
            int x = index % _textureSize;
            int y = index / _textureSize;

            // Calculate wave vector.
            float ky = (2.0f * Mathf.PI * (y - _textureSize * 0.5f)) / _patchSize;
            float kx = (2.0f * Mathf.PI * (x - _textureSize * 0.5f)) / _patchSize;

            // Calculate wavelength of wave vector.
            float wavelength = Mathf.Sqrt(kx * kx + ky * ky);

            // Calculate time factor.
            float exponent = _waveDispersion[index] * _time;

            // Calculate ~h0 and conjugate using stationary spectrum, seed based on time.
            Complex tildeH0 = _stationarySpectrum[2 * index] * new Complex(Mathf.Cos(exponent), Mathf.Sin(exponent));
            Complex tildeH0Conjugate = _stationarySpectrum[2 * index + 1] * new Complex(Mathf.Cos(exponent), Mathf.Sin(-exponent));

            // Calculate ~h from composite variables.
            Complex tildeH = tildeH0 + tildeH0Conjugate;

            // Set height component of time dependent spectrum.
            _dy[index] = tildeH;

            // Make sure wavelength is a non-zero number, avoid division error.
            if (wavelength < 0.000001f)
            {
                // Set displacement to zero.
                _dx[index] = new Complex(0.0f, 0.0f);
                _dz[index] = new Complex(0.0f, 0.0f);
            }
            else
            {
                // Set x and y components of time spectrum.
                _dx[index] = tildeH * new Complex(0.0f, -kx / wavelength);
                _dz[index] = tildeH * new Complex(0.0f, -ky / wavelength);
            }
        }
    }

    private void UpdateDisplacement(float time)
    {
        // Check if texture exists.
        if (_displacementTexture)
        {
            // Create arrays for time dependent spectrum.
            Complex[] dx = new Complex[_textureSize * _textureSize];
            Complex[] dy = new Complex[_textureSize * _textureSize];
            Complex[] dz = new Complex[_textureSize * _textureSize];

            // Create native array for input data.
            NativeArray<float> nativeWaveDispersion = new NativeArray<float>(_waveDispersion.Length, Allocator.TempJob);
            NativeArray<Complex> nativeFrequencySpectrum = new NativeArray<Complex>(_stationarySpectrum.Length, Allocator.TempJob);

            // Create native array for output data.
            NativeArray<Complex> nativeX = new NativeArray<Complex>(_textureSize * _textureSize, Allocator.TempJob);
            NativeArray<Complex> nativeY = new NativeArray<Complex>(_textureSize * _textureSize, Allocator.TempJob);
            NativeArray<Complex> nativeZ = new NativeArray<Complex>(_textureSize * _textureSize, Allocator.TempJob);

            // Copy source array into native array.
            for (int i = 0; i < _waveDispersion.Length; i++)
            {
                nativeWaveDispersion[i] = _waveDispersion[i];
                nativeFrequencySpectrum[2 * i] = _stationarySpectrum[2 * i];
                nativeFrequencySpectrum[2 * i + 1] = _stationarySpectrum[2 * i + 1];
            }

            // Create new job.
            TimeSpectrumJob job = new TimeSpectrumJob()
            {
                _waveDispersion = nativeWaveDispersion,
                _stationarySpectrum = nativeFrequencySpectrum,
                _dx = nativeX,
                _dy = nativeY,
                _dz = nativeZ,
                _time = time,
                _textureSize = _textureSize,
                _patchSize = _patchSize,
            };

            // Create job handle.
            JobHandle jobHandle;

            // Perform horizontal operation.
            jobHandle = job.Schedule(_textureSize * _textureSize, 1);

            // Wait until all threads are finished.
            jobHandle.Complete();

            // Copy output array to source array.
            for (int i = 0; i < _textureSize * _textureSize; i++)
            {
                dx[i] = nativeX[i];
                dy[i] = nativeY[i];
                dz[i] = nativeZ[i];
            }

            // Dispose of input array structures.
            nativeWaveDispersion.Dispose();
            nativeFrequencySpectrum.Dispose();

            // Dispose of output array structures.
            nativeX.Dispose();
            nativeY.Dispose();
            nativeZ.Dispose();

            // Perform FFT operations.
            FastFourierTransform.Compute(dx, _textureSize, _reversedIndex);
            FastFourierTransform.Compute(dy, _textureSize, _reversedIndex);
            FastFourierTransform.Compute(dz, _textureSize, _reversedIndex);

            float[] signs = { 1.0f, -1.0f };

            // Loop through each pixel in the texture.
            for (int y = 0; y < _textureSize; y++)
            {
                for (int x = 0; x < _textureSize; x++)
                {
                    // Determine sign of current sample.
                    float sign = signs[(x + y) % 2];

                    // Calculate final displacement values.
                    float r = sign * (float)dx[x + y * _textureSize].Real * -1.0f * _steepness;
                    float g = sign * (float)dy[x + y * _textureSize].Real;
                    float b = sign * (float)dz[x + y * _textureSize].Real * -1.0f * _steepness;

                    // Store output in displacement texture.
                    //_displacementTexture.SetPixel(x, y, new Color(r, g, b, 1.0f));
                    _displacement[x + y * _textureSize] = new Vector3(r, g, b);
                }
            }

            // Apply changes to texture.
            //_displacementTexture.Apply();
        }
    }

    // Returns a displacement vector based on a given world position.
    public Vector3 GetDisplacement(Vector3 position)
    {
        // Calculate relative pixel coordinates of sample.
        float x = _textureSize * (position.x / _patchSize + 0.5f) % _textureSize;
        float y = _textureSize * (position.z / _patchSize + 0.5f) % _textureSize;

        // Calculate minimum and maximum values for sample.
        int minX = (int)x;
        int minY = (int)y;

        int maxX = (minX + 1) % _textureSize;
        int maxY = (minY + 1) % _textureSize;

        // Extract decimal value from sample coordinates, used for interpolation.
        float lerpX = x - minX;
        float lerpY = y - minY;

        // Interpolate between minimum and maximum values of x.
        Vector3 top = Vector3.Lerp(_displacement[minX + minY * _textureSize], _displacement[maxX + minY * _textureSize], lerpX);
        Vector3 bottom = Vector3.Lerp(_displacement[minX + maxY * _textureSize], _displacement[maxX + maxY * _textureSize], lerpX);

        // Interpolate between top and bottom samples.
        Vector3 sample = Vector3.Lerp(top, bottom, lerpY);

        // Return sample as a vector.
        return new Vector3(sample.x, sample.y, sample.z);
    }

    public void SetGaussianNoise(Complex[] noiseTexture)
    {
        _gaussianNoise = noiseTexture;
    }
}