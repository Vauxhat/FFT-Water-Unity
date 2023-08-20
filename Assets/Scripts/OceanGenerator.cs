using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Vauxhat;

using Complex = System.Numerics.Complex;

public class OceanGenerator : MonoBehaviour
{
    [System.Serializable]
    public class WaveData
    {
        [Min(0.0f)] public int patchSize = 128;
        [Min(0.0f)] public float repeatTime = 200.0f;
        [Min(0.0f)] public float minWavelength = 0.001f;
        [Min(0.0f)] public float phillipsConstant = 0.0002f;
    }

    public WaveData[] _waveData = new WaveData[1];

    // Public variables for handling wave simulation.
    [Range(0.0f, 360.0f)] public float _windAngle = 0.0f;
    [Min(0.0f)] public float _windSpeed = 6.0f;
    [Min(0.0f)] public float _gravity = 9.81f;
    [Min(0.0f)] public float _steepness = 0.0f;

    private int _textureSizeGPU = 128;
    private int _textureSizeCPU = 128;

    private FourierWaterCPU _fourierWaterCPU;
    private FourierWaterGPU _fourierWaterGPU;

    // Compute shaders.
    public ComputeShader _stationarySpectrumShader;
    public ComputeShader _fourierTransformShader;
    public ComputeShader _timeSpectrumShader;
    public ComputeShader _dispersionShader;
    public ComputeShader _textureMergeShader;

    private Material _material;

    public bool _renderToggle = false;

    // Start is called before the first frame update
    void Start()
    {
        _fourierWaterGPU = new FourierWaterGPU();
        _fourierWaterCPU = new FourierWaterCPU();

        // Generate gaussian noise texture.
        InitialiseGaussianNoise();

        // Initialise GPU water simulation.
        _fourierWaterGPU._windAngle = _windAngle;
        _fourierWaterGPU._windSpeed = _windSpeed;
        _fourierWaterGPU._gravity = _gravity;
        _fourierWaterGPU._steepness = _steepness;
        _fourierWaterGPU._patchSize = _waveData[0].patchSize;
        _fourierWaterGPU._repeatTime = _waveData[0].repeatTime;
        _fourierWaterGPU._minWavelength = _waveData[0].minWavelength;
        _fourierWaterGPU._phillipsConstant = _waveData[0].phillipsConstant;
        _fourierWaterGPU._textureSize = _textureSizeGPU;

        _fourierWaterGPU._stationarySpectrumShader = _stationarySpectrumShader;
        _fourierWaterGPU._fourierTransformShader = _fourierTransformShader;
        _fourierWaterGPU._timeSpectrumShader = _timeSpectrumShader;
        _fourierWaterGPU._dispersionShader = _dispersionShader;
        _fourierWaterGPU._textureMergeShader = _textureMergeShader;

        _fourierWaterGPU.Initialise();

        // Initialise CPU water simulation.
        _fourierWaterCPU._windAngle = _windAngle;
        _fourierWaterCPU._windSpeed = _windSpeed;
        _fourierWaterCPU._gravity = _gravity;
        _fourierWaterCPU._steepness = _steepness;
        _fourierWaterCPU._patchSize = _waveData[0].patchSize;
        _fourierWaterCPU._repeatTime = _waveData[0].repeatTime;
        _fourierWaterCPU._minWavelength = _waveData[0].minWavelength;
        _fourierWaterCPU._phillipsConstant = _waveData[0].phillipsConstant;
        _fourierWaterCPU._textureSize = _textureSizeCPU;

        _fourierWaterCPU.Initialise();

        // Store local reference to object material.
        _material = this.GetComponent<Renderer>().sharedMaterial;
    }

    // Update is called once per frame
    void Update()
    {
        float time = Time.time;

        _fourierWaterGPU.Update(time);
        _fourierWaterCPU.Update(time);

        // Check if material exists.
        if (_material)
        {
            if (_renderToggle)
            {
                // Pass heightfield texture to material shader.
                _material.SetTexture("_MainTex", _fourierWaterCPU._displacementTexture);
            }
            else
            {
                // Pass heightfield texture to material shader.
                _material.SetTexture("_MainTex", _fourierWaterGPU._displacement);
            }
        }
    }

    void InitialiseGaussianNoise()
    {
        // Initialise array of gaussian values, generating two per pixel.
        Complex[] gaussianNoiseGPU = new Complex[_textureSizeGPU * _textureSizeGPU * 2];
        Texture2D gaussianNoiseTexture = new Texture2D(_textureSizeGPU, _textureSizeGPU, TextureFormat.RGBAFloat, false);

        // Fill array with independent gaussian values.
        for (int i = 0; i < _textureSizeGPU * _textureSizeGPU; i++)
        {
            int x = i % _textureSizeGPU;
            int y = i / _textureSizeGPU;

            float r, g, b, a;

            // Generate two pairs of gaussian noise.
            MathsExt.GaussianRandom(out r, out g);
            MathsExt.GaussianRandom(out b, out a);

            // Store gaussian values in array.
            gaussianNoiseGPU[2 * i] = new Complex(r, g);
            gaussianNoiseGPU[2 * i + 1] = new Complex(b, a);

            // Store gaussian variables in texture.
            gaussianNoiseTexture.SetPixel(x, y, new Color(r, g, b, a));
        }

        gaussianNoiseTexture.Apply();

        // Send noise to GPU.
        _fourierWaterGPU.SetGaussianNoise(gaussianNoiseTexture);

        // Create seperate array for CPU values.
        Complex[] gaussianNoiseCPU = new Complex[_textureSizeCPU * _textureSizeCPU * 2];

        // Calculate scale from CPU array to GPU array.
        int scale = _textureSizeGPU / _textureSizeCPU;

        // Loop through each pixel in the texture.
        for (int y = 0; y < _textureSizeCPU; y++)
        {
            for (int x = 0; x < _textureSizeCPU; x++)
            {
                // Calculate index in CPU array and relative position in GPU array.
                int indexCPU = x + y * _textureSizeCPU;
                //int indexGPU = scale * (x + y * _textureSizeCPU);

                int scaledX = scale * x;
                int scaledY = scale * y;

                int indexGPU = scaledX + scaledY * _textureSizeCPU;

                // Copy data from GPU texture over to CPU.
                gaussianNoiseCPU[2 * indexCPU] = gaussianNoiseGPU[2 * indexGPU];
                gaussianNoiseCPU[2 * indexCPU + 1] = gaussianNoiseGPU[2 * indexGPU + 1];
            }
        }

        // Send noise to CPU.
        _fourierWaterCPU.SetGaussianNoise(gaussianNoiseCPU);
    }

    public Vector3 GetDisplacement(Vector3 position)
    {
        return _fourierWaterCPU.GetDisplacement(position);
    }
}
