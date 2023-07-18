using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class WaterGenerator : MonoBehaviour
{
    [System.Serializable]
    public class Wave
    {
        [Range(0.0f, 1.0f)] public float steepness = 0.0f;
        public float amplitude = 1.0f;
        [Range(0.0f, 360.0f)] public float angle = 0.0f;
        public float speed = 1.0f;
        [Min(0.001f)] public float wavelength = 1.0f;

        public float frequency = 1.0f;
    }

    public Wave[] waves = new Wave[128];
    public int numWaves = 128;

    [Min(0.001f)] public float medianWavelength = 1.0f;
    [Min(0.001f)] public float medianAmplitude = 1.0f;
    [Range(0.0f, 360.0f)] public float windDirection = 0.0f;
    [Range(0.0f, 360.0f)] public float windCone = 30.0f;

    // Start is called before the first frame update
    void Start()
    {
        //GenerateRandomWaves();

        //Wave[] test = new Wave[128];

        float maxHeight = 0;

        for (int i = 0; i < waves.Length; i++)
        {
            waves[i].steepness = Random.Range(1.0f, 1.0f);
            waves[i].amplitude = Random.Range(0.01f, 0.05f);
            waves[i].angle = Random.Range(60.0f, 120.0f);
            waves[i].speed = Random.Range(1.0f, 8.0f);
            waves[i].wavelength = Random.Range(0.4f, 64.0f);

            //waves[i].frequency = Mathf.Sqrt(gravity * (2.0f * Mathf.PI) / waves[i].wavelength);
            waves[i].frequency = (2.0f * Mathf.PI) / waves[i].wavelength;

            // Sanitise angle value between 0 and 360.
            if (waves[i].angle < 0.0f) waves[i].angle += 360.0f;
            else waves[i].angle = waves[i].angle % 360.0f;

            // Add current amplitude to max height.
            maxHeight += waves[i].amplitude;
        }

        // Get material from renderer.
        Material material = this.gameObject.GetComponent<Renderer>().material;

        // Check if material was found.
        if (material)
        {
            material.SetFloat("_GradientHeight", maxHeight);
        }

        /*for (int i = 0; i < waves.Length; i++)
        {
            float wavelength = Random.Range(medianWavelength * 0.5f, medianWavelength * 2.0f);
            float amplitude = Random.Range(medianAmplitude * 0.5f, medianAmplitude * 2.0f);

            float angle = Random.Range(windDirection - windCone * 0.5f, windDirection + windCone * 0.5f);
            
            if (angle < 0.0f) angle += 360.0f;
            else angle = angle % 360.0f;

            const float gravity = 9.8f;

            float frequency = Mathf.Sqrt(gravity * (2.0f * Mathf.PI) / wavelength);
        }*/
    }

    // Update is called once per frame
    void Update()
    {
        // Get material from renderer.
        Material material = this.gameObject.GetComponent<Renderer>().material;

        if (material)
        {
            float[] steepness = new float[128];
            float[] amplitude = new float[128];
            float[] x = new float[128];
            float[] y = new float[128];
            float[] speed = new float[128];
            float[] wavelength = new float[128];
            float[] frequency = new float[128];

            for (int i = 0; i < waves.Length; i++)
            {
                steepness[i] = waves[i].steepness;
                amplitude[i] = waves[i].amplitude;
                x[i] = Mathf.Sin(Mathf.Deg2Rad * waves[i].angle);
                y[i] = Mathf.Cos(Mathf.Deg2Rad * waves[i].angle);
                speed[i] = waves[i].speed;
                wavelength[i] = waves[i].wavelength;
                frequency[i] = waves[i].frequency;
            }

            material.SetInt("_NumWaves", waves.Length);

            material.SetFloatArray("_Steepness", steepness);
            material.SetFloatArray("_Amplitude", amplitude);
            material.SetFloatArray("_DirX", x);
            material.SetFloatArray("_DirY", y);
            material.SetFloatArray("_Speed", speed);
            material.SetFloatArray("_Wavelength", wavelength);
            material.SetFloatArray("_Frequency", frequency);

            // Pass wave variables to shader.
            //material.SetFloat("_Steepness", waves[0].steepness);
            //material.SetFloat("_Amplitude", waves[0].amplitude);
            //material.SetFloat("_DirX", Mathf.Sin(Mathf.Deg2Rad * waves[0].angle));
            //material.SetFloat("_DirY", Mathf.Cos(Mathf.Deg2Rad * waves[0].angle));
            //material.SetFloat("_Speed", waves[0].speed);
            //material.SetFloat("_Wavelength", waves[0].wavelength);
        }
    }

    // Returns a vector offset based on an input position.
    public Vector3 GetWaveOffset(Vector3 position, float cuttoff = 0.0f)
    {
        // Initialise offset.
        Vector3 offset = new Vector3(0.0f, 0.0f, 0.0f);

        foreach (Wave wave in waves)
        {
            // Check if frequency is within cutoff.
            if (wave.wavelength >= cuttoff)
            {
                // Calculate frequency (w) and phase based on input values.
                float w = (2.0f * Mathf.PI) / wave.wavelength;
                float phase = wave.speed * w;
                float WA = w * wave.amplitude;
                float Q = wave.steepness / (w * wave.amplitude * waves.Length);

                // Convert wave angle into normalised direction.
                Vector2 direction = new Vector2(Mathf.Sin(Mathf.Deg2Rad * wave.angle), Mathf.Cos(Mathf.Deg2Rad * wave.angle));

                // Precalculate sine and cosine values.
                float cosine = Mathf.Cos(w * Vector2.Dot(direction, new Vector2(position.x, position.z)) + phase * Time.time);
                float sine = Mathf.Sin(w * Vector2.Dot(direction, new Vector2(position.x, position.z)) + phase * Time.time);

                // Apply current vector to total offset.
                offset.x += Q * wave.amplitude * direction.x * cosine;
                offset.y += wave.amplitude * sine;
                offset.z += Q * wave.amplitude * direction.y * cosine;
            }
            else
            {
                continue;
            }
        }

        // Return final value.
        return offset;
    }

    // Returns a vector offset based on an input position.
    public Vector3 GetWaveNormal(Vector3 position, float cutoff = 0.0f)
    {
        // Initialise normal.
        Vector3 normal = new Vector3(0.0f, 0.0f, 0.0f);

        foreach (Wave wave in waves)
        {
            // Skip wave if wavelength is less than cutoff.
            if (wave.wavelength < cutoff) continue;

            // Calculate frequency (w) and phase based on input values.
            float w = (2.0f * Mathf.PI) / wave.wavelength;
            float phase = wave.speed * w;
            float WA = w * wave.amplitude;
            float Q = wave.steepness / (w * wave.amplitude * waves.Length);

            // Convert wave angle into normalised direction.
            Vector2 direction = new Vector2(Mathf.Sin(Mathf.Deg2Rad * wave.angle), Mathf.Cos(Mathf.Deg2Rad * wave.angle));

            // Precalculate sine and cosine values.
            float cosine = Mathf.Cos(w * Vector2.Dot(direction, new Vector2(position.x, position.z)) + phase * Time.time);
            float sine = Mathf.Sin(w * Vector2.Dot(direction, new Vector2(position.x, position.z)) + phase * Time.time);

            // Calculate surface normal.
            normal.x -= (direction.x * WA * cosine);
            normal.y += Q * WA * sine;
            normal.z -= (direction.y * WA * cosine);
        }

        normal.y = 1.0f - normal.y;
        normal.Normalize();

        // Return final value.
        return normal;
    }

    public void GenerateRandomWaves()
    {
        waves = new Wave[numWaves];

        for(int i = 0; i < numWaves; i++)
        {
            waves[i].steepness = Random.Range(0.0f, 1.0f);
            waves[i].amplitude = Random.Range(0.1f, 0.5f);
            waves[i].angle = Random.Range(85.0f, 95.0f);
            waves[i].speed = Random.Range(12.0f, 30.0f);
            waves[i].wavelength = Random.Range(2.0f, 16.0f);
        }
    }
}
