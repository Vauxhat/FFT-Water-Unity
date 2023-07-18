using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using UnityEngine;

namespace Vauxhat
{
    // A collection off extended maths functionality not found in Unity.
    struct MathsExt
    {
        // Returns the bit-reversed integer of a given input.
        static public uint BitReverse(uint number, int bits = 32)
        {
            // Initialise output.
            uint output = 0;

            // Loop thorugh each bit in the number. Unsigned integers are 32 bits in size.
            for (uint i = 0; i < bits; i++)
            {
                // Shift output one to the left.
                output <<= 1;

                // Add first bit from number to the front of the sequence.
                output = output | (number & 1);

                // Move number one to the left, setting up for the next iteration.
                number = number >> 1;
            }

            // Return output.
            return output;
        }

        /// <summary>
        /// Returns a random value with normal (gaussian) distribution using a Box-Muller transform. Mean and standard deviation can be recieved as inputs
        /// </summary>
        /// <param name="mean">The mean or average generated value.</param>
        /// <param name="deviation">The standard deviation of the generated spectrum.</param>
        static public float BoxMullerTransform(float mean = 0, float deviation = 1)
        {
            // Declare local variables.
            float u, v, s;

            // Generate first random value.
            do { u = Random.Range(0.0f, 1.0f); } while (u <= float.Epsilon);

            // Generate second random value.
            v = Random.Range(0.0f, 1.0f);

            // Calculate factor for gaussian distribution.
            s = Mathf.Sqrt(-2.0f * Mathf.Log(u)) * Mathf.Cos(2.0f * Mathf.PI * v);

            // Return output value based on mean, deviation and random distribution.
            return (mean + deviation * s);
        }

        // Returns a random value with normal distribution.
        static public void GaussianRandom(float mean, float deviation, out float a, out float b)
        {
            float u1, u2;

            do
            {
                u1 = Random.Range(0.0f, 1.0f);
            }
            while (u1 <= float.Epsilon);

            u2 = Random.Range(0.0f, 1.0f);

            float mag = deviation * Mathf.Sqrt(-2.0f * Mathf.Log(u1));
            a = mag * Mathf.Cos(2.0f * Mathf.PI * u2) + mean;
            b = mag * Mathf.Sin(2.0f * Mathf.PI * u2) + mean;
        }

        static public float NextGaussianDouble()
        {
            float u, v, S;

            do
            {
                u = 2.0f * Random.value - 1.0f;
                v = 2.0f * Random.value - 1.0f;
                S = u * u + v * v;
            }
            while (S >= 1.0);

            float fac = Mathf.Sqrt(-2.0f * Mathf.Log(S) / S);
            return u * fac;
        }

        static public Complex[] CooleyTukeyFFT(Complex[] a)
        {
            Complex J = new Complex(0, 1);

            int log2n = (int)Mathf.Log(a.Length, 2);

            Complex[] b = new Complex[a.Length];

            int n = 1 << log2n;
            for (uint i = 0; i < n; ++i)
            {
                b[BitReverse(i, log2n)] = a[i];
            }
            for (int s = 1; s <= log2n; ++s)
            {
                int m = 1 << s;
                int m2 = m >> 1;
                Complex w = new Complex(1, 0);
                Complex wm = -J * new Complex(Mathf.Cos(Mathf.PI / m2), Mathf.Sin(Mathf.PI / m2));
                for (int j = 0; j < m2; ++j)
                {
                    for (int k = j; k < n; k += m)
                    {
                        Complex t = w * b[k + m2];
                        Complex u = b[k];
                        b[k] = u + t;
                        b[k + m2] = u - t;
                    }
                    w *= wm;
                }
            }

            return b;
        }

        static public Complex[] FFT(Complex[] a)
        {
            // Store length of array as local variable.
            int N = a.Length;

            // Initialise array for bit reversal.
            Complex[] output = new Complex[N];

            // Calculate the number of bits needed to store the number.
            int bits = (int)(Mathf.Log(N, 2));

            // Copy elements from one array into another, using a bit reversal permutation.
            for (uint i = 0; i < N - 1; i++)
            {
                output[BitReverse(i, bits)] = a[i];
            }

            //int m = 2;

            for (int s = 1; s < bits; s++)
            {
                int m = (int)Mathf.Pow(2, s);

                float factor = (-2.0f * Mathf.PI) / m;
                Complex w_m = new Complex(Mathf.Cos(factor), Mathf.Sin(factor));

                for (int k = 0; k < N - 1; k += m)
                {
                    Complex w = 1;

                    for (int j = 0; j < (m / 2) - 1; j++)
                    {
                        Complex t = w * output[k + j + (m / 2)];
                        Complex u = a[k + j];

                        // Set values in output array.
                        output[k + j] = u + t;
                        output[k + j + (m / 2)] = u - t;

                        w *= w_m;
                    }
                }
            }

            return output;
        }
    }
}
