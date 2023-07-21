using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using UnityEngine;

namespace Vauxhat
{
    // A collection of extended maths functionality not found in Unity.
    struct MathsExt
    {
        // Returns the bit-reversed integer of a given input.
        static public int BitReverse(int number, int bits = 32)
        {
            // Initialise output.
            int output = 0;

            // Loop thorugh each bit in the number. Unsigned integers are 32 bits in size.
            for (int i = 0; i < bits; i++)
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
    }
}
