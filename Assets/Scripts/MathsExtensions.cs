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

        // Returns a random value with normal distribution.
        static public void GaussianRandom(out float a, out float b, float mean = 0, float deviation = 1)
        {
            // Declare local variables.
            float u, v, s;

            // Generate first random value.
            do { u = Random.Range(0.0f, 1.0f); } while (u <= float.Epsilon);

            // Generate second random value.
            v = Random.Range(0.0f, 1.0f);

            // Calculate factor for gaussian distribution.
            s = deviation * Mathf.Sqrt(-2.0f * Mathf.Log(u));

            // Generate two outputs using sine and cosine.
            a = s * Mathf.Cos(2.0f * Mathf.PI * v) + mean;
            b = s * Mathf.Sin(2.0f * Mathf.PI * v) + mean;
        }
    }
}
