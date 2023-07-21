using System.Collections;
using System.Collections.Generic;
using Unity.Jobs;
using Unity.Collections;
using UnityEngine;
using Vauxhat;

using Complex = System.Numerics.Complex;

public struct FastFourierTransform
{
    private struct FFTJob : IJobParallelFor
    {
        // Input buffer.
        [NativeDisableParallelForRestriction] public NativeArray<Complex> _input;

        // Read only input variables.
        [ReadOnly] public NativeArray<int> _reversedIndex;
        [ReadOnly] public int _textureSize;
        [ReadOnly] public int _fourierStages;
        [ReadOnly] public bool _direction;

        public void Execute(int index)
        {
            // Calculate FFT based on direction.
            if (_direction)
            {
                // Vertical Pass.
                FastFourierTransform(_textureSize, index);
            }
            else
            {
                // Horizontal Pass.
                FastFourierTransform(1, index * _textureSize);
            }
        }

        private void FastFourierTransform(int stride, int offset)
        {
            // Create a pingpong buffer of size N.
            Complex[,] buffer = new Complex[2, _textureSize];

            // Copy current row of input into buffer using bit-reversal permutation.
            for (int i = 0; i < _textureSize; i++)
            {
                buffer[0, i] = _input[_reversedIndex[i] * stride + offset];
            }

            // Initialise pingpong state.
            int pingpong = 0;

            // Perform FFT for each stage (zero to Log N).
            for (int s = 0; s < _fourierStages; s++)
            {
                // Calculate 2 to the power of s.
                int m = 2 << s;

                float twiddle = (-2.0f * Mathf.PI) / m;
                Complex wm = new Complex(Mathf.Cos(twiddle), Mathf.Sin(twiddle));

                // Loop through each index in the line, in groups of size m.
                for (int k = 0; k < _textureSize; k += m)
                {
                    Complex w = new Complex(1.0f, 0.0f);

                    // loop through the first half of indices in the group.
                    for (int j = 0; j < (m / 2); j++)
                    {
                        // Get top and bottom index from sample.
                        int topIndex = k + j;
                        int bottomIndex = k + j + (m / 2);

                        // Get complex variables from previous stage.
                        Complex a = buffer[pingpong, topIndex];
                        Complex b = buffer[pingpong, bottomIndex] * w;

                        // Calculate new values for top and bottom index.
                        buffer[(pingpong + 1) % 2, topIndex] = a + b;
                        buffer[(pingpong + 1) % 2, bottomIndex] = a - b;

                        w *= wm;
                    }
                }

                // Switch to next array in buffer.
                pingpong = (pingpong + 1) % 2;
            }

            // Copy final FFT for the current sequence into the output buffer.
            for (int i = 0; i < _textureSize; i++)
            {
                _input[i * stride + offset] = buffer[pingpong, i];
            }
        }
    }

    // Calculates an array of reversed indices in the range 0 - textureSize.
    public static int[] PrecomputeReversedIndex(int textureSize)
    {
        // Initialise output array.
        int[] output = new int[textureSize];

        // Calculate the number of bits used to represent indices.
        int bits = (int)Mathf.Log(textureSize, 2);

        // Loop through the first half of the array.
        for (int i = 0; i < output.Length / 2; i++)
        {
            // Calculate reversed index.
            int reversed = MathsExt.BitReverse(i, bits);

            // Update values at current and reversed indices.
            output[i] = reversed;
            output[reversed] = i;
        }

        return output;
    }

    // Computes the 2D FFT of the input, intermediate variables will be calculated at run-time.
    public static void Compute(Complex[] input, int textureSize)
    {
        // Create array for reversed index.
        int[] reversedIndex = new int[textureSize];

        // Calculate number of stages.
        int stages = (int)Mathf.Log(textureSize, 2);

        // Initialise array.
        for (int i = 0; i < textureSize; i++)
        {
            reversedIndex[i] = MathsExt.BitReverse(i, stages);
        }

        // Compute FFT.
        Compute(input, textureSize, reversedIndex);
    }

    // Computes the 2D FFT of the input, intermediate variables are passed into the function.
    public static void Compute(Complex[] input, int textureSize, int[] reversedIndex)
    {
        // Create native array for input data.
        NativeArray<Complex> nativeInput = new NativeArray<Complex>(input.Length, Allocator.TempJob);
        NativeArray<int> nativeReversedIndex = new NativeArray<int>(textureSize, Allocator.TempJob);

        // Calculate number of stages.
        int stages = (int)Mathf.Log(textureSize, 2);

        // Copy reverse index to native array.
        for (int i = 0; i < textureSize; i++)
        {
            nativeReversedIndex[i] = reversedIndex[i];
        }

        // Copy source array into native array.
        for (int i = 0; i < input.Length; i++)
        {
            nativeInput[i] = input[i];
        }

        // Create new job.
        FFTJob job = new FFTJob()
        {
            _textureSize = textureSize,
            _fourierStages = stages,
            _input = nativeInput,
            _reversedIndex = nativeReversedIndex
        };

        // Create job handle.
        JobHandle jobHandle;

        // Perform horizontal operation.
        job._direction = false;
        jobHandle = job.Schedule(textureSize, 1);

        // Wait until all threads are finished.
        jobHandle.Complete();

        // Perform vertical operation.
        job._direction = true;
        jobHandle = job.Schedule(textureSize, 1);

        // Wait until all threads are finished.
        jobHandle.Complete();

        // Copy output to source array.
        for (int i = 0; i < input.Length; i++)
        {
            input[i] = nativeInput[i];
        }

        // Dispose of native array structure.
        nativeInput.Dispose();
        nativeReversedIndex.Dispose();
    }
}
