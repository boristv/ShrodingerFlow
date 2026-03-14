using UnityEngine;

namespace ComputeShaderSF
{
    public class CSFFT
    {
        private ComputeShader _shader;
        private int _bitReverseKernel;
        private int _butterflyKernel;
        private int _resX, _resY, _resZ;

        public void Init(ComputeShader shader, int resX, int resY, int resZ)
        {
            _shader = shader;
            _resX = resX;
            _resY = resY;
            _resZ = resZ;
            _bitReverseKernel = shader.FindKernel("BitReverse");
            _butterflyKernel = shader.FindKernel("Butterfly");
        }

        public void FFT3D(ComputeBuffer data, bool inverse)
        {
            FFT1DAlongAxis(data, _resZ, 1,
                _resX * _resY, 1, _resZ, inverse);

            FFT1DAlongAxis(data, _resY, _resZ,
                _resX * _resZ, _resZ, _resY * _resZ, inverse);

            FFT1DAlongAxis(data, _resX, _resY * _resZ,
                _resY * _resZ, _resY * _resZ, _resX * _resY * _resZ, inverse);
        }

        private static int IntLog2(int n)
        {
            int result = 0;
            while (n > 1) { n >>= 1; result++; }
            return result;
        }

        private void FFT1DAlongAxis(ComputeBuffer data, int length, int stride,
            int count, int innerSize, int outerSize, bool inverse)
        {
            int log2N = IntLog2(length);

            _shader.SetInt("_FFTLength", length);
            _shader.SetInt("_FFTStride", stride);
            _shader.SetInt("_FFTCount", count);
            _shader.SetInt("_FFTInnerSize", innerSize);
            _shader.SetInt("_FFTOuterSize", outerSize);
            _shader.SetInt("_FFTLog2N", log2N);
            _shader.SetFloat("_FFTDirection", inverse ? 1.0f : -1.0f);

            _shader.SetBuffer(_bitReverseKernel, "_Data", data);
            int bitRevGroups = (length * count + 255) / 256;
            _shader.Dispatch(_bitReverseKernel, bitRevGroups, 1, 1);

            _shader.SetBuffer(_butterflyKernel, "_Data", data);
            int butterflyGroups = ((length / 2) * count + 255) / 256;
            for (int s = 0; s < log2N; s++)
            {
                _shader.SetInt("_FFTStage", s);
                _shader.Dispatch(_butterflyKernel, butterflyGroups, 1, 1);
            }
        }
    }
}
