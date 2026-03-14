using UnityEngine;

namespace ComputeShaderSF
{
    public class CSVelocity : System.IDisposable
    {
        public ComputeBuffer vx;
        public ComputeBuffer vy;
        public ComputeBuffer vz;

        public CSVelocity(int rx, int ry, int rz)
        {
            int count = rx * ry * rz;
            vx = new ComputeBuffer(count, sizeof(float));
            vy = new ComputeBuffer(count, sizeof(float));
            vz = new ComputeBuffer(count, sizeof(float));
        }

        public void Dispose()
        {
            vx?.Release();
            vy?.Release();
            vz?.Release();
        }
    }
}
