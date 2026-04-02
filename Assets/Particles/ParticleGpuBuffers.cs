using UnityEngine;

namespace ShrodingerFlow.Particles
{
    /// <summary>
    /// GPU-буферы позиций/скоростей для отрисовки через ParticleDisplay3D (instanced indirect).
    /// </summary>
    public class ParticleGpuBuffers : MonoBehaviour
    {
        private ComputeBuffer _positions;
        private ComputeBuffer _velocities;
        private ComputeBuffer _debug;

        public ComputeBuffer PositionBuffer => _positions;
        public ComputeBuffer VelocityBuffer => _velocities;
        public ComputeBuffer DebugBuffer => _debug;

        public int ActiveCount { get; private set; }
        public int Capacity { get; private set; }

        private static readonly int StrideFloat3 = sizeof(float) * 3;

        public void EnsureCapacity(int capacity)
        {
            if (capacity <= 0) capacity = 1;
            if (Capacity >= capacity && _positions != null) return;

            ReleaseBuffers();
            Capacity = capacity;

            _positions = new ComputeBuffer(Capacity, StrideFloat3);
            _velocities = new ComputeBuffer(Capacity, StrideFloat3);
            _debug = new ComputeBuffer(Capacity, StrideFloat3);

            var zero = new Vector3[Capacity];
            _debug.SetData(zero);
        }

        /// <summary> Загрузить первые <paramref name="count"/> элементов (массивы не короче count). </summary>
        public void Upload(Vector3[] positions, Vector3[] velocities, int count)
        {
            if (_positions == null || count <= 0)
            {
                ActiveCount = 0;
                return;
            }

            ActiveCount = Mathf.Min(count, Capacity);
            _positions.SetData(positions, 0, 0, ActiveCount);
            _velocities.SetData(velocities, 0, 0, ActiveCount);
        }

        private void ReleaseBuffers()
        {
            _positions?.Release();
            _velocities?.Release();
            _debug?.Release();
            _positions = null;
            _velocities = null;
            _debug = null;
            Capacity = 0;
            ActiveCount = 0;
        }

        private void OnDestroy()
        {
            ReleaseBuffers();
        }
    }
}
