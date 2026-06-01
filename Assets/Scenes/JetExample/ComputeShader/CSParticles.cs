using System;
using UnityEngine;

namespace ComputeShaderSF
{
    public class CSParticles : IDisposable
    {
        private ComputeShader _shader;
        private int _addKernel, _interpKernel, _rk4Kernel, _wrapKernel;
        private int _clampKernel, _pushSolidKernel;
        private int _disperseKernel;
        private int _sourceDriftKernel;
        private int _ventDriftKernel;
        private int _deflectWallKernel;

        private ComputeBuffer _x, _y, _z;
        private ComputeBuffer _k1x, _k1y, _k1z;
        private ComputeBuffer _k2x, _k2y, _k2z;
        private ComputeBuffer _k3x, _k3y, _k3z;
        private ComputeBuffer _k4x, _k4y, _k4z;

        private ComputeBuffer _stagingX, _stagingY, _stagingZ;
        private int _stagingCapacity;

        private int _size;
        private int _writeHead;
        private int _maxCnt;
        private float _dt;
        private int _torResX, _torResY, _torResZ;
        private float _torDX, _torDY, _torDZ;

        public void Init(ComputeShader shader, int maxParticles, CSISF isf)
        {
            _shader = shader;
            _maxCnt = maxParticles;
            _size = 0;

            _addKernel = shader.FindKernel("AddParticles");
            _interpKernel = shader.FindKernel("InterpolateVelocity");
            _rk4Kernel = shader.FindKernel("RK4Update");
            _wrapKernel = shader.FindKernel("WrapPositions");
            _clampKernel = -1;
            if (shader.HasKernel("ClampPositionsToVolume"))
                _clampKernel = shader.FindKernel("ClampPositionsToVolume");
            _pushSolidKernel = -1;
            if (shader.HasKernel("PushParticlesOutOfSolidMask"))
                _pushSolidKernel = shader.FindKernel("PushParticlesOutOfSolidMask");
            _disperseKernel = -1;
            if (shader.HasKernel("DisperseMazeParticles"))
                _disperseKernel = shader.FindKernel("DisperseMazeParticles");
            _sourceDriftKernel = -1;
            if (shader.HasKernel("DriftMazeSourceParticles"))
                _sourceDriftKernel = shader.FindKernel("DriftMazeSourceParticles");
            _ventDriftKernel = -1;
            if (shader.HasKernel("DriftMazeTowardVent"))
                _ventDriftKernel = shader.FindKernel("DriftMazeTowardVent");
            _deflectWallKernel = -1;
            if (shader.HasKernel("DeflectMazeAtWallParticles"))
                _deflectWallKernel = shader.FindKernel("DeflectMazeAtWallParticles");

            _x = new ComputeBuffer(maxParticles, sizeof(float));
            _y = new ComputeBuffer(maxParticles, sizeof(float));
            _z = new ComputeBuffer(maxParticles, sizeof(float));

            _k1x = new ComputeBuffer(maxParticles, sizeof(float));
            _k1y = new ComputeBuffer(maxParticles, sizeof(float));
            _k1z = new ComputeBuffer(maxParticles, sizeof(float));
            _k2x = new ComputeBuffer(maxParticles, sizeof(float));
            _k2y = new ComputeBuffer(maxParticles, sizeof(float));
            _k2z = new ComputeBuffer(maxParticles, sizeof(float));
            _k3x = new ComputeBuffer(maxParticles, sizeof(float));
            _k3y = new ComputeBuffer(maxParticles, sizeof(float));
            _k3z = new ComputeBuffer(maxParticles, sizeof(float));
            _k4x = new ComputeBuffer(maxParticles, sizeof(float));
            _k4y = new ComputeBuffer(maxParticles, sizeof(float));
            _k4z = new ComputeBuffer(maxParticles, sizeof(float));

            _dt = isf.dt;
            _torResX = isf.resX;
            _torResY = isf.resY;
            _torResZ = isf.resZ;
            _torDX = isf.dx;
            _torDY = isf.dy;
            _torDZ = isf.dz;
        }

        public int Size => _size;

        public void AddParticles(float[] xx, float[] yy, float[] zz, int count,
            bool ring = false)
        {
            if (count == 0) return;

            if (_writeHead + count > _maxCnt)
            {
                if (!ring) return;
                _writeHead = 0;
            }

            EnsureStagingCapacity(count);
            _stagingX.SetData(xx, 0, 0, count);
            _stagingY.SetData(yy, 0, 0, count);
            _stagingZ.SetData(zz, 0, 0, count);

            _shader.SetInt("_ParticleCount", count);
            _shader.SetInt("_ParticleOffset", _writeHead);

            _shader.SetBuffer(_addKernel, "_PosX", _x);
            _shader.SetBuffer(_addKernel, "_PosY", _y);
            _shader.SetBuffer(_addKernel, "_PosZ", _z);
            _shader.SetBuffer(_addKernel, "_NewX", _stagingX);
            _shader.SetBuffer(_addKernel, "_NewY", _stagingY);
            _shader.SetBuffer(_addKernel, "_NewZ", _stagingZ);
            _shader.Dispatch(_addKernel, (count + 255) / 256, 1, 1);

            _writeHead += count;
            if (_size < _writeHead)
                _size = Mathf.Min(_writeHead, _maxCnt);
        }

        private void EnsureStagingCapacity(int needed)
        {
            if (_stagingCapacity >= needed) return;

            _stagingX?.Release();
            _stagingY?.Release();
            _stagingZ?.Release();

            _stagingCapacity = Mathf.Max(needed, 256);
            _stagingX = new ComputeBuffer(_stagingCapacity, sizeof(float));
            _stagingY = new ComputeBuffer(_stagingCapacity, sizeof(float));
            _stagingZ = new ComputeBuffer(_stagingCapacity, sizeof(float));
        }

        public void CalculateMovement(CSVelocity vel)
        {
            CalculateMovement(vel, false);
        }

        public void CalculateMovement(CSVelocity vel, bool clampSampling)
        {
            if (_size == 0) return;

            SetTorusUniforms();
            _shader.SetInt("_ClampSampling", clampSampling ? 1 : 0);
            _shader.SetInt("_ParticleCount", _size);
            int groups = (_size + 255) / 256;

            RunRK4Step(vel, _k1x, _k1y, _k1z, _k1x, _k1y, _k1z, 0f, groups);
            RunRK4Step(vel, _k1x, _k1y, _k1z, _k2x, _k2y, _k2z, _dt * 0.5f, groups);
            RunRK4Step(vel, _k2x, _k2y, _k2z, _k3x, _k3y, _k3z, _dt * 0.5f, groups);
            RunRK4Step(vel, _k3x, _k3y, _k3z, _k4x, _k4y, _k4z, _dt, groups);

            _shader.SetFloat("_DT", _dt);
            _shader.SetBuffer(_rk4Kernel, "_PosX", _x);
            _shader.SetBuffer(_rk4Kernel, "_PosY", _y);
            _shader.SetBuffer(_rk4Kernel, "_PosZ", _z);
            _shader.SetBuffer(_rk4Kernel, "_K1X", _k1x);
            _shader.SetBuffer(_rk4Kernel, "_K1Y", _k1y);
            _shader.SetBuffer(_rk4Kernel, "_K1Z", _k1z);
            _shader.SetBuffer(_rk4Kernel, "_K2X", _k2x);
            _shader.SetBuffer(_rk4Kernel, "_K2Y", _k2y);
            _shader.SetBuffer(_rk4Kernel, "_K2Z", _k2z);
            _shader.SetBuffer(_rk4Kernel, "_K3X", _k3x);
            _shader.SetBuffer(_rk4Kernel, "_K3Y", _k3y);
            _shader.SetBuffer(_rk4Kernel, "_K3Z", _k3z);
            _shader.SetBuffer(_rk4Kernel, "_K4X", _k4x);
            _shader.SetBuffer(_rk4Kernel, "_K4Y", _k4y);
            _shader.SetBuffer(_rk4Kernel, "_K4Z", _k4z);
            _shader.Dispatch(_rk4Kernel, groups, 1, 1);
        }

        private void RunRK4Step(CSVelocity vel,
            ComputeBuffer shiftX, ComputeBuffer shiftY, ComputeBuffer shiftZ,
            ComputeBuffer outX, ComputeBuffer outY, ComputeBuffer outZ,
            float shiftFactor, int groups)
        {
            _shader.SetFloat("_ShiftFactor", shiftFactor);

            _shader.SetBuffer(_interpKernel, "_PosX", _x);
            _shader.SetBuffer(_interpKernel, "_PosY", _y);
            _shader.SetBuffer(_interpKernel, "_PosZ", _z);
            _shader.SetBuffer(_interpKernel, "_ShiftX", shiftX);
            _shader.SetBuffer(_interpKernel, "_ShiftY", shiftY);
            _shader.SetBuffer(_interpKernel, "_ShiftZ", shiftZ);
            _shader.SetBuffer(_interpKernel, "_VelFieldX", vel.vx);
            _shader.SetBuffer(_interpKernel, "_VelFieldY", vel.vy);
            _shader.SetBuffer(_interpKernel, "_VelFieldZ", vel.vz);
            _shader.SetBuffer(_interpKernel, "_OutX", outX);
            _shader.SetBuffer(_interpKernel, "_OutY", outY);
            _shader.SetBuffer(_interpKernel, "_OutZ", outZ);
            _shader.Dispatch(_interpKernel, groups, 1, 1);
        }

        private void SetTorusUniforms()
        {
            _shader.SetInt("_TorResX", _torResX);
            _shader.SetInt("_TorResY", _torResY);
            _shader.SetInt("_TorResZ", _torResZ);
            _shader.SetFloat("_TorDX", _torDX);
            _shader.SetFloat("_TorDY", _torDY);
            _shader.SetFloat("_TorDZ", _torDZ);
        }

        public void WrapPositions(float volSizeX, float volSizeY, float volSizeZ)
        {
            if (_size == 0) return;
            _shader.SetInt("_ParticleCount", _size);
            _shader.SetFloat("_VolSizeX", volSizeX);
            _shader.SetFloat("_VolSizeY", volSizeY);
            _shader.SetFloat("_VolSizeZ", volSizeZ);
            _shader.SetBuffer(_wrapKernel, "_PosX", _x);
            _shader.SetBuffer(_wrapKernel, "_PosY", _y);
            _shader.SetBuffer(_wrapKernel, "_PosZ", _z);
            _shader.Dispatch(_wrapKernel, (_size + 255) / 256, 1, 1);
        }

        public void ClampPositionsToVolume(float volSizeX, float volSizeY, float volSizeZ)
        {
            if (_size == 0 || _clampKernel < 0) return;
            _shader.SetInt("_ParticleCount", _size);
            _shader.SetFloat("_VolSizeX", volSizeX);
            _shader.SetFloat("_VolSizeY", volSizeY);
            _shader.SetFloat("_VolSizeZ", volSizeZ);
            _shader.SetBuffer(_clampKernel, "_PosX", _x);
            _shader.SetBuffer(_clampKernel, "_PosY", _y);
            _shader.SetBuffer(_clampKernel, "_PosZ", _z);
            _shader.Dispatch(_clampKernel, (_size + 255) / 256, 1, 1);
        }

        public void PushOutOfSolidMask(ComputeBuffer solidMask, int searchCells, Vector3 ventTarget, float ventBiasWeight)
        {
            if (_size == 0 || solidMask == null || _pushSolidKernel < 0) return;
            SetTorusUniforms();
            _shader.SetInt("_ParticleCount", _size);
            _shader.SetInt("_PushSearchCells", Mathf.Max(1, searchCells));
            _shader.SetVector("_VentTarget", ventTarget);
            _shader.SetFloat("_VentBiasWeight", ventBiasWeight);
            _shader.SetBuffer(_pushSolidKernel, "_PosX", _x);
            _shader.SetBuffer(_pushSolidKernel, "_PosY", _y);
            _shader.SetBuffer(_pushSolidKernel, "_PosZ", _z);
            _shader.SetBuffer(_pushSolidKernel, "_SolidMask", solidMask);
            _shader.Dispatch(_pushSolidKernel, (_size + 255) / 256, 1, 1);
        }

        /// <summary>Случайный шаг в плоскости XZ — турбулентное рассеивание; у стен усилено через маску.</summary>
        public void DisperseMaze(ComputeBuffer solidMask, float sigma, float wallBoost, int seed)
        {
            if (_size == 0 || _disperseKernel < 0 || sigma <= 0f) return;
            SetTorusUniforms();
            _shader.SetInt("_ParticleCount", _size);
            _shader.SetFloat("_DisperseSigma", sigma);
            _shader.SetFloat("_DisperseWallBoost", wallBoost);
            _shader.SetInt("_DisperseSeed", seed);
            if (solidMask != null)
                _shader.SetBuffer(_disperseKernel, "_SolidMask", solidMask);
            _shader.SetBuffer(_disperseKernel, "_PosX", _x);
            _shader.SetBuffer(_disperseKernel, "_PosY", _y);
            _shader.SetBuffer(_disperseKernel, "_PosZ", _z);
            _shader.Dispatch(_disperseKernel, (_size + 255) / 256, 1, 1);
        }

        public void DriftInSourceBox(Vector3 center, Vector3 half, Vector3 velocity, float stepDt)
        {
            if (_size == 0 || _sourceDriftKernel < 0 || velocity.sqrMagnitude < 1e-12f) return;
            SetTorusUniforms();
            _shader.SetInt("_ParticleCount", _size);
            _shader.SetVector("_SourceCenter", center);
            _shader.SetVector("_SourceHalf", half);
            _shader.SetVector("_SourceDriftV", velocity);
            _shader.SetFloat("_SourceDriftDT", stepDt);
            _shader.SetBuffer(_sourceDriftKernel, "_PosX", _x);
            _shader.SetBuffer(_sourceDriftKernel, "_PosY", _y);
            _shader.SetBuffer(_sourceDriftKernel, "_PosZ", _z);
            _shader.Dispatch(_sourceDriftKernel, (_size + 255) / 256, 1, 1);
        }

        public void DriftTowardVent(Vector3 ventCenter, float strength, float stepDt)
        {
            if (_size == 0 || _ventDriftKernel < 0 || strength <= 0f) return;
            SetTorusUniforms();
            _shader.SetInt("_ParticleCount", _size);
            _shader.SetVector("_VentDriftTarget", ventCenter);
            _shader.SetFloat("_VentDriftStrength", strength);
            _shader.SetFloat("_VentDriftDT", stepDt);
            _shader.SetBuffer(_ventDriftKernel, "_PosX", _x);
            _shader.SetBuffer(_ventDriftKernel, "_PosY", _y);
            _shader.SetBuffer(_ventDriftKernel, "_PosZ", _z);
            _shader.Dispatch(_ventDriftKernel, (_size + 255) / 256, 1, 1);
        }

        public void DeflectAtWall(ComputeBuffer solidMask, float stepWorld, int searchCells)
        {
            if (_size == 0 || solidMask == null || _deflectWallKernel < 0 || stepWorld <= 0f) return;
            SetTorusUniforms();
            _shader.SetInt("_ParticleCount", _size);
            _shader.SetFloat("_DeflectStep", stepWorld);
            _shader.SetInt("_DeflectSearchCells", Mathf.Max(4, searchCells));
            _shader.SetBuffer(_deflectWallKernel, "_SolidMask", solidMask);
            _shader.SetBuffer(_deflectWallKernel, "_PosX", _x);
            _shader.SetBuffer(_deflectWallKernel, "_PosY", _y);
            _shader.SetBuffer(_deflectWallKernel, "_PosZ", _z);
            _shader.Dispatch(_deflectWallKernel, (_size + 255) / 256, 1, 1);
        }

        public void ReadPositions(float[] outX, float[] outY, float[] outZ)
        {
            if (_size == 0) return;
            _x.GetData(outX, 0, 0, _size);
            _y.GetData(outY, 0, 0, _size);
            _z.GetData(outZ, 0, 0, _size);
        }

        public void WritePositions(float[] px, float[] py, float[] pz, int count)
        {
            if (count <= 0 || count > _size) return;
            _x.SetData(px, 0, 0, count);
            _y.SetData(py, 0, 0, count);
            _z.SetData(pz, 0, 0, count);
        }

        /// <param name="reorderAux1">
        /// Опционально: массив той же логической длины, что и частицы (напр. прошлые мировые позицы для цвета скорости).
        /// При уплотнении переставляется тем же образом, что и позиции — иначе индекс перепутывает историю после компакта.
        /// </param>
        /// <param name="reorderAux2">Опционально второй массив (напр. сглаженная скорость для отображения).</param>
        public int CompactParticles(float[] px, float[] py, float[] pz,
            float maxX, float maxY, float maxZ,
            Vector3[] reorderAux1 = null, Vector3[] reorderAux2 = null)
        {
            if (_size == 0) return 0;

            _x.GetData(px, 0, 0, _size);
            _y.GetData(py, 0, 0, _size);
            _z.GetData(pz, 0, 0, _size);

            int alive = 0;
            for (int i = 0; i < _size; i++)
            {
                if (px[i] < 0f || px[i] > maxX ||
                    py[i] < 0f || py[i] > maxY ||
                    pz[i] < 0f || pz[i] > maxZ)
                    continue;

                if (alive != i)
                {
                    px[alive] = px[i];
                    py[alive] = py[i];
                    pz[alive] = pz[i];
                    if (reorderAux1 != null)
                        reorderAux1[alive] = reorderAux1[i];
                    if (reorderAux2 != null)
                        reorderAux2[alive] = reorderAux2[i];
                }
                alive++;
            }

            if (alive < _size)
            {
                _x.SetData(px, 0, 0, alive);
                _y.SetData(py, 0, 0, alive);
                _z.SetData(pz, 0, 0, alive);
                _size = alive;
                _writeHead = alive;
            }

            return alive;
        }

        public void Dispose()
        {
            _x?.Release(); _y?.Release(); _z?.Release();
            _k1x?.Release(); _k1y?.Release(); _k1z?.Release();
            _k2x?.Release(); _k2y?.Release(); _k2z?.Release();
            _k3x?.Release(); _k3y?.Release(); _k3z?.Release();
            _k4x?.Release(); _k4y?.Release(); _k4z?.Release();
            _stagingX?.Release(); _stagingY?.Release(); _stagingZ?.Release();
        }
    }
}
