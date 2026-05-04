using System;
using UnityEngine;

namespace ComputeShaderSF
{
    public class CSParticles : IDisposable
    {
        private ComputeShader _shader;
        private ComputeShader _chiConstrainShader;
        private int _addKernel, _interpKernel, _rk4Kernel, _wrapKernel, _clampKernel;
        private int _constrainLiquidChiKernel;

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

        /// <param name="chiConstrainShader">Опционально: отдельный compute только с <c>ConstrainParticlesToLiquidChi</c> (не смешивать с основным шейдером частиц).</param>
        public void Init(ComputeShader shader, int maxParticles, CSISF isf,
            ComputeShader chiConstrainShader = null)
        {
            _shader = shader;
            _chiConstrainShader = chiConstrainShader;
            _maxCnt = maxParticles;
            _size = 0;

            _addKernel = shader.FindKernel("AddParticles");
            _interpKernel = shader.FindKernel("InterpolateVelocity");
            _rk4Kernel = shader.FindKernel("RK4Update");
            _wrapKernel = shader.FindKernel("WrapPositions");
            _clampKernel = shader.FindKernel("ClampPositionsToVolume");

            _constrainLiquidChiKernel = -1;
            if (_chiConstrainShader != null
                && _chiConstrainShader.HasKernel("ConstrainParticlesToLiquidChi"))
                _constrainLiquidChiKernel = _chiConstrainShader.FindKernel("ConstrainParticlesToLiquidChi");

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

        /// <param name="clampVelocityGrid">
        /// Для непериодического домена (ёмкость): индексы сетки без wrap по модулям, иначе у границ
        /// частицы читают скорость с «противоположной» стороны → мигание цвета и почти нулевой дрейф.
        /// </param>
        public void CalculateMovement(CSVelocity vel, bool clampVelocityGrid = false)
        {
            if (_size == 0) return;

            SetTorusUniforms();
            _shader.SetInt("_VelInterpClampGrid", clampVelocityGrid ? 1 : 0);
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

        /// <param name="jitterSigmaInCells">Множитель к min(dx,dy,dz): лёгкий разброс после клампа (трассы без объёма не давят друг друга).</param>
        public void ClampPositionsToVolume(float volSizeX, float volSizeY, float volSizeZ,
            float jitterSigmaInCells = 0f, int jitterSeed = 0)
        {
            if (_size == 0) return;
            float mcell = Mathf.Min(_torDX, Mathf.Min(_torDY, _torDZ));
            float margin = Mathf.Max(1e-5f, 0.35f * mcell);
            float halfMin = 0.5f * Mathf.Min(volSizeX, Mathf.Min(volSizeY, volSizeZ));
            if (margin >= halfMin - 1e-4f)
                margin = Mathf.Max(1e-5f, 0.2f * halfMin);
            float jitter = jitterSigmaInCells > 0f ? jitterSigmaInCells * mcell : 0f;

            _shader.SetFloat("_ClampMargin", margin);
            _shader.SetFloat("_JitterSigma", jitter);
            _shader.SetInt("_JitterSeed", jitterSeed);

            _shader.SetInt("_ParticleCount", _size);
            _shader.SetFloat("_VolSizeX", volSizeX);
            _shader.SetFloat("_VolSizeY", volSizeY);
            _shader.SetFloat("_VolSizeZ", volSizeZ);
            _shader.SetBuffer(_clampKernel, "_PosX", _x);
            _shader.SetBuffer(_clampKernel, "_PosY", _y);
            _shader.SetBuffer(_clampKernel, "_PosZ", _z);
            _shader.Dispatch(_clampKernel, (_size + 255) / 256, 1, 1);
        }

        /// <summary>
        /// Визуальная привязка трассеров к жидкой фазе χ: в газе шаг вдоль ∇χ (сэмпл χ трилинейно) плюс слабое смещение вдоль gravityDir.
        /// </summary>
        public void ConstrainToLiquidChi(ComputeBuffer liquidChi, float threshold,
            float strength = 0f, Vector3 gravityDir = default, float gravityTermPerCell = 0f,
            float chiSoftMargin = 0.1f)
        {
            if (_size == 0 || liquidChi == null || _constrainLiquidChiKernel < 0
                || _chiConstrainShader == null)
                return;

            _chiConstrainShader.SetInt("_TorResX", _torResX);
            _chiConstrainShader.SetInt("_TorResY", _torResY);
            _chiConstrainShader.SetInt("_TorResZ", _torResZ);
            _chiConstrainShader.SetFloat("_TorDX", _torDX);
            _chiConstrainShader.SetFloat("_TorDY", _torDY);
            _chiConstrainShader.SetFloat("_TorDZ", _torDZ);
            _chiConstrainShader.SetInt("_ParticleCount", _size);
            _chiConstrainShader.SetFloat("_ChiThresholdParticles", threshold);
            _chiConstrainShader.SetFloat("_ChiConstrainStrength", Mathf.Clamp01(strength));
            float maxMargin = Mathf.Max(0f, threshold - 0.03f);
            _chiConstrainShader.SetFloat("_ChiConstrainSoftMargin", Mathf.Clamp(chiSoftMargin, 0.05f, maxMargin));
            _chiConstrainShader.SetVector("_ChiGravityDir", gravityDir);
            _chiConstrainShader.SetFloat("_ChiGravityScale", gravityTermPerCell);
            _chiConstrainShader.SetBuffer(_constrainLiquidChiKernel, "_PosX", _x);
            _chiConstrainShader.SetBuffer(_constrainLiquidChiKernel, "_PosY", _y);
            _chiConstrainShader.SetBuffer(_constrainLiquidChiKernel, "_PosZ", _z);
            _chiConstrainShader.SetBuffer(_constrainLiquidChiKernel, "_LiquidChi", liquidChi);
            _chiConstrainShader.Dispatch(_constrainLiquidChiKernel, (_size + 255) / 256, 1, 1);
        }

        public void ReadPositions(float[] outX, float[] outY, float[] outZ)
        {
            if (_size == 0) return;
            _x.GetData(outX, 0, 0, _size);
            _y.GetData(outY, 0, 0, _size);
            _z.GetData(outZ, 0, 0, _size);
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
