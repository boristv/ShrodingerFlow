using System;
using UnityEngine;

namespace ComputeShaderSF
{
    public class CSParticles : IDisposable
    {
        private ComputeShader _shader;
        private ComputeShader _wallShader;
        private ComputeShader _chiConstrainShader;
        private int _addKernel, _interpKernel, _rk4Kernel, _wrapKernel;
        private int _wallInterpKernel, _wallClampKernel;
        private int _wallPushSolidKernel, _wallKillVentKernel;
        private int _constrainLiquidChiKernel;
        private int _buoyantDriftKernel;
        private int _killLowChiKernel;
        private int _disperseKernel;

        private static bool s_warnedWallMissing;

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
        /// <param name="wallParticleShader">Опционально: <c>SFComputeParticlesWall</c> — интерполяция u без wrap по сетке и кламп у границ объёма (ёмкость).</param>
        public void Init(ComputeShader shader, int maxParticles, CSISF isf,
            ComputeShader chiConstrainShader = null,
            ComputeShader wallParticleShader = null)
        {
            _shader = shader;
            _wallShader = wallParticleShader;
            _chiConstrainShader = chiConstrainShader;
            _maxCnt = maxParticles;
            _size = 0;

            _addKernel = shader.FindKernel("AddParticles");
            _interpKernel = shader.FindKernel("InterpolateVelocity");
            _rk4Kernel = shader.FindKernel("RK4Update");
            _wrapKernel = shader.FindKernel("WrapPositions");

            _wallInterpKernel = -1;
            _wallClampKernel = -1;
            _wallPushSolidKernel = -1;
            _wallKillVentKernel = -1;
            if (_wallShader != null)
            {
                _wallInterpKernel = _wallShader.FindKernel("InterpolateVelocity");
                _wallClampKernel = _wallShader.FindKernel("ClampPositionsToVolume");
                if (_wallShader.HasKernel("PushParticlesOutOfSolid"))
                    _wallPushSolidKernel = _wallShader.FindKernel("PushParticlesOutOfSolid");
                if (_wallShader.HasKernel("KillParticlesInVent"))
                    _wallKillVentKernel = _wallShader.FindKernel("KillParticlesInVent");
            }

            _constrainLiquidChiKernel = -1;
            _buoyantDriftKernel = -1;
            if (_chiConstrainShader != null
                && _chiConstrainShader.HasKernel("ConstrainParticlesToLiquidChi"))
                _constrainLiquidChiKernel = _chiConstrainShader.FindKernel("ConstrainParticlesToLiquidChi");
            if (_chiConstrainShader != null
                && _chiConstrainShader.HasKernel("AddBuoyantDriftToParticles"))
                _buoyantDriftKernel = _chiConstrainShader.FindKernel("AddBuoyantDriftToParticles");
            _killLowChiKernel = -1;
            if (_chiConstrainShader != null
                && _chiConstrainShader.HasKernel("KillParticlesLowChi"))
                _killLowChiKernel = _chiConstrainShader.FindKernel("KillParticlesLowChi");
            _disperseKernel = -1;
            if (_chiConstrainShader != null
                && _chiConstrainShader.HasKernel("DisperseParticlesInChi"))
                _disperseKernel = _chiConstrainShader.FindKernel("DisperseParticlesInChi");

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

            bool useWallInterp = clampVelocityGrid && _wallShader != null && _wallInterpKernel >= 0;
            if (clampVelocityGrid && !useWallInterp && !s_warnedWallMissing)
            {
                Debug.LogWarning(
                    "[CSParticles] clampVelocityGrid без SFComputeParticlesWall — остаётся периодическая интерполяция u (артефакты у стёнок ёмкости). Задайте wall-шейдер в SFUnifiedCS / соберите сцену с ассетом.");
                s_warnedWallMissing = true;
            }

            SetTorusUniforms(_shader);
            SetTorusUniforms(_wallShader);
            _shader.SetInt("_ParticleCount", _size);
            if (_wallShader != null)
                _wallShader.SetInt("_ParticleCount", _size);
            int groups = (_size + 255) / 256;

            RunRK4Step(vel, _k1x, _k1y, _k1z, _k1x, _k1y, _k1z, 0f, groups, useWallInterp);
            RunRK4Step(vel, _k1x, _k1y, _k1z, _k2x, _k2y, _k2z, _dt * 0.5f, groups, useWallInterp);
            RunRK4Step(vel, _k2x, _k2y, _k2z, _k3x, _k3y, _k3z, _dt * 0.5f, groups, useWallInterp);
            RunRK4Step(vel, _k3x, _k3y, _k3z, _k4x, _k4y, _k4z, _dt, groups, useWallInterp);

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
            float shiftFactor, int groups, bool useWallInterp)
        {
            ComputeShader interpSh = useWallInterp && _wallShader != null ? _wallShader : _shader;
            int interpK = useWallInterp && _wallShader != null ? _wallInterpKernel : _interpKernel;

            interpSh.SetFloat("_ShiftFactor", shiftFactor);

            interpSh.SetBuffer(interpK, "_PosX", _x);
            interpSh.SetBuffer(interpK, "_PosY", _y);
            interpSh.SetBuffer(interpK, "_PosZ", _z);
            interpSh.SetBuffer(interpK, "_ShiftX", shiftX);
            interpSh.SetBuffer(interpK, "_ShiftY", shiftY);
            interpSh.SetBuffer(interpK, "_ShiftZ", shiftZ);
            interpSh.SetBuffer(interpK, "_VelFieldX", vel.vx);
            interpSh.SetBuffer(interpK, "_VelFieldY", vel.vy);
            interpSh.SetBuffer(interpK, "_VelFieldZ", vel.vz);
            interpSh.SetBuffer(interpK, "_OutX", outX);
            interpSh.SetBuffer(interpK, "_OutY", outY);
            interpSh.SetBuffer(interpK, "_OutZ", outZ);
            interpSh.Dispatch(interpK, groups, 1, 1);
        }

        private void SetTorusUniforms(ComputeShader s)
        {
            if (s == null) return;
            s.SetInt("_TorResX", _torResX);
            s.SetInt("_TorResY", _torResY);
            s.SetInt("_TorResZ", _torResZ);
            s.SetFloat("_TorDX", _torDX);
            s.SetFloat("_TorDY", _torDY);
            s.SetFloat("_TorDZ", _torDZ);
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
            if (_size == 0 || _wallShader == null || _wallClampKernel < 0) return;
            float mcell = Mathf.Min(_torDX, Mathf.Min(_torDY, _torDZ));
            float margin = Mathf.Max(1e-5f, 0.35f * mcell);
            float halfMin = 0.5f * Mathf.Min(volSizeX, Mathf.Min(volSizeY, volSizeZ));
            if (margin >= halfMin - 1e-4f)
                margin = Mathf.Max(1e-5f, 0.2f * halfMin);
            float jitter = jitterSigmaInCells > 0f ? jitterSigmaInCells * mcell : 0f;

            _wallShader.SetFloat("_ClampMargin", margin);
            _wallShader.SetFloat("_JitterSigma", jitter);
            _wallShader.SetInt("_JitterSeed", jitterSeed);

            _wallShader.SetInt("_ParticleCount", _size);
            _wallShader.SetFloat("_VolSizeX", volSizeX);
            _wallShader.SetFloat("_VolSizeY", volSizeY);
            _wallShader.SetFloat("_VolSizeZ", volSizeZ);
            _wallShader.SetBuffer(_wallClampKernel, "_PosX", _x);
            _wallShader.SetBuffer(_wallClampKernel, "_PosY", _y);
            _wallShader.SetBuffer(_wallClampKernel, "_PosZ", _z);
            _wallShader.Dispatch(_wallClampKernel, (_size + 255) / 256, 1, 1);
        }

        /// <summary>Выталкивание трассеров из твёрдых ячеек (стены/перегородка) по маске солида. Только при наличии wall-шейдера.</summary>
        public void PushOutOfSolid(ComputeBuffer solidMask, int searchCells = 4)
        {
            if (_size == 0 || _wallShader == null || _wallPushSolidKernel < 0 || solidMask == null)
                return;
            SetTorusUniforms(_wallShader);
            _wallShader.SetInt("_ParticleCount", _size);
            _wallShader.SetInt("_PushSearchCells", Mathf.Max(1, searchCells));
            _wallShader.SetBuffer(_wallPushSolidKernel, "_SolidMask", solidMask);
            _wallShader.SetBuffer(_wallPushSolidKernel, "_PosX", _x);
            _wallShader.SetBuffer(_wallPushSolidKernel, "_PosY", _y);
            _wallShader.SetBuffer(_wallPushSolidKernel, "_PosZ", _z);
            _wallShader.Dispatch(_wallPushSolidKernel, (_size + 255) / 256, 1, 1);
        }

        /// <summary>
        /// Дрейф всплытия трассеров дыма (drift-flux): частицы внутри дыма (χ ≥ порога) поднимаются вместе с дымом
        /// сквозь спокойный воздух, где поле скорости замаскировано нулём. Требует χ-constrain шейдер.
        /// </summary>
        public void AddBuoyantDrift(ComputeBuffer liquidChi, Vector3 drift, float dt,
            float chiLo, float chiHi)
        {
            if (_size == 0 || liquidChi == null || _buoyantDriftKernel < 0 || _chiConstrainShader == null)
                return;
            _chiConstrainShader.SetInt("_TorResX", _torResX);
            _chiConstrainShader.SetInt("_TorResY", _torResY);
            _chiConstrainShader.SetInt("_TorResZ", _torResZ);
            _chiConstrainShader.SetFloat("_TorDX", _torDX);
            _chiConstrainShader.SetFloat("_TorDY", _torDY);
            _chiConstrainShader.SetFloat("_TorDZ", _torDZ);
            _chiConstrainShader.SetInt("_ParticleCount", _size);
            _chiConstrainShader.SetVector("_Drift", drift);
            _chiConstrainShader.SetFloat("_DriftDT", dt);
            _chiConstrainShader.SetFloat("_DriftChiLo", chiLo);
            _chiConstrainShader.SetFloat("_DriftChiHi", chiHi);
            _chiConstrainShader.SetBuffer(_buoyantDriftKernel, "_PosX", _x);
            _chiConstrainShader.SetBuffer(_buoyantDriftKernel, "_PosY", _y);
            _chiConstrainShader.SetBuffer(_buoyantDriftKernel, "_PosZ", _z);
            _chiConstrainShader.SetBuffer(_buoyantDriftKernel, "_LiquidChi", liquidChi);
            _chiConstrainShader.Dispatch(_buoyantDriftKernel, (_size + 255) / 256, 1, 1);
        }

        /// <summary>Удалить трассеры, ушедшие из дыма (χ ниже порога): убирает статичный «замёрзший» шар в неподвижном воздухе.</summary>
        public void KillLowChi(ComputeBuffer liquidChi, float threshold)
        {
            if (_size == 0 || liquidChi == null || _killLowChiKernel < 0 || _chiConstrainShader == null)
                return;
            _chiConstrainShader.SetInt("_TorResX", _torResX);
            _chiConstrainShader.SetInt("_TorResY", _torResY);
            _chiConstrainShader.SetInt("_TorResZ", _torResZ);
            _chiConstrainShader.SetFloat("_TorDX", _torDX);
            _chiConstrainShader.SetFloat("_TorDY", _torDY);
            _chiConstrainShader.SetFloat("_TorDZ", _torDZ);
            _chiConstrainShader.SetInt("_ParticleCount", _size);
            _chiConstrainShader.SetFloat("_KillChiThreshold", threshold);
            _chiConstrainShader.SetBuffer(_killLowChiKernel, "_PosX", _x);
            _chiConstrainShader.SetBuffer(_killLowChiKernel, "_PosY", _y);
            _chiConstrainShader.SetBuffer(_killLowChiKernel, "_PosZ", _z);
            _chiConstrainShader.SetBuffer(_killLowChiKernel, "_LiquidChi", liquidChi);
            _chiConstrainShader.Dispatch(_killLowChiKernel, (_size + 255) / 256, 1, 1);
        }

        /// <summary>Турбулентная дисперсия трассеров внутри дыма (случайное блуждание, взвешено по χ) — расширяет султан.</summary>
        public void DisperseInChi(ComputeBuffer liquidChi, float sigma, int seed, float threshold)
        {
            if (_size == 0 || liquidChi == null || _disperseKernel < 0 || _chiConstrainShader == null
                || sigma <= 0f)
                return;
            _chiConstrainShader.SetInt("_TorResX", _torResX);
            _chiConstrainShader.SetInt("_TorResY", _torResY);
            _chiConstrainShader.SetInt("_TorResZ", _torResZ);
            _chiConstrainShader.SetFloat("_TorDX", _torDX);
            _chiConstrainShader.SetFloat("_TorDY", _torDY);
            _chiConstrainShader.SetFloat("_TorDZ", _torDZ);
            _chiConstrainShader.SetInt("_ParticleCount", _size);
            _chiConstrainShader.SetFloat("_DisperseSigma", sigma);
            _chiConstrainShader.SetInt("_DisperseSeed", seed);
            _chiConstrainShader.SetFloat("_DisperseChiThreshold", threshold);
            _chiConstrainShader.SetBuffer(_disperseKernel, "_PosX", _x);
            _chiConstrainShader.SetBuffer(_disperseKernel, "_PosY", _y);
            _chiConstrainShader.SetBuffer(_disperseKernel, "_PosZ", _z);
            _chiConstrainShader.SetBuffer(_disperseKernel, "_LiquidChi", liquidChi);
            _chiConstrainShader.Dispatch(_disperseKernel, (_size + 255) / 256, 1, 1);
        }

        /// <summary>Пометить трассеры в зоне вытяжки как «мёртвые» (вынос за границы); удаляются последующим <see cref="CompactParticles"/>.</summary>
        public void KillInVent(Vector3 ventMin, Vector3 ventMax)
        {
            if (_size == 0 || _wallShader == null || _wallKillVentKernel < 0)
                return;
            _wallShader.SetInt("_ParticleCount", _size);
            _wallShader.SetFloat("_VentMinX", ventMin.x);
            _wallShader.SetFloat("_VentMinY", ventMin.y);
            _wallShader.SetFloat("_VentMinZ", ventMin.z);
            _wallShader.SetFloat("_VentMaxX", ventMax.x);
            _wallShader.SetFloat("_VentMaxY", ventMax.y);
            _wallShader.SetFloat("_VentMaxZ", ventMax.z);
            _wallShader.SetBuffer(_wallKillVentKernel, "_PosX", _x);
            _wallShader.SetBuffer(_wallKillVentKernel, "_PosY", _y);
            _wallShader.SetBuffer(_wallKillVentKernel, "_PosZ", _z);
            _wallShader.Dispatch(_wallKillVentKernel, (_size + 255) / 256, 1, 1);
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
