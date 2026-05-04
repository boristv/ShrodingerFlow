using System;
using UnityEngine;

namespace ComputeShaderSF
{
    public class CSISF : IDisposable
    {
        public int resX, resY, resZ, num;
        public float dx, dy, dz;
        public int sizeX, sizeY, sizeZ;
        public float hbar, dt;
        public float Cs = 0.07f;
        public float filterFac = 4.0f;
        /// <summary>Кинематическая вязкость ν: ∂u += ν∇²u·dt на поле скорости после VelocityOneForm; при LES складывается с ν_t.</summary>
        public float kinematicViscosity;
        /// <summary>
        /// Условия на границе сетки без периодического wrap для VelocityOne и Div (домен-ёмкость).
        /// Иначе сосед «+X» у правой стенки = левая стенка → ложная скорость и мигание частиц.
        /// </summary>
        public bool clampGridBorders;

        /// <summary>Маркер жидкости χ∈[0,1] на сетке: адвекция по u и «вакуум» ψ в газе (отдельно от нормировки |ψ| по ячейке).</summary>
        public bool useLiquidChiField;

        /// <summary>Ячейки с χ ниже порога — газ; в них ψ принудительно подменяется на малый вакуум после нормировки и фазовых шагов.</summary>
        public float liquidChiThreshold = 0.5f;

        public float[] pxCPU, pyCPU, pzCPU;
        public ComputeBuffer psi1, psi2;

        private ComputeBuffer _mask, _fac;
        private ComputeBuffer _px, _py, _pz;
        private ComputeBuffer _divResult, _poissonTemp;
        private CSVelocity _velCurrent, _velFiltered, _velTemp;

        private ComputeShader _kernels;
        private CSFFT _fft;
        private CSLES _les;

        private int _normalizeK, _gaugeK, _shiftK, _mulEachK;
        private int _copyR2CK, _fftNormK, _velOneK;
        private int _staggeredK, _divK, _jetK;
        private int _gravK, _heatK, _addGravVelK;
        private int _advectLiquidChiK, _gasVacuumPsiK;
        private int _maskVelChiK;
        private int _sumHorizVelK, _subHorizMeanK;

        private ComputeBuffer _liquidChi;
        private ComputeBuffer _liquidChiTemp;
        private ComputeBuffer _partialVelXZSum;
        private Vector2[] _partialVelXZCpu;
        private int _velHorizMeanGroupCount;

        public void Init(ComputeShader kernels, ComputeShader fftShader,
            ComputeShader lesShader, int[] volSize, int[] volRes, float hbar, float dt)
        {
            sizeX = volSize[0]; sizeY = volSize[1]; sizeZ = volSize[2];
            int rx0 = volRes[0], ry0 = volRes[1], rz0 = volRes[2];
            resX = FloorToPowerOfTwo(rx0);
            resY = FloorToPowerOfTwo(ry0);
            resZ = FloorToPowerOfTwo(rz0);
            if (resX != rx0 || resY != ry0 || resZ != rz0)
            {
                volRes[0] = resX;
                volRes[1] = resY;
                volRes[2] = resZ;
                Debug.LogWarning(
                    $"[CSISF] SFComputeFFT radix-2: vol_res приведено к степеням двойки (было {rx0},{ry0},{rz0} → {resX},{resY},{resZ}). " +
                    "Иначе ψ/FFT ломаются и частицы «пропадают». Сохрани сцену, чтобы зафиксировать значения.");
            }
            num = resX * resY * resZ;
            dx = sizeX / (float)resX;
            dy = sizeY / (float)resY;
            dz = sizeZ / (float)resZ;
            this.hbar = hbar;
            this.dt = dt;

            _kernels = kernels;

            _normalizeK = kernels.FindKernel("Normalize");
            _gaugeK = kernels.FindKernel("Gauge");
            _shiftK = kernels.FindKernel("Shift");
            _mulEachK = kernels.FindKernel("MulEach");
            _copyR2CK = kernels.FindKernel("CopyRealToComplex");
            _fftNormK = kernels.FindKernel("FFTNorm");
            _velOneK = kernels.FindKernel("VelocityOne");
            _staggeredK = kernels.FindKernel("StaggeredSharp");
            _divK = kernels.FindKernel("Div");
            _jetK = kernels.FindKernel("ApplyJetBoundary");
            _gravK = kernels.FindKernel("GravityPsi2");
            _heatK = kernels.FindKernel("HeatSinkPsi1");
            _addGravVelK = kernels.HasKernel("AddGravityToVelocity")
                ? kernels.FindKernel("AddGravityToVelocity")
                : -1;
            _advectLiquidChiK = kernels.FindKernel("AdvectLiquidChi");
            _gasVacuumPsiK = kernels.FindKernel("ApplyGasVacuumPsi");
            _maskVelChiK = kernels.FindKernel("MaskVelocityByLiquidChi");
            _sumHorizVelK = kernels.FindKernel("SumHorizontalVelocityPartial");
            _subHorizMeanK = kernels.FindKernel("SubtractHorizontalVelocityMean");

            psi1 = new ComputeBuffer(num, sizeof(float) * 2);
            psi2 = new ComputeBuffer(num, sizeof(float) * 2);
            _mask = new ComputeBuffer(num, sizeof(float) * 2);
            _fac = new ComputeBuffer(num, sizeof(float) * 2);
            _divResult = new ComputeBuffer(num, sizeof(float));
            _poissonTemp = new ComputeBuffer(num, sizeof(float) * 2);

            _px = new ComputeBuffer(num, sizeof(float));
            _py = new ComputeBuffer(num, sizeof(float));
            _pz = new ComputeBuffer(num, sizeof(float));

            _velCurrent = new CSVelocity(resX, resY, resZ);
            _velFiltered = new CSVelocity(resX, resY, resZ);
            _velTemp = new CSVelocity(resX, resY, resZ);

            _liquidChi = new ComputeBuffer(num, sizeof(float));
            _liquidChiTemp = new ComputeBuffer(num, sizeof(float));

            _velHorizMeanGroupCount = (num + 255) / 256;
            _partialVelXZSum = new ComputeBuffer(_velHorizMeanGroupCount, sizeof(float) * 2);
            _partialVelXZCpu = new Vector2[_velHorizMeanGroupCount];
            _fft = new CSFFT();
            _fft.Init(fftShader, resX, resY, resZ);

            _les = new CSLES();
            _les.Init(lesShader, resX, resY, resZ);

            BuildPositionGrids();
            BuildMask();
            BuildFac();
        }

        private static bool IsPowerOfTwo(int n) => n > 0 && (n & (n - 1)) == 0;

        /// <summary>Наибольшая степень двойки ≤ n (для radix-2 FFT; например 192 → 128).</summary>
        private static int FloorToPowerOfTwo(int n)
        {
            if (n <= 0) return 1;
            if (IsPowerOfTwo(n)) return n;
            int p = 1;
            while (p * 2 <= n) p *= 2;
            return p;
        }

        private void BuildPositionGrids()
        {
            pxCPU = new float[num];
            pyCPU = new float[num];
            pzCPU = new float[num];
            for (int i = 0; i < resX; i++)
                for (int j = 0; j < resY; j++)
                    for (int k = 0; k < resZ; k++)
                    {
                        int idx = i * resY * resZ + j * resZ + k;
                        pxCPU[idx] = i * dx;
                        pyCPU[idx] = j * dy;
                        pzCPU[idx] = k * dz;
                    }
            _px.SetData(pxCPU);
            _py.SetData(pyCPU);
            _pz.SetData(pzCPU);
        }

        private void BuildMask()
        {
            float fac = -4f * Mathf.PI * Mathf.PI * hbar;
            var data = new Vector2[num];
            for (int i = 0; i < resX; i++)
                for (int j = 0; j < resY; j++)
                    for (int k = 0; k < resZ; k++)
                    {
                        float kx = (i - resX / 2f) / sizeX;
                        float ky = (j - resY / 2f) / sizeY;
                        float kz = (k - resZ / 2f) / sizeZ;
                        float lambda = fac * (kx * kx + ky * ky + kz * kz);
                        float phase = lambda * dt / 2f;
                        int idx = i * resY * resZ + j * resZ + k;
                        data[idx] = new Vector2(Mathf.Cos(phase), Mathf.Sin(phase));
                    }
            _mask.SetData(data);
        }

        private void BuildFac()
        {
            var data = new Vector2[num];
            for (int i = 0; i < resX; i++)
                for (int j = 0; j < resY; j++)
                    for (int k = 0; k < resZ; k++)
                    {
                        float sx = Mathf.Sin(Mathf.PI * i / resX) / dx;
                        float sy = Mathf.Sin(Mathf.PI * j / resY) / dy;
                        float sz = Mathf.Sin(Mathf.PI * k / resZ) / dz;
                        float denom = sx * sx + sy * sy + sz * sz;
                        int idx = i * resY * resZ + j * resZ + k;
                        data[idx] = (i == 0 && j == 0 && k == 0)
                            ? Vector2.zero
                            : new Vector2(-0.25f / denom, 0);
                    }
            _fac.SetData(data);
        }

        private int Groups1D => (num + 255) / 256;

        private void SetCommonUniforms()
        {
            _kernels.SetInt("_ResX", resX);
            _kernels.SetInt("_ResY", resY);
            _kernels.SetInt("_ResZ", resZ);
            _kernels.SetInt("_Num", num);
            _kernels.SetFloat("_DX", dx);
            _kernels.SetFloat("_DY", dy);
            _kernels.SetFloat("_DZ", dz);
            _kernels.SetFloat("_Hbar", hbar);
            _kernels.SetInt("_ClampGridBorders", clampGridBorders ? 1 : 0);
        }

        public void Normalize()
        {
            SetCommonUniforms();
            _kernels.SetBuffer(_normalizeK, "_Psi1", psi1);
            _kernels.SetBuffer(_normalizeK, "_Psi2", psi2);
            _kernels.Dispatch(_normalizeK, Groups1D, 1, 1);
        }

        /// <summary>Текущее поле χ (после последней адвекции); null, если <see cref="useLiquidChiField"/> выключен — для визуализации.</summary>
        public ComputeBuffer LiquidChiBuffer => useLiquidChiField ? _liquidChi : null;

        /// <summary>Загрузить χ на GPU; копия и во временный буфер для первого шага адвекции.</summary>
        public void UploadLiquidChi(float[] chi)
        {
            if (chi == null || chi.Length != num)
                throw new ArgumentException($"[CSISF] UploadLiquidChi: ожидалось {num} значений.", nameof(chi));
            _liquidChi.SetData(chi);
            _liquidChiTemp.SetData(chi);
        }

        /// <summary>
        /// Адвекция χ полулагранжевски по полю u до маски «газа» (в ёмкости — после <see cref="ApplyVelocityGravity"/>,
        /// если она используется), иначе интерфейс χ не движется; затем вызывайте <see cref="ApplyLiquidChiVelocityMask"/>.
        /// CFL ограничен в шейдере.
        /// </summary>
        public void AdvectLiquidChi(CSVelocity vel)
        {
            if (!useLiquidChiField)
                return;
            SetCommonUniforms();
            _kernels.SetFloat("_DT", dt);
            _kernels.SetFloat("_BoundX", sizeX);
            _kernels.SetFloat("_BoundY", sizeY);
            _kernels.SetFloat("_BoundZ", sizeZ);
            _kernels.SetBuffer(_advectLiquidChiK, "_ChiIn", _liquidChi);
            _kernels.SetBuffer(_advectLiquidChiK, "_ChiOut", _liquidChiTemp);
            _kernels.SetBuffer(_advectLiquidChiK, "_VX", vel.vx);
            _kernels.SetBuffer(_advectLiquidChiK, "_VY", vel.vy);
            _kernels.SetBuffer(_advectLiquidChiK, "_VZ", vel.vz);
            _kernels.SetBuffer(_advectLiquidChiK, "_PX", _px);
            _kernels.SetBuffer(_advectLiquidChiK, "_PY", _py);
            _kernels.SetBuffer(_advectLiquidChiK, "_PZ", _pz);
            _kernels.Dispatch(_advectLiquidChiK, Groups1D, 1, 1);
            var swap = _liquidChi;
            _liquidChi = _liquidChiTemp;
            _liquidChiTemp = swap;
        }

        private void ApplyGasVacuumFromChi()
        {
            if (!useLiquidChiField)
                return;
            SetCommonUniforms();
            _kernels.SetFloat("_ChiLiqThreshold", liquidChiThreshold);
            const float e1 = 1e-6f, e2 = 1e-7f;
            _kernels.SetFloat("_Vac1R", e1);
            _kernels.SetFloat("_Vac1I", 0f);
            _kernels.SetFloat("_Vac2R", e2);
            _kernels.SetFloat("_Vac2I", 0f);
            _kernels.SetBuffer(_gasVacuumPsiK, "_ChiField", _liquidChi);
            _kernels.SetBuffer(_gasVacuumPsiK, "_Psi1", psi1);
            _kernels.SetBuffer(_gasVacuumPsiK, "_Psi2", psi2);
            _kernels.Dispatch(_gasVacuumPsiK, Groups1D, 1, 1);
        }

        /// <summary>
        /// Обнуляет компоненты u, если соответствующая разность ψ тянется через «газ» (низкий χ).
        /// Иначе |ψ_жидк|≫|ψ_газ| после нормировки + вакуума даёт взрывные фазовые градиенты и унос χ/частиц.
        /// </summary>
        private void MaskVelocityByLiquidChi(CSVelocity v)
        {
            if (!useLiquidChiField)
                return;
            SetCommonUniforms();
            _kernels.SetFloat("_ChiLiqThreshold", liquidChiThreshold);
            _kernels.SetBuffer(_maskVelChiK, "_ChiField", _liquidChi);
            _kernels.SetBuffer(_maskVelChiK, "_VX", v.vx);
            _kernels.SetBuffer(_maskVelChiK, "_VY", v.vy);
            _kernels.SetBuffer(_maskVelChiK, "_VZ", v.vz);
            _kernels.Dispatch(_maskVelChiK, Groups1D, 1, 1);
        }

        /// <summary>
        /// При непериодических границах VelocityOne даёт несимметричный шаблон (ix=0 vs ix=ResX-1) → ненулевое среднее Vx/Vz и дрейф «вбок».
        /// Снимаем постоянную составляющую по X и Z (по Y не трогаем — гравитация). Только при <see cref="clampGridBorders"/>.
        /// </summary>
        private void RemoveMeanHorizontalVelocity(CSVelocity v)
        {
            if (!clampGridBorders)
                return;
            SetCommonUniforms();
            _kernels.SetBuffer(_sumHorizVelK, "_VX", v.vx);
            _kernels.SetBuffer(_sumHorizVelK, "_VZ", v.vz);
            _kernels.SetBuffer(_sumHorizVelK, "_PartialVelXZSum", _partialVelXZSum);
            _kernels.Dispatch(_sumHorizVelK, _velHorizMeanGroupCount, 1, 1);

            _partialVelXZSum.GetData(_partialVelXZCpu);
            double sx = 0.0, sz = 0.0;
            for (int i = 0; i < _velHorizMeanGroupCount; i++)
            {
                sx += _partialVelXZCpu[i].x;
                sz += _partialVelXZCpu[i].y;
            }
            float mx = (float)(sx / num);
            float mz = (float)(sz / num);
            _kernels.SetFloat("_MeanSubX", mx);
            _kernels.SetFloat("_MeanSubZ", mz);
            _kernels.SetBuffer(_subHorizMeanK, "_VX", v.vx);
            _kernels.SetBuffer(_subHorizMeanK, "_VZ", v.vz);
            _kernels.Dispatch(_subHorizMeanK, Groups1D, 1, 1);
        }

        private void SchoedingerFlow()
        {
            SetCommonUniforms();

            _fft.FFT3D(psi1, false);
            _fft.FFT3D(psi2, false);

            ShiftBuffer(psi1);
            ShiftBuffer(psi2);

            MulEachBuffers(psi1, _mask);
            MulEachBuffers(psi2, _mask);

            ShiftBuffer(psi1);
            ShiftBuffer(psi2);

            _fft.FFT3D(psi1, true);
            _fft.FFT3D(psi2, true);

            FFTNormBuffer(psi1);
            FFTNormBuffer(psi2);
        }

        private void ShiftBuffer(ComputeBuffer buf)
        {
            _kernels.SetBuffer(_shiftK, "_BufComplex", buf);
            _kernels.Dispatch(_shiftK, Groups1D, 1, 1);
        }

        private void MulEachBuffers(ComputeBuffer a, ComputeBuffer b)
        {
            _kernels.SetBuffer(_mulEachK, "_BufComplex", a);
            _kernels.SetBuffer(_mulEachK, "_BufComplex2", b);
            _kernels.Dispatch(_mulEachK, Groups1D, 1, 1);
        }

        private void FFTNormBuffer(ComputeBuffer buf)
        {
            _kernels.SetBuffer(_fftNormK, "_BufComplex", buf);
            _kernels.Dispatch(_fftNormK, Groups1D, 1, 1);
        }

        /// <summary>
        /// Скорость из градиента фазы Madelunga: множитель должен совпадать с ℏ во всех местах (LES, PP, частицы).
        /// Раньше здесь было 1.0 — проекция давления видела поле в ~1/ℏ раз «не то», чем то, что извлекается для частиц → поломанный gauge и дрейф (в т.ч. вбок).
        /// </summary>
        private void VelocityOneForm(CSVelocity v)
        {
            VelocityOneForm(v, hbar);
        }

        private void VelocityOneForm(CSVelocity v, float h)
        {
            SetCommonUniforms();
            _kernels.SetFloat("_Hbar", h);
            _kernels.SetBuffer(_velOneK, "_Psi1", psi1);
            _kernels.SetBuffer(_velOneK, "_Psi2", psi2);
            _kernels.SetBuffer(_velOneK, "_VX", v.vx);
            _kernels.SetBuffer(_velOneK, "_VY", v.vy);
            _kernels.SetBuffer(_velOneK, "_VZ", v.vz);
            _kernels.Dispatch(_velOneK, Groups1D, 1, 1);
        }

        private void StaggeredSharp(CSVelocity vel)
        {
            SetCommonUniforms();

            _kernels.SetBuffer(_staggeredK, "_BufFloat", vel.vx);
            _kernels.SetFloat("_Param", dx);
            _kernels.Dispatch(_staggeredK, Groups1D, 1, 1);

            _kernels.SetBuffer(_staggeredK, "_BufFloat", vel.vy);
            _kernels.SetFloat("_Param", dy);
            _kernels.Dispatch(_staggeredK, Groups1D, 1, 1);

            _kernels.SetBuffer(_staggeredK, "_BufFloat", vel.vz);
            _kernels.SetFloat("_Param", dz);
            _kernels.Dispatch(_staggeredK, Groups1D, 1, 1);
        }

        private void Div(CSVelocity v, ComputeBuffer result)
        {
            SetCommonUniforms();
            _kernels.SetBuffer(_divK, "_VX", v.vx);
            _kernels.SetBuffer(_divK, "_VY", v.vy);
            _kernels.SetBuffer(_divK, "_VZ", v.vz);
            _kernels.SetBuffer(_divK, "_BufFloat", result);
            _kernels.Dispatch(_divK, Groups1D, 1, 1);
        }

        private void PoissonSolve(ComputeBuffer f, ComputeBuffer result)
        {
            SetCommonUniforms();

            _kernels.SetBuffer(_copyR2CK, "_BufComplex", result);
            _kernels.SetBuffer(_copyR2CK, "_BufFloat", f);
            _kernels.Dispatch(_copyR2CK, Groups1D, 1, 1);

            _fft.FFT3D(result, false);
            MulEachBuffers(result, _fac);
            _fft.FFT3D(result, true);
        }

        public void PressureProject(CSVelocity v)
        {
            Div(v, _divResult);
            PoissonSolve(_divResult, _poissonTemp);

            SetCommonUniforms();
            _kernels.SetBuffer(_gaugeK, "_Psi1", psi1);
            _kernels.SetBuffer(_gaugeK, "_Psi2", psi2);
            _kernels.SetBuffer(_gaugeK, "_BufComplex", _poissonTemp);
            _kernels.Dispatch(_gaugeK, Groups1D, 1, 1);
        }

        public void PressureProject()
        {
            VelocityOneForm(_velTemp);
            MaskVelocityByLiquidChi(_velTemp);
            RemoveMeanHorizontalVelocity(_velTemp);
            PressureProject(_velTemp);
        }

        private void LESStep()
        {
            VelocityOneForm(_velCurrent);
            MaskVelocityByLiquidChi(_velCurrent);
            _les.Filter(_velCurrent, _velFiltered, dx, dy, dz);
            _les.ComputeNuT(_velFiltered, dx, dy, dz, Cs, filterFac);
            _les.ApplyViscosity(_velCurrent, dx, dy, dz, dt, kinematicViscosity, 1f);
            MaskVelocityByLiquidChi(_velCurrent);
            RemoveMeanHorizontalVelocity(_velCurrent);
        }

        private void LaminarViscosityStep()
        {
            if (kinematicViscosity <= 1e-20f)
                return;
            VelocityOneForm(_velCurrent);
            MaskVelocityByLiquidChi(_velCurrent);
            _les.ApplyViscosity(_velCurrent, dx, dy, dz, dt, kinematicViscosity, 0f);
            MaskVelocityByLiquidChi(_velCurrent);
            RemoveMeanHorizontalVelocity(_velCurrent);
        }

        private bool UsesVelocityFieldViscosity =>
            kinematicViscosity > 1e-20f;

        /// <param name="gravityRotatePsi1Too">
        /// Если true — та же фаза на ψ₁ и ψ₂ (нужно для «тяжёлой» жидкости: иначе |ψ₁|≫|ψ₂| и VelocityOne почти не видит градиент от гравитации).
        /// Cigarette / hip: оставить false (только ψ₂).
        /// </param>
        public void UpdateSpace(bool useLES = true, Vector3? gravityPsi2 = null,
            bool gravityRotatePsi1Too = false)
        {
            SetCommonUniforms();
            SchoedingerFlow();
            if (useLES)
                LESStep();
            else
                LaminarViscosityStep();
            Normalize();
            ApplyGasVacuumFromChi();
            if (gravityPsi2.HasValue && gravityPsi2.Value.sqrMagnitude > 1e-20f)
                ApplyGravityPsi2(gravityPsi2.Value, gravityRotatePsi1Too);
            ApplyGasVacuumFromChi();
            bool ppFromVel = useLES || UsesVelocityFieldViscosity;
            if (ppFromVel)
                PressureProject(_velCurrent);
            else
                PressureProject();
            ApplyGasVacuumFromChi();
        }

        /// <summary>
        /// Цепочка как в example_cigarette.hip (DOP): после Шрёдингера — нормировка, фаза g·P на ψ₂,
        /// давление, обнуление ψ₁ в маске «сигареты», снова нормировка и второй PP.
        /// </summary>
        public void UpdateCigaretteSpace(bool useLES, Vector3 gravityG, ComputeBuffer isJetMask)
        {
            SetCommonUniforms();
            SchoedingerFlow();
            if (useLES)
                LESStep();
            else
                LaminarViscosityStep();
            Normalize();
            ApplyGasVacuumFromChi();
            ApplyGravityPsi2(gravityG);
            ApplyGasVacuumFromChi();
            bool ppFromVel = useLES || UsesVelocityFieldViscosity;
            if (ppFromVel)
                PressureProject(_velCurrent);
            else
                PressureProject();
            ApplyGasVacuumFromChi();
            ApplyHeatSinkPsi1(isJetMask);
            Normalize();
            ApplyGasVacuumFromChi();
            if (ppFromVel)
                PressureProject(_velCurrent);
            else
                PressureProject();
            ApplyGasVacuumFromChi();
        }

        private void ApplyGravityPsi2(Vector3 g, bool rotatePsi1Too = false)
        {
            SetCommonUniforms();
            _kernels.SetInt("_GravityBothPsi", rotatePsi1Too ? 1 : 0);
            _kernels.SetFloat("_GX", g.x);
            _kernels.SetFloat("_GY", g.y);
            _kernels.SetFloat("_GZ", g.z);
            _kernels.SetFloat("_DT", dt);
            _kernels.SetBuffer(_gravK, "_Psi1", psi1);
            _kernels.SetBuffer(_gravK, "_Psi2", psi2);
            _kernels.SetBuffer(_gravK, "_PX", _px);
            _kernels.SetBuffer(_gravK, "_PY", _py);
            _kernels.SetBuffer(_gravK, "_PZ", _pz);
            _kernels.Dispatch(_gravK, Groups1D, 1, 1);
        }

        private void ApplyHeatSinkPsi1(ComputeBuffer isJet)
        {
            SetCommonUniforms();
            _kernels.SetBuffer(_heatK, "_Psi1", psi1);
            _kernels.SetBuffer(_heatK, "_IsJet", isJet);
            _kernels.Dispatch(_heatK, Groups1D, 1, 1);
        }

        /// <summary>Пересчёт u из ψ (без маски по χ): для адвекции χ нужна полная скорость на границе жидкость/газ.</summary>
        public void UpdateVelocities(CSVelocity vel)
        {
            VelocityOneForm(vel);
            StaggeredSharp(vel);
            RemoveMeanHorizontalVelocity(vel);
        }

        /// <summary>
        /// Гравитация на уровне поля скорости: добавляет g·dt к каждой компоненте, затем <see cref="PressureProject(CSVelocity)"/>.
        /// В отличие от GravityPsi2, не накапливает фазу в ψ. После Gauge снова вызывается <see cref="ApplyGasVacuumFromChi"/>,
        /// чтобы в «газе» χ не остались искажённые ψ после фазового сдвига.
        /// Маску по χ НЕ применяет — её вызывает владелец после <see cref="AdvectLiquidChi"/> (чтобы χ двигался на интерфейсе).
        /// </summary>
        public void ApplyVelocityGravity(Vector3 g, CSVelocity vel)
        {
            if (_addGravVelK < 0)
            {
                Debug.LogError("[CSISF] Шейдер ISF без kernel AddGravityToVelocity — обновите SFComputeKernels.compute / asset.");
                return;
            }
            SetCommonUniforms();
            _kernels.SetFloat("_GX", g.x * dt);
            _kernels.SetFloat("_GY", g.y * dt);
            _kernels.SetFloat("_GZ", g.z * dt);
            _kernels.SetBuffer(_addGravVelK, "_VX", vel.vx);
            _kernels.SetBuffer(_addGravVelK, "_VY", vel.vy);
            _kernels.SetBuffer(_addGravVelK, "_VZ", vel.vz);
            _kernels.Dispatch(_addGravVelK, Groups1D, 1, 1);
            RemoveMeanHorizontalVelocity(vel);
            PressureProject(vel);
            ApplyGasVacuumFromChi();
        }

        /// <summary>Обнулить скорость в «газе» по текущему χ — только при <see cref="useLiquidChiField"/>; после <see cref="AdvectLiquidChi"/> и перед трассерами.</summary>
        public void ApplyLiquidChiVelocityMask(CSVelocity vel)
        {
            if (!useLiquidChiField)
                return;
            MaskVelocityByLiquidChi(vel);
        }

        public void ApplyJetBoundary(ComputeBuffer isJet,
            float kvecX, float kvecY, float kvecZ, float phaseOffset)
        {
            SetCommonUniforms();
            _kernels.SetFloat("_KVecX", kvecX);
            _kernels.SetFloat("_KVecY", kvecY);
            _kernels.SetFloat("_KVecZ", kvecZ);
            _kernels.SetFloat("_PhaseOffset", phaseOffset);

            _kernels.SetBuffer(_jetK, "_Psi1", psi1);
            _kernels.SetBuffer(_jetK, "_Psi2", psi2);
            _kernels.SetBuffer(_jetK, "_IsJet", isJet);
            _kernels.SetBuffer(_jetK, "_PX", _px);
            _kernels.SetBuffer(_jetK, "_PY", _py);
            _kernels.SetBuffer(_jetK, "_PZ", _pz);
            _kernels.Dispatch(_jetK, Groups1D, 1, 1);
        }

        public void Dispose()
        {
            psi1?.Release();
            psi2?.Release();
            _mask?.Release();
            _fac?.Release();
            _px?.Release();
            _py?.Release();
            _pz?.Release();
            _divResult?.Release();
            _poissonTemp?.Release();
            _velCurrent?.Dispose();
            _velFiltered?.Dispose();
            _velTemp?.Dispose();
            _liquidChi?.Release();
            _liquidChiTemp?.Release();
            _partialVelXZSum?.Release();
            _les?.Dispose();
        }
    }
}
