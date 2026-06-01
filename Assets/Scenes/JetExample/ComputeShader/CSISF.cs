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

        /// <summary>
        /// χ влияет на ψ (вакуум в «газе») — нужно для жидкость/газ (ёмкость). Для пассивного скаляра дыма выключить:
        /// тогда χ адвектируется как маркер концентрации, но ψ и поле скорости не «обнуляются» вне дыма.
        /// </summary>
        public bool chiAffectsPsiVacuum = true;

        /// <summary>Обнулять u в ячейках с низким χ. Для жидкость/газ — true; для дыма (среда заполняет весь объём) — false.</summary>
        public bool maskVelocityWithChi = true;

        public float[] pxCPU, pyCPU, pzCPU;
        public ComputeBuffer psi1, psi2;

        private ComputeBuffer _mask, _fac;
        private ComputeBuffer _px, _py, _pz;
        private ComputeBuffer _divResult, _poissonTemp;
        private CSVelocity _velCurrent, _velFiltered, _velTemp;

        /// <summary>Базовый ISF (Jet, кольца, …): без χ/стенок, тот же bytecode что и до ветки ёмкости.</summary>
        private ComputeShader _kCore;
        /// <summary>Расширенный ISF (стенки сетки, χ, среднее u): только при <see cref="clampGridBorders"/>.</summary>
        private ComputeShader _kExt;
        private CSFFT _fft;
        private CSLES _les;

        private struct IsfKernelTable
        {
            public int Normalize, Gauge, Shift, MulEach, CopyR2C, FFTNorm, VelOne, Staggered, Div, Jet, Grav, Heat;
            public int AddGravVel, AdvectChi, GasVacuum, MaskVelChi, SumHoriz, SubHorizMean;
            public int Buoyancy, InjectChi, ChiSink, DiffuseChi, CurlNoise;
        }

        private IsfKernelTable _coreIds, _extIds;
        private ComputeShader _activeSh;
        private IsfKernelTable _activeIds;
        private bool _warnedMissingContainerExt;

        private ComputeBuffer _liquidChi;
        private ComputeBuffer _liquidChiTemp;
        private ComputeBuffer _partialVelXZSum;
        private Vector2[] _partialVelXZCpu;
        private int _velHorizMeanGroupCount;

        public void Init(ComputeShader coreKernels, ComputeShader fftShader,
            ComputeShader lesShader, int[] volSize, int[] volRes, float hbar, float dt,
            ComputeShader containerKernels = null)
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

            _kCore = coreKernels;
            _kExt = containerKernels;
            BindKernelTable(_kCore, ref _coreIds, extended: false);
            if (_kExt != null)
                BindKernelTable(_kExt, ref _extIds, extended: true);
            else
                ClearExtKernelSlots(ref _extIds);

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

        private static void BindKernelTable(ComputeShader s, ref IsfKernelTable t, bool extended)
        {
            t.Normalize = s.FindKernel("Normalize");
            t.Gauge = s.FindKernel("Gauge");
            t.Shift = s.FindKernel("Shift");
            t.MulEach = s.FindKernel("MulEach");
            t.CopyR2C = s.FindKernel("CopyRealToComplex");
            t.FFTNorm = s.FindKernel("FFTNorm");
            t.VelOne = s.FindKernel("VelocityOne");
            t.Staggered = s.FindKernel("StaggeredSharp");
            t.Div = s.FindKernel("Div");
            t.Jet = s.FindKernel("ApplyJetBoundary");
            t.Grav = s.FindKernel("GravityPsi2");
            t.Heat = s.FindKernel("HeatSinkPsi1");
            if (extended)
            {
                t.AddGravVel = s.FindKernel("AddGravityToVelocity");
                t.AdvectChi = s.FindKernel("AdvectLiquidChi");
                t.GasVacuum = s.FindKernel("ApplyGasVacuumPsi");
                t.MaskVelChi = s.FindKernel("MaskVelocityByLiquidChi");
                t.SumHoriz = s.FindKernel("SumHorizontalVelocityPartial");
                t.SubHorizMean = s.FindKernel("SubtractHorizontalVelocityMean");
                t.Buoyancy = s.FindKernel("AddBuoyancyToVelocity");
                t.InjectChi = s.FindKernel("InjectChiSource");
                t.ChiSink = s.FindKernel("ChiSinkVent");
                t.DiffuseChi = s.FindKernel("DiffuseChi");
                t.CurlNoise = s.FindKernel("AddCurlNoiseTurbulence");
            }
            else
                ClearExtKernelSlots(ref t);
        }

        private static void ClearExtKernelSlots(ref IsfKernelTable t)
        {
            t.AddGravVel = -1;
            t.AdvectChi = -1;
            t.GasVacuum = -1;
            t.MaskVelChi = -1;
            t.SumHoriz = -1;
            t.SubHorizMean = -1;
            t.Buoyancy = -1;
            t.InjectChi = -1;
            t.ChiSink = -1;
            t.DiffuseChi = -1;
            t.CurlNoise = -1;
        }

        private void PickActive()
        {
            if (clampGridBorders && _kExt == null && !_warnedMissingContainerExt)
            {
                Debug.LogWarning(
                    "[CSISF] clampGridBorders включён, но не задан SFComputeKernelsContainer — используется базовый ISF с периодическими границами (ψ/u и трассеры будут вести себя неверно в ёмкости).");
                _warnedMissingContainerExt = true;
            }
            bool useExt = clampGridBorders && _kExt != null;
            _activeSh = useExt ? _kExt : _kCore;
            _activeIds = useExt ? _extIds : _coreIds;
        }

        private void SetUniformsBase(ComputeShader s)
        {
            s.SetInt("_ResX", resX);
            s.SetInt("_ResY", resY);
            s.SetInt("_ResZ", resZ);
            s.SetInt("_Num", num);
            s.SetFloat("_DX", dx);
            s.SetFloat("_DY", dy);
            s.SetFloat("_DZ", dz);
            s.SetFloat("_Hbar", hbar);
        }

        private void SetExtLayoutUniforms(ComputeShader s)
        {
            s.SetInt("_ClampGridBorders", clampGridBorders ? 1 : 0);
            s.SetInt("_LiquidChiNormalizeFallback", useLiquidChiField ? 1 : 0);
        }

        private void BindExtKernelUniforms()
        {
            if (_kExt == null)
                return;
            SetUniformsBase(_kExt);
            SetExtLayoutUniforms(_kExt);
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
            PickActive();
            SetUniformsBase(_activeSh);
            if (ReferenceEquals(_activeSh, _kExt) && _kExt != null)
                SetExtLayoutUniforms(_activeSh);
        }

        public void Normalize()
        {
            SetCommonUniforms();
            _activeSh.SetBuffer(_activeIds.Normalize, "_Psi1", psi1);
            _activeSh.SetBuffer(_activeIds.Normalize, "_Psi2", psi2);
            _activeSh.Dispatch(_activeIds.Normalize, Groups1D, 1, 1);
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
            AdvectLiquidChi(vel, Vector3.zero);
        }

        /// <param name="drift">Скорость всплытия дыма относительно воздуха (drift-flux): добавляется к u при переносе χ.</param>
        public void AdvectLiquidChi(CSVelocity vel, Vector3 drift)
        {
            if (!useLiquidChiField || _kExt == null || _extIds.AdvectChi < 0)
                return;
            BindExtKernelUniforms();
            _kExt.SetFloat("_ChiDriftX", drift.x);
            _kExt.SetFloat("_ChiDriftY", drift.y);
            _kExt.SetFloat("_ChiDriftZ", drift.z);
            _kExt.SetFloat("_DT", dt);
            _kExt.SetFloat("_BoundX", sizeX);
            _kExt.SetFloat("_BoundY", sizeY);
            _kExt.SetFloat("_BoundZ", sizeZ);
            _kExt.SetBuffer(_extIds.AdvectChi, "_ChiIn", _liquidChi);
            _kExt.SetBuffer(_extIds.AdvectChi, "_ChiOut", _liquidChiTemp);
            _kExt.SetBuffer(_extIds.AdvectChi, "_VX", vel.vx);
            _kExt.SetBuffer(_extIds.AdvectChi, "_VY", vel.vy);
            _kExt.SetBuffer(_extIds.AdvectChi, "_VZ", vel.vz);
            _kExt.SetBuffer(_extIds.AdvectChi, "_PX", _px);
            _kExt.SetBuffer(_extIds.AdvectChi, "_PY", _py);
            _kExt.SetBuffer(_extIds.AdvectChi, "_PZ", _pz);
            _kExt.Dispatch(_extIds.AdvectChi, Groups1D, 1, 1);
            var swap = _liquidChi;
            _liquidChi = _liquidChiTemp;
            _liquidChiTemp = swap;
        }

        private void ApplyGasVacuumFromChi()
        {
            if (!useLiquidChiField || !chiAffectsPsiVacuum || _kExt == null || _extIds.GasVacuum < 0)
                return;
            BindExtKernelUniforms();
            _kExt.SetFloat("_ChiLiqThreshold", liquidChiThreshold);
            const float e1 = 1e-6f, e2 = 1e-7f;
            _kExt.SetFloat("_Vac1R", e1);
            _kExt.SetFloat("_Vac1I", 0f);
            _kExt.SetFloat("_Vac2R", e2);
            _kExt.SetFloat("_Vac2I", 0f);
            _kExt.SetBuffer(_extIds.GasVacuum, "_ChiField", _liquidChi);
            _kExt.SetBuffer(_extIds.GasVacuum, "_Psi1", psi1);
            _kExt.SetBuffer(_extIds.GasVacuum, "_Psi2", psi2);
            _kExt.Dispatch(_extIds.GasVacuum, Groups1D, 1, 1);
        }

        /// <summary>
        /// Обнуляет компоненты u, если соответствующая разность ψ тянется через «газ» (низкий χ).
        /// Иначе |ψ_жидк|≫|ψ_газ| после нормировки + вакуума даёт взрывные фазовые градиенты и унос χ/частиц.
        /// </summary>
        private void MaskVelocityByLiquidChi(CSVelocity v)
        {
            if (!useLiquidChiField || !maskVelocityWithChi || _kExt == null || _extIds.MaskVelChi < 0)
                return;
            BindExtKernelUniforms();
            _kExt.SetFloat("_ChiLiqThreshold", liquidChiThreshold);
            _kExt.SetBuffer(_extIds.MaskVelChi, "_ChiField", _liquidChi);
            _kExt.SetBuffer(_extIds.MaskVelChi, "_VX", v.vx);
            _kExt.SetBuffer(_extIds.MaskVelChi, "_VY", v.vy);
            _kExt.SetBuffer(_extIds.MaskVelChi, "_VZ", v.vz);
            _kExt.Dispatch(_extIds.MaskVelChi, Groups1D, 1, 1);
        }

        /// <summary>
        /// При непериодических границах VelocityOne даёт несимметричный шаблон (ix=0 vs ix=ResX-1) → ненулевое среднее Vx/Vz и дрейф «вбок».
        /// Снимаем постоянную составляющую по X и Z (по Y не трогаем — гравитация). Только при <see cref="clampGridBorders"/>.
        /// </summary>
        private void RemoveMeanHorizontalVelocity(CSVelocity v)
        {
            if (!clampGridBorders || _kExt == null || _extIds.SumHoriz < 0)
                return;
            BindExtKernelUniforms();
            _kExt.SetBuffer(_extIds.SumHoriz, "_VX", v.vx);
            _kExt.SetBuffer(_extIds.SumHoriz, "_VZ", v.vz);
            _kExt.SetBuffer(_extIds.SumHoriz, "_PartialVelXZSum", _partialVelXZSum);
            _kExt.Dispatch(_extIds.SumHoriz, _velHorizMeanGroupCount, 1, 1);

            _partialVelXZSum.GetData(_partialVelXZCpu);
            double sx = 0.0, sz = 0.0;
            for (int i = 0; i < _velHorizMeanGroupCount; i++)
            {
                sx += _partialVelXZCpu[i].x;
                sz += _partialVelXZCpu[i].y;
            }
            float mx = (float)(sx / num);
            float mz = (float)(sz / num);
            _kExt.SetFloat("_MeanSubX", mx);
            _kExt.SetFloat("_MeanSubZ", mz);
            _kExt.SetBuffer(_extIds.SubHorizMean, "_VX", v.vx);
            _kExt.SetBuffer(_extIds.SubHorizMean, "_VZ", v.vz);
            _kExt.Dispatch(_extIds.SubHorizMean, Groups1D, 1, 1);
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
            _activeSh.SetBuffer(_activeIds.Shift, "_BufComplex", buf);
            _activeSh.Dispatch(_activeIds.Shift, Groups1D, 1, 1);
        }

        private void MulEachBuffers(ComputeBuffer a, ComputeBuffer b)
        {
            _activeSh.SetBuffer(_activeIds.MulEach, "_BufComplex", a);
            _activeSh.SetBuffer(_activeIds.MulEach, "_BufComplex2", b);
            _activeSh.Dispatch(_activeIds.MulEach, Groups1D, 1, 1);
        }

        private void FFTNormBuffer(ComputeBuffer buf)
        {
            _activeSh.SetBuffer(_activeIds.FFTNorm, "_BufComplex", buf);
            _activeSh.Dispatch(_activeIds.FFTNorm, Groups1D, 1, 1);
        }

        /// <summary>
        /// Скорость из фазовых градиентов ψ. Uniform <c>_Hbar</c> в kernel задаётся как <b>1</b> (как до ветки χ/ёмкости);
        /// подставка ℏ ломала масштаб u, проекцию давления и трассеры в Jet/кольцах.
        /// </summary>
        private void VelocityOneForm(CSVelocity v)
        {
            VelocityOneForm(v, 1.0f);
        }

        private void VelocityOneForm(CSVelocity v, float h)
        {
            SetCommonUniforms();
            _activeSh.SetFloat("_Hbar", h);
            _activeSh.SetBuffer(_activeIds.VelOne, "_Psi1", psi1);
            _activeSh.SetBuffer(_activeIds.VelOne, "_Psi2", psi2);
            _activeSh.SetBuffer(_activeIds.VelOne, "_VX", v.vx);
            _activeSh.SetBuffer(_activeIds.VelOne, "_VY", v.vy);
            _activeSh.SetBuffer(_activeIds.VelOne, "_VZ", v.vz);
            _activeSh.Dispatch(_activeIds.VelOne, Groups1D, 1, 1);
        }

        private void StaggeredSharp(CSVelocity vel)
        {
            SetCommonUniforms();

            _activeSh.SetBuffer(_activeIds.Staggered, "_BufFloat", vel.vx);
            _activeSh.SetFloat("_Param", dx);
            _activeSh.Dispatch(_activeIds.Staggered, Groups1D, 1, 1);

            _activeSh.SetBuffer(_activeIds.Staggered, "_BufFloat", vel.vy);
            _activeSh.SetFloat("_Param", dy);
            _activeSh.Dispatch(_activeIds.Staggered, Groups1D, 1, 1);

            _activeSh.SetBuffer(_activeIds.Staggered, "_BufFloat", vel.vz);
            _activeSh.SetFloat("_Param", dz);
            _activeSh.Dispatch(_activeIds.Staggered, Groups1D, 1, 1);
        }

        private void Div(CSVelocity v, ComputeBuffer result)
        {
            SetCommonUniforms();
            _activeSh.SetBuffer(_activeIds.Div, "_VX", v.vx);
            _activeSh.SetBuffer(_activeIds.Div, "_VY", v.vy);
            _activeSh.SetBuffer(_activeIds.Div, "_VZ", v.vz);
            _activeSh.SetBuffer(_activeIds.Div, "_BufFloat", result);
            _activeSh.Dispatch(_activeIds.Div, Groups1D, 1, 1);
        }

        private void PoissonSolve(ComputeBuffer f, ComputeBuffer result)
        {
            SetCommonUniforms();

            _activeSh.SetBuffer(_activeIds.CopyR2C, "_BufComplex", result);
            _activeSh.SetBuffer(_activeIds.CopyR2C, "_BufFloat", f);
            _activeSh.Dispatch(_activeIds.CopyR2C, Groups1D, 1, 1);

            _fft.FFT3D(result, false);
            MulEachBuffers(result, _fac);
            _fft.FFT3D(result, true);
        }

        public void PressureProject(CSVelocity v)
        {
            Div(v, _divResult);
            PoissonSolve(_divResult, _poissonTemp);

            SetCommonUniforms();
            _activeSh.SetBuffer(_activeIds.Gauge, "_Psi1", psi1);
            _activeSh.SetBuffer(_activeIds.Gauge, "_Psi2", psi2);
            _activeSh.SetBuffer(_activeIds.Gauge, "_BufComplex", _poissonTemp);
            _activeSh.Dispatch(_activeIds.Gauge, Groups1D, 1, 1);
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
            if (ReferenceEquals(_activeSh, _kExt) && _kExt != null)
                _activeSh.SetInt("_GravityBothPsi", rotatePsi1Too ? 1 : 0);
            _activeSh.SetFloat("_GX", g.x);
            _activeSh.SetFloat("_GY", g.y);
            _activeSh.SetFloat("_GZ", g.z);
            _activeSh.SetFloat("_DT", dt);
            _activeSh.SetBuffer(_activeIds.Grav, "_Psi1", psi1);
            _activeSh.SetBuffer(_activeIds.Grav, "_Psi2", psi2);
            _activeSh.SetBuffer(_activeIds.Grav, "_PX", _px);
            _activeSh.SetBuffer(_activeIds.Grav, "_PY", _py);
            _activeSh.SetBuffer(_activeIds.Grav, "_PZ", _pz);
            _activeSh.Dispatch(_activeIds.Grav, Groups1D, 1, 1);
        }

        private void ApplyHeatSinkPsi1(ComputeBuffer isJet)
        {
            SetCommonUniforms();
            _activeSh.SetBuffer(_activeIds.Heat, "_Psi1", psi1);
            _activeSh.SetBuffer(_activeIds.Heat, "_IsJet", isJet);
            _activeSh.Dispatch(_activeIds.Heat, Groups1D, 1, 1);
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
            if (_kExt == null || _extIds.AddGravVel < 0)
            {
                Debug.LogError("[CSISF] AddGravityToVelocity только в SFComputeKernelsContainer — задай второй ISF compute для сценария ёмкости.");
                return;
            }
            BindExtKernelUniforms();
            _kExt.SetFloat("_GX", g.x * dt);
            _kExt.SetFloat("_GY", g.y * dt);
            _kExt.SetFloat("_GZ", g.z * dt);
            _kExt.SetBuffer(_extIds.AddGravVel, "_VX", vel.vx);
            _kExt.SetBuffer(_extIds.AddGravVel, "_VY", vel.vy);
            _kExt.SetBuffer(_extIds.AddGravVel, "_VZ", vel.vz);
            _kExt.Dispatch(_extIds.AddGravVel, Groups1D, 1, 1);
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

        /// <summary>
        /// Плавучесть (Буссинеск): u += β·χ·dir·dt по полю скорости, затем <see cref="PressureProject(CSVelocity)"/>.
        /// Сила направлена вдоль <paramref name="dir"/> (обычно +Y) и пропорциональна концентрации дыма χ.
        /// Требует <see cref="useLiquidChiField"/> (χ — носитель дыма) и расширенный compute (ёмкость/комната).
        /// </summary>
        public void ApplyBuoyancy(Vector3 dir, float beta, CSVelocity vel)
        {
            if (!useLiquidChiField || _kExt == null || _extIds.Buoyancy < 0)
            {
                Debug.LogError("[CSISF] AddBuoyancyToVelocity только в SFComputeKernelsContainer при включённом χ-поле.");
                return;
            }
            Vector3 a = dir.normalized * (beta * dt);
            BindExtKernelUniforms();
            _kExt.SetFloat("_BuoyAX", a.x);
            _kExt.SetFloat("_BuoyAY", a.y);
            _kExt.SetFloat("_BuoyAZ", a.z);
            _kExt.SetBuffer(_extIds.Buoyancy, "_ChiField", _liquidChi);
            _kExt.SetBuffer(_extIds.Buoyancy, "_VX", vel.vx);
            _kExt.SetBuffer(_extIds.Buoyancy, "_VY", vel.vy);
            _kExt.SetBuffer(_extIds.Buoyancy, "_VZ", vel.vz);
            _kExt.Dispatch(_extIds.Buoyancy, Groups1D, 1, 1);
            // Горизонтальный снос убираем, вертикаль (плавучесть) — нет.
            RemoveMeanHorizontalVelocity(vel);
            PressureProject(vel);
            ApplyGasVacuumFromChi();
        }

        /// <summary>Источник дыма: χ → max(χ, value) в ячейках маски. Требует <see cref="useLiquidChiField"/>.</summary>
        public void InjectChi(ComputeBuffer mask, float value)
        {
            if (!useLiquidChiField || _kExt == null || _extIds.InjectChi < 0)
                return;
            BindExtKernelUniforms();
            _kExt.SetFloat("_ChiInjectValue", value);
            _kExt.SetBuffer(_extIds.InjectChi, "_ChiRW", _liquidChi);
            _kExt.SetBuffer(_extIds.InjectChi, "_IsJet", mask);
            _kExt.Dispatch(_extIds.InjectChi, Groups1D, 1, 1);
        }

        /// <summary>Вытяжка: χ *= decay в ячейках маски (сток дыма). Требует <see cref="useLiquidChiField"/>.</summary>
        public void VentChiSink(ComputeBuffer mask, float decay)
        {
            if (!useLiquidChiField || _kExt == null || _extIds.ChiSink < 0)
                return;
            BindExtKernelUniforms();
            _kExt.SetFloat("_ChiVentDecay", Mathf.Clamp01(decay));
            _kExt.SetBuffer(_extIds.ChiSink, "_ChiRW", _liquidChi);
            _kExt.SetBuffer(_extIds.ChiSink, "_IsJet", mask);
            _kExt.Dispatch(_extIds.ChiSink, Groups1D, 1, 1);
        }

        /// <summary>Турбулентная диффузия χ: расширяет султан и даёт растекание дыма. alpha = D·dt/h² ∈ [0, ~0.16].</summary>
        public void DiffuseChi(float alpha)
        {
            if (!useLiquidChiField || _kExt == null || _extIds.DiffuseChi < 0 || alpha <= 0f)
                return;
            BindExtKernelUniforms();
            _kExt.SetFloat("_ChiDiffAlpha", Mathf.Clamp(alpha, 0f, 0.16f));
            _kExt.SetBuffer(_extIds.DiffuseChi, "_ChiIn", _liquidChi);
            _kExt.SetBuffer(_extIds.DiffuseChi, "_ChiOut", _liquidChiTemp);
            _kExt.Dispatch(_extIds.DiffuseChi, Groups1D, 1, 1);
            var swap = _liquidChi;
            _liquidChi = _liquidChiTemp;
            _liquidChiTemp = swap;
        }

        /// <summary>Curl-noise турбулентность: добавляет к скорости бездивергентное вихревое поле внутри дыма (клубление, боковое вовлечение).</summary>
        public void AddCurlTurbulence(CSVelocity vel, float amp, float scale, Vector3 timeOffset, float chiLo, float chiHi)
        {
            if (!useLiquidChiField || _kExt == null || _extIds.CurlNoise < 0 || amp <= 0f)
                return;
            BindExtKernelUniforms();
            _kExt.SetFloat("_CurlAmp", amp);
            _kExt.SetFloat("_CurlScale", Mathf.Max(scale, 1e-4f));
            _kExt.SetVector("_CurlTimeOffset", timeOffset);
            _kExt.SetFloat("_CurlChiLo", chiLo);
            _kExt.SetFloat("_CurlChiHi", chiHi);
            _kExt.SetBuffer(_extIds.CurlNoise, "_ChiField", _liquidChi);
            _kExt.SetBuffer(_extIds.CurlNoise, "_VX", vel.vx);
            _kExt.SetBuffer(_extIds.CurlNoise, "_VY", vel.vy);
            _kExt.SetBuffer(_extIds.CurlNoise, "_VZ", vel.vz);
            _kExt.Dispatch(_extIds.CurlNoise, Groups1D, 1, 1);
        }

        public void ApplyJetBoundary(ComputeBuffer isJet,
            float kvecX, float kvecY, float kvecZ, float phaseOffset)
        {
            SetCommonUniforms();
            _activeSh.SetFloat("_KVecX", kvecX);
            _activeSh.SetFloat("_KVecY", kvecY);
            _activeSh.SetFloat("_KVecZ", kvecZ);
            _activeSh.SetFloat("_PhaseOffset", phaseOffset);

            _activeSh.SetBuffer(_activeIds.Jet, "_Psi1", psi1);
            _activeSh.SetBuffer(_activeIds.Jet, "_Psi2", psi2);
            _activeSh.SetBuffer(_activeIds.Jet, "_IsJet", isJet);
            _activeSh.SetBuffer(_activeIds.Jet, "_PX", _px);
            _activeSh.SetBuffer(_activeIds.Jet, "_PY", _py);
            _activeSh.SetBuffer(_activeIds.Jet, "_PZ", _pz);
            _activeSh.Dispatch(_activeIds.Jet, Groups1D, 1, 1);
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
