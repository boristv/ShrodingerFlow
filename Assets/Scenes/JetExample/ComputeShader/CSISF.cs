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
        private int _gravK, _heatK, _uniformForceK;

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
            _uniformForceK = kernels.FindKernel("UniformForce");

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
        }

        public void Normalize()
        {
            SetCommonUniforms();
            _kernels.SetBuffer(_normalizeK, "_Psi1", psi1);
            _kernels.SetBuffer(_normalizeK, "_Psi2", psi2);
            _kernels.Dispatch(_normalizeK, Groups1D, 1, 1);
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

        private void VelocityOneForm(CSVelocity v)
        {
            VelocityOneForm(v, 1.0f);
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
            PressureProject(_velTemp);
        }

        private void LESStep()
        {
            VelocityOneForm(_velCurrent);
            _les.Filter(_velCurrent, _velFiltered, dx, dy, dz);
            _les.ComputeNuT(_velFiltered, dx, dy, dz, Cs, filterFac);
            _les.ApplyViscosity(_velCurrent, dx, dy, dz, dt, kinematicViscosity, 1f);
        }

        private void LaminarViscosityStep()
        {
            if (kinematicViscosity <= 1e-20f)
                return;
            VelocityOneForm(_velCurrent);
            _les.ApplyViscosity(_velCurrent, dx, dy, dz, dt, kinematicViscosity, 0f);
        }

        private bool UsesVelocityFieldViscosity =>
            kinematicViscosity > 1e-20f;

        /// <param name="gravityPsi2">
        /// Опционально: после нормировки применить к ψ₂ фазу exp(i (g·x) dt / ℏ) — тот же шаг, что в
        /// <see cref="UpdateCigaretteSpace"/> (гравитация / «сила» через потенциал, а не через v).
        /// </param>
        public void UpdateSpace(bool useLES = true, Vector3? gravityPsi2 = null)
        {
            SetCommonUniforms();
            SchoedingerFlow();
            if (useLES)
                LESStep();
            else
                LaminarViscosityStep();
            Normalize();
            if (gravityPsi2.HasValue && gravityPsi2.Value.sqrMagnitude > 1e-20f)
                ApplyGravityPsi2(gravityPsi2.Value);
            bool ppFromVel = useLES || UsesVelocityFieldViscosity;
            if (ppFromVel)
                PressureProject(_velCurrent);
            else
                PressureProject();
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
            ApplyGravityPsi2(gravityG);
            bool ppFromVel = useLES || UsesVelocityFieldViscosity;
            if (ppFromVel)
                PressureProject(_velCurrent);
            else
                PressureProject();
            ApplyHeatSinkPsi1(isJetMask);
            Normalize();
            if (ppFromVel)
                PressureProject(_velCurrent);
            else
                PressureProject();
        }

        private void ApplyGravityPsi2(Vector3 g)
        {
            SetCommonUniforms();
            _kernels.SetFloat("_GX", g.x);
            _kernels.SetFloat("_GY", g.y);
            _kernels.SetFloat("_GZ", g.z);
            _kernels.SetFloat("_DT", dt);
            _kernels.SetBuffer(_gravK, "_Psi2", psi2);
            _kernels.SetBuffer(_gravK, "_PX", _px);
            _kernels.SetBuffer(_gravK, "_PY", _py);
            _kernels.SetBuffer(_gravK, "_PZ", _pz);
            _kernels.Dispatch(_gravK, Groups1D, 1, 1);
        }

        /// <summary>
        /// Однородная объёмная сила (постоянный градиент давления / «ветер») на обе компоненты ψ.
        /// Прибавляет скорость f·dt всей жидкости; после вызова нужен PressureProject.
        /// Для SmokeMaze2D: гонит поток источник → лабиринт → вытяжка (локальный outflow этого не делает).
        /// </summary>
        public void ApplyUniformForce(Vector3 f)
        {
            SetCommonUniforms();
            _kernels.SetFloat("_GX", f.x);
            _kernels.SetFloat("_GY", f.y);
            _kernels.SetFloat("_GZ", f.z);
            _kernels.SetFloat("_DT", dt);
            _kernels.SetBuffer(_uniformForceK, "_Psi1", psi1);
            _kernels.SetBuffer(_uniformForceK, "_Psi2", psi2);
            _kernels.SetBuffer(_uniformForceK, "_PX", _px);
            _kernels.SetBuffer(_uniformForceK, "_PY", _py);
            _kernels.SetBuffer(_uniformForceK, "_PZ", _pz);
            _kernels.Dispatch(_uniformForceK, Groups1D, 1, 1);
        }

        private void ApplyHeatSinkPsi1(ComputeBuffer isJet)
        {
            SetCommonUniforms();
            _kernels.SetBuffer(_heatK, "_Psi1", psi1);
            _kernels.SetBuffer(_heatK, "_IsJet", isJet);
            _kernels.Dispatch(_heatK, Groups1D, 1, 1);
        }

        public void UpdateVelocities(CSVelocity vel)
        {
            VelocityOneForm(vel, hbar);
            StaggeredSharp(vel);
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
            _les?.Dispose();
        }
    }
}
