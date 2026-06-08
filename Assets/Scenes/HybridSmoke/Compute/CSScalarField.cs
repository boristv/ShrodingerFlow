using System;
using UnityEngine;

namespace ComputeShaderSF
{
    /// <summary>
    /// Плотностно-фазовое поле α гибридной волново-плотностной модели (диссертация, гл. 3.1, 5.1.5).
    /// Переносится скоростью ũ, восстановленной из ψ и стабилизированной LES:
    /// ∂α/∂t + ũ·∇α = D_α∇²α. α — носитель плотности дыма / интерфейса / оптики (поле состояния, не частицы).
    /// </summary>
    public class CSScalarField : IDisposable
    {
        public int resX, resY, resZ, num;
        public float dx, dy, dz;

        /// <summary>Текущее поле α (актуальный буфер после Step).</summary>
        public ComputeBuffer Alpha => _a;

        private ComputeBuffer _a, _b, _c;
        private readonly ComputeShader _shader;
        private readonly int _clearK, _advectK, _diffuseK, _boundaryK, _buoyK, _mcK;

        public CSScalarField(ComputeShader shader, int rx, int ry, int rz, float dx, float dy, float dz)
        {
            _shader = shader;
            resX = rx; resY = ry; resZ = rz; num = rx * ry * rz;
            this.dx = dx; this.dy = dy; this.dz = dz;

            _clearK = shader.FindKernel("ClearScalar");
            _advectK = shader.FindKernel("AdvectScalar");
            _diffuseK = shader.FindKernel("DiffuseScalar");
            _boundaryK = shader.FindKernel("ScalarBoundary");
            _buoyK = shader.FindKernel("BuoyancyPotential");
            _mcK = shader.FindKernel("MacCormackCombine");

            _a = new ComputeBuffer(num, sizeof(float));
            _b = new ComputeBuffer(num, sizeof(float));
            _c = new ComputeBuffer(num, sizeof(float));
            Clear();
        }

        private int Groups => (num + 255) / 256;

        private void SetGrid()
        {
            _shader.SetInt("_ResX", resX);
            _shader.SetInt("_ResY", resY);
            _shader.SetInt("_ResZ", resZ);
            _shader.SetInt("_Num", num);
            _shader.SetFloat("_DX", dx);
            _shader.SetFloat("_DY", dy);
            _shader.SetFloat("_DZ", dz);
        }

        public void Clear()
        {
            SetGrid();
            _shader.SetBuffer(_clearK, "_Alpha", _a);
            _shader.Dispatch(_clearK, Groups, 1, 1);
            _shader.SetBuffer(_clearK, "_Alpha", _b);
            _shader.Dispatch(_clearK, Groups, 1, 1);
        }

        /// <summary>
        /// Шаг переноса α: advection (semi-Lagrangian) → опц. диффузия → граничные условия по единым маскам.
        /// Скорость <paramref name="vel"/> — стабилизированная ũ из ISF+LES (CSISF.UpdateVelocities).
        /// </summary>
        public void Step(CSVelocity vel, ComputeBuffer solidMask, ComputeBuffer sourceMask, ComputeBuffer sinkMask,
            float dt, float diffusion, float sourceValue, float sinkFactor, int diffuseIters = 0, float decay = 0f,
            bool macCormack = false)
        {
            SetGrid();
            _shader.SetFloat("_DT", dt);

            if (!macCormack)
            {
                // Полулагранжев перенос: _a -> _b
                Advect(_a, _b, vel, solidMask);
                Swap();
            }
            else
            {
                // MacCormack: forward _a->_b, backward _b->_c (−dt), коррекция в _a (ин-плейс).
                Advect(_a, _b, vel, solidMask);
                _shader.SetFloat("_DT", -dt);
                Advect(_b, _c, vel, solidMask);
                _shader.SetFloat("_DT", dt);
                _shader.SetBuffer(_mcK, "_Alpha", _a);
                _shader.SetBuffer(_mcK, "_AlphaOut", _b);
                _shader.SetBuffer(_mcK, "_MBack", _c);
                _shader.SetBuffer(_mcK, "_SolidMask", solidMask);
                _shader.Dispatch(_mcK, Groups, 1, 1);
            }

            // Diffusion: несколько явных шагов (по устойчивости D·dt/dx² ≤ 1/6)
            if (diffusion > 0f && diffuseIters > 0)
            {
                _shader.SetFloat("_DAlpha", diffusion);
                for (int it = 0; it < diffuseIters; it++)
                {
                    _shader.SetBuffer(_diffuseK, "_AlphaIn", _a);
                    _shader.SetBuffer(_diffuseK, "_AlphaOut", _b);
                    _shader.SetBuffer(_diffuseK, "_SolidMask", solidMask);
                    _shader.Dispatch(_diffuseK, Groups, 1, 1);
                    Swap();
                }
            }

            // Boundary (in-place на актуальном _a)
            _shader.SetFloat("_SourceValue", sourceValue);
            _shader.SetFloat("_SinkFactor", sinkFactor);
            _shader.SetFloat("_Decay", decay);
            _shader.SetBuffer(_boundaryK, "_Alpha", _a);
            _shader.SetBuffer(_boundaryK, "_SolidMask", solidMask);
            _shader.SetBuffer(_boundaryK, "_SourceMask", sourceMask);
            _shader.SetBuffer(_boundaryK, "_SinkMask", sinkMask);
            _shader.Dispatch(_boundaryK, Groups, 1, 1);
        }

        /// <summary>Потенциал плавучести: вертикальный (по Y) интеграл текущего α в <paramref name="outB"/>.</summary>
        public void ComputeBuoyancyPotential(ComputeBuffer outB)
        {
            SetGrid();
            _shader.SetBuffer(_buoyK, "_Alpha", _a);
            _shader.SetBuffer(_buoyK, "_BuoyOut", outB);
            int cols = resX * resZ;
            _shader.Dispatch(_buoyK, (cols + 63) / 64, 1, 1);
        }

        private void Advect(ComputeBuffer src, ComputeBuffer dst, CSVelocity vel, ComputeBuffer solidMask)
        {
            _shader.SetBuffer(_advectK, "_AlphaIn", src);
            _shader.SetBuffer(_advectK, "_AlphaOut", dst);
            _shader.SetBuffer(_advectK, "_VX", vel.vx);
            _shader.SetBuffer(_advectK, "_VY", vel.vy);
            _shader.SetBuffer(_advectK, "_VZ", vel.vz);
            _shader.SetBuffer(_advectK, "_SolidMask", solidMask);
            _shader.Dispatch(_advectK, Groups, 1, 1);
        }

        private void Swap()
        {
            (_a, _b) = (_b, _a);
        }

        public void Dispose()
        {
            _a?.Release();
            _b?.Release();
            _c?.Release();
        }
    }
}
