using System;
using UnityEngine;

namespace ComputeShaderSF
{
    /// <summary>
    /// Vorticity confinement (Fedkiw) на транспортной скорости ũ: f = ε·(N×ω), N = ∇|ω|/|∇|ω||.
    /// Возвращает мелкие завихрения, размытые численной диффузией переноса → клубящийся факел.
    /// Действует только на скорость, несущую α/T/трассеры; ψ-динамика ISF не затрагивается.
    /// </summary>
    public class CSVorticityConfine : IDisposable
    {
        private readonly int rx, ry, rz, num;
        private readonly float dx, dy, dz;
        private readonly ComputeShader _sh;
        private readonly int _magK, _confK, _copyK;
        private ComputeBuffer _w, _ox, _oy, _oz;

        public CSVorticityConfine(ComputeShader shader, int rx, int ry, int rz, float dx, float dy, float dz)
        {
            _sh = shader;
            this.rx = rx; this.ry = ry; this.rz = rz; num = rx * ry * rz;
            this.dx = dx; this.dy = dy; this.dz = dz;
            _magK = shader.FindKernel("VorticityMag");
            _confK = shader.FindKernel("VorticityConfine");
            _copyK = shader.FindKernel("CopyVel");
            _w = new ComputeBuffer(num, sizeof(float));
            _ox = new ComputeBuffer(num, sizeof(float));
            _oy = new ComputeBuffer(num, sizeof(float));
            _oz = new ComputeBuffer(num, sizeof(float));
        }

        private int Groups => (num + 255) / 256;

        private void Bind(int k, CSVelocity vel)
        {
            _sh.SetBuffer(k, "_VelX", vel.vx);
            _sh.SetBuffer(k, "_VelY", vel.vy);
            _sh.SetBuffer(k, "_VelZ", vel.vz);
            _sh.SetBuffer(k, "_Wmag", _w);
            _sh.SetBuffer(k, "_OX", _ox);
            _sh.SetBuffer(k, "_OY", _oy);
            _sh.SetBuffer(k, "_OZ", _oz);
        }

        public void Apply(CSVelocity vel, float strength, float dt)
        {
            if (strength <= 0f) return;

            _sh.SetInt("_ResX", rx); _sh.SetInt("_ResY", ry); _sh.SetInt("_ResZ", rz); _sh.SetInt("_Num", num);
            _sh.SetFloat("_DX", dx); _sh.SetFloat("_DY", dy); _sh.SetFloat("_DZ", dz);
            _sh.SetFloat("_DT", dt);
            _sh.SetFloat("_Eps", strength * Mathf.Min(dx, Mathf.Min(dy, dz)));

            Bind(_magK, vel); _sh.Dispatch(_magK, Groups, 1, 1);
            Bind(_confK, vel); _sh.Dispatch(_confK, Groups, 1, 1);
            Bind(_copyK, vel); _sh.Dispatch(_copyK, Groups, 1, 1);
        }

        public void Dispose()
        {
            _w?.Release(); _ox?.Release(); _oy?.Release(); _oz?.Release();
        }
    }
}
