using UnityEngine;

namespace ComputeShaderSF
{
    public class CSLES : System.IDisposable
    {
        private ComputeShader _shader;
        private int _boxFilterKernel;
        private int _computeStrainKernel;
        private int _addTurbViscKernel;

        private ComputeBuffer _nuT;
        private int _nx, _ny, _nz;

        public void Init(ComputeShader shader, int nx, int ny, int nz)
        {
            _shader = shader;
            _nx = nx;
            _ny = ny;
            _nz = nz;

            _nuT = new ComputeBuffer(nx * ny * nz, sizeof(float));

            _boxFilterKernel = shader.FindKernel("BoxFilter");
            _computeStrainKernel = shader.FindKernel("ComputeStrain");
            _addTurbViscKernel = shader.FindKernel("AddTurbViscosity");
        }

        private void SetGridUniforms()
        {
            _shader.SetInt("_NX", _nx);
            _shader.SetInt("_NY", _ny);
            _shader.SetInt("_NZ", _nz);
        }

        private void Dispatch3D(int kernel)
        {
            _shader.Dispatch(kernel,
                (_nx + 7) / 8, (_ny + 7) / 8, (_nz + 7) / 8);
        }

        public void Filter(CSVelocity src, CSVelocity dst,
            float dx, float dy, float dz)
        {
            SetGridUniforms();

            _shader.SetBuffer(_boxFilterKernel, "_Dst", dst.vx);
            _shader.SetBuffer(_boxFilterKernel, "_Src", src.vx);
            Dispatch3D(_boxFilterKernel);

            _shader.SetBuffer(_boxFilterKernel, "_Dst", dst.vy);
            _shader.SetBuffer(_boxFilterKernel, "_Src", src.vy);
            Dispatch3D(_boxFilterKernel);

            _shader.SetBuffer(_boxFilterKernel, "_Dst", dst.vz);
            _shader.SetBuffer(_boxFilterKernel, "_Src", src.vz);
            Dispatch3D(_boxFilterKernel);
        }

        public void ComputeNuT(CSVelocity filtered,
            float dx, float dy, float dz, float Cs, float filterFac)
        {
            SetGridUniforms();
            float delta = filterFac * Mathf.Min(dx, Mathf.Min(dy, dz));

            _shader.SetFloat("_DX", dx);
            _shader.SetFloat("_DY", dy);
            _shader.SetFloat("_DZ", dz);
            _shader.SetFloat("_Cs", Cs);
            _shader.SetFloat("_Delta", delta);

            _shader.SetBuffer(_computeStrainKernel, "_NuT", _nuT);
            _shader.SetBuffer(_computeStrainKernel, "_U", filtered.vx);
            _shader.SetBuffer(_computeStrainKernel, "_V", filtered.vy);
            _shader.SetBuffer(_computeStrainKernel, "_W", filtered.vz);
            Dispatch3D(_computeStrainKernel);
        }

        /// <param name="nuMol">Молекулярная (ламинарная) кинематическая вязкость, добавляется к ν.</param>
        /// <param name="turbNuScale">1 — использовать ν_t из буфера (Smagorinsky); 0 — только νMol.</param>
        public void ApplyViscosity(CSVelocity vel,
            float dx, float dy, float dz, float dt, float nuMol, float turbNuScale)
        {
            SetGridUniforms();

            _shader.SetFloat("_DX", dx);
            _shader.SetFloat("_DY", dy);
            _shader.SetFloat("_DZ", dz);
            _shader.SetFloat("_DT", dt);
            _shader.SetFloat("_NuMol", nuMol);
            _shader.SetFloat("_TurbNuScale", turbNuScale);

            _shader.SetBuffer(_addTurbViscKernel, "_NuT", _nuT);
            _shader.SetBuffer(_addTurbViscKernel, "_U", vel.vx);
            _shader.SetBuffer(_addTurbViscKernel, "_V", vel.vy);
            _shader.SetBuffer(_addTurbViscKernel, "_W", vel.vz);
            Dispatch3D(_addTurbViscKernel);
        }

        public void Dispose()
        {
            _nuT?.Release();
        }
    }
}
