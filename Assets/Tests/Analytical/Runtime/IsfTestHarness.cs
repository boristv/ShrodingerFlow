using System;
using ComputeShaderSF;
using UnityEngine;

namespace ShrodingerFlow.AnalyticalTests
{
    /// <summary>
    /// Декларативное описание вихревого кольца: используется <see cref="IsfTestHarness.InitWithVortexRings"/>.
    /// Соответствует ядру <c>AddCircle</c> из <c>SFUnifiedCS</c>.
    /// </summary>
    [Serializable]
    public struct VortexRingDef
    {
        [Tooltip("Центр кольца в координатах объёма [0..vol_size].")]
        public Vector3 Center;
        [Tooltip("Нормаль (направление переноса вихря: фактическое движение — против +n).")]
        public Vector3 Normal;
        [Tooltip("Радиус кольца (расстояние от оси до сечения трубки).")]
        public float Radius;
        [Tooltip("Толщина переходного слоя d. Если <= 0 — будет взято 5·dx.")]
        public float Thickness;

        public static VortexRingDef Default(Vector3 center, Vector3 normal, float radius) =>
            new VortexRingDef { Center = center, Normal = normal, Radius = radius, Thickness = 0f };
    }

    /// <summary>Параметры одного шага симуляции.</summary>
    public struct StepOptions
    {
        public bool UseLES;
        public float KinematicViscosity;
        public Vector3? GravityPsi2;
        public float CsSmagorinsky;

        public static StepOptions Default => new StepOptions
        {
            UseLES = false,
            KinematicViscosity = 0f,
            GravityPsi2 = null,
            CsSmagorinsky = 0.07f
        };
    }

    /// <summary>
    /// Тонкая обёртка над <see cref="CSISF"/> и <see cref="CSVelocity"/> для аналитических тестов.
    /// Не зависит от <c>SFUnifiedCS</c>, частиц и рендера — только физика на сетке.
    /// </summary>
    public sealed class IsfTestHarness : IDisposable
    {
        private readonly CSISF _isf;
        private readonly CSVelocity _velocity;
        private readonly int[] _volSize;
        private readonly int[] _volRes;
        private float _simTime;

        public CSISF Isf => _isf;
        public CSVelocity Velocity => _velocity;
        public int[] VolSize => _volSize;
        public int[] VolRes => _volRes;
        public float Dt => _isf.dt;
        public float Hbar => _isf.hbar;
        public float SimTime => _simTime;

        /// <summary>
        /// Инициализация ISF с сеткой и параметрами. <paramref name="volRes"/> при необходимости будет
        /// приведено к ближайшим степеням двойки <see cref="CSISF.Init"/> — это поведение CSISF.
        /// </summary>
        public IsfTestHarness(ComputeShader kernels, ComputeShader fft, ComputeShader les,
            int[] volSize, int[] volRes, float hbar, float dt)
        {
            _volSize = (int[])volSize.Clone();
            _volRes = (int[])volRes.Clone();
            _isf = new CSISF();
            _isf.Init(kernels, fft, les, _volSize, _volRes, hbar, dt);
            _velocity = new CSVelocity(_isf.resX, _isf.resY, _isf.resZ);
        }

        /// <summary>ψ₁ = 1+0i, ψ₂ = (psi2Initial.x + i·psi2Initial.y) — слабая «вторая компонента» как в исходных тестах.</summary>
        public void InitUniformPsi(Vector2? psi2Initial = null)
        {
            int num = _isf.num;
            var tmp1 = new Vector2[num];
            var tmp2 = new Vector2[num];
            Vector2 p2 = psi2Initial ?? new Vector2(0.01f, 0f);
            for (int i = 0; i < num; i++)
            {
                tmp1[i] = new Vector2(1f, 0f);
                tmp2[i] = p2;
            }
            _isf.psi1.SetData(tmp1);
            _isf.psi2.SetData(tmp2);
            _isf.Normalize();
        }

        /// <summary>
        /// Гауссов волновой пакет: ψ₁ = A·exp(-|x-c|²/(2σ²))·exp(i k·x), k = v/ℏ.
        /// Плотность ρ = |ψ₁|² — локализованное Гауссово пятно; при <paramref name="groupVelocity"/> ≠ 0
        /// центр масс пятна движется с этой скоростью.
        /// </summary>
        public void InitGaussianBlob(Vector3 center, float sigma, Vector3 groupVelocity, float psi2Amplitude = 0.01f)
        {
            int num = _isf.num;
            float kx = groupVelocity.x / _isf.hbar;
            float ky = groupVelocity.y / _isf.hbar;
            float kz = groupVelocity.z / _isf.hbar;
            float invTwoSigSq = sigma > 0f ? 1f / (2f * sigma * sigma) : 1e6f;

            var tmp1 = new Vector2[num];
            var tmp2 = new Vector2[num];
            for (int i = 0; i < num; i++)
            {
                float ex = _isf.pxCPU[i] - center.x;
                float ey = _isf.pyCPU[i] - center.y;
                float ez = _isf.pzCPU[i] - center.z;
                float envelope = Mathf.Exp(-(ex * ex + ey * ey + ez * ez) * invTwoSigSq);
                float phase = kx * _isf.pxCPU[i] + ky * _isf.pyCPU[i] + kz * _isf.pzCPU[i];
                float c = Mathf.Cos(phase), s = Mathf.Sin(phase);
                tmp1[i] = new Vector2(envelope * c, envelope * s);
                tmp2[i] = new Vector2(psi2Amplitude * c, psi2Amplitude * s);
            }
            _isf.psi1.SetData(tmp1);
            _isf.psi2.SetData(tmp2);
            _isf.Normalize();
        }

        /// <summary>Плоская волна ψ = exp(i k·x), k = U/ℏ (как Set_background_flow в example_*.hip).</summary>
        public void InitPlaneWave(Vector3 backgroundU, float psi2Amplitude = 0.01f)
        {
            int num = _isf.num;
            float kx = backgroundU.x / _isf.hbar;
            float ky = backgroundU.y / _isf.hbar;
            float kz = backgroundU.z / _isf.hbar;

            var tmp1 = new Vector2[num];
            var tmp2 = new Vector2[num];
            for (int i = 0; i < num; i++)
            {
                float phase = kx * _isf.pxCPU[i] + ky * _isf.pyCPU[i] + kz * _isf.pzCPU[i];
                float c = Mathf.Cos(phase), s = Mathf.Sin(phase);
                tmp1[i] = new Vector2(c, s);
                tmp2[i] = new Vector2(c * psi2Amplitude, s * psi2Amplitude);
            }
            _isf.psi1.SetData(tmp1);
            _isf.psi2.SetData(tmp2);
            _isf.Normalize();
        }

        /// <summary>Добавить N вихревых колец на текущий ψ₁ (как AddCircle в SFUnifiedCS).</summary>
        public void AddVortexRings(params VortexRingDef[] rings)
        {
            if (rings == null || rings.Length == 0)
                return;

            int num = _isf.num;
            var psi = new Vector2[num];
            _isf.psi1.GetData(psi);
            float defaultD = _isf.dx * 5f;

            foreach (var ring in rings)
            {
                float d = ring.Thickness > 0f ? ring.Thickness : defaultD;
                AddCircleInPlace(psi, ring.Center, ring.Normal, ring.Radius, d);
            }

            _isf.psi1.SetData(psi);
        }

        private void AddCircleInPlace(Vector2[] psi, Vector3 center, Vector3 normal, float radius, float d)
        {
            normal = normal.normalized;
            int rx = _isf.resX, ry = _isf.resY, rz = _isf.resZ;
            float r2 = radius * radius;

            for (int i = 0; i < rx; i++)
            {
                for (int j = 0; j < ry; j++)
                {
                    for (int k = 0; k < rz; k++)
                    {
                        int idx = i * ry * rz + j * rz + k;
                        float ex = _isf.pxCPU[idx] - center.x;
                        float ey = _isf.pyCPU[idx] - center.y;
                        float ez = _isf.pzCPU[idx] - center.z;

                        float z = ex * normal.x + ey * normal.y + ez * normal.z;
                        float rPerp2 = ex * ex + ey * ey + ez * ez - z * z;

                        float alpha = 0f;
                        if (rPerp2 < r2)
                        {
                            if (z > 0f && z <= d * 0.5f)
                                alpha = -Mathf.PI * (2f * z / d - 1f);
                            else if (z <= 0f && z >= -d * 0.5f)
                                alpha = -Mathf.PI * (2f * z / d + 1f);
                        }

                        if (alpha != 0f)
                        {
                            float ca = Mathf.Cos(alpha), sa = Mathf.Sin(alpha);
                            float pr = psi[idx].x, pi = psi[idx].y;
                            psi[idx] = new Vector2(pr * ca - pi * sa, pr * sa + pi * ca);
                        }
                    }
                }
            }
        }

        /// <summary>После закладки полей — обязательная пара Normalize + PressureProject (зануляет дивергенцию).</summary>
        public void NormalizeAndProject()
        {
            _isf.Normalize();
            _isf.PressureProject();
        }

        /// <summary>Один шаг ISF: <see cref="CSISF.UpdateSpace"/> + обновление скоростного поля для диагностики.</summary>
        public void Step(StepOptions options)
        {
            _isf.Cs = options.CsSmagorinsky;
            _isf.kinematicViscosity = options.KinematicViscosity;
            _isf.UpdateSpace(options.UseLES, options.GravityPsi2);
            _isf.UpdateVelocities(_velocity);
            _simTime += _isf.dt;
        }

        public void Dispose()
        {
            _velocity?.Dispose();
            _isf?.Dispose();
        }
    }
}
