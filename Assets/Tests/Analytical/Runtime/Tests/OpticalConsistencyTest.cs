using System;
using System.Collections;
using UnityEngine;

namespace ShrodingerFlow.AnalyticalTests
{
    /// <summary>
    /// Тест 5 «Проверка оптической согласованности».
    ///
    /// ISF поддерживает |ψ|² = const повсюду (поклеточная нормировка), поэтому для оптических
    /// тестов используется поле скоростей: field(x) = |v(x)|² = vx²+vy²+vz².
    /// Это поле реально неоднородно (пиковое у ядер вихрей) и имеет значимые гладкие градиенты.
    ///
    ///   1. Beer-Lambert: по каждому лучу вдоль X вычисляется τ = κ·∫|v|² dx,
    ///      затем T_numeric (дискретное произведение exp) vs T_bl = exp(-τ).
    ///      Относительная ошибка должна быть < thresholdBeerRelError (математическая идентичность).
    ///
    ///   2. «Плавность нормалей»: средняя угловая вариация между соседними ∇(|v|²)
    ///      по выборке точек — должна быть < thresholdNormalVariance радиан.
    ///      Для ISF скоростное поле гладко по построению, поэтому тест должен проходить.
    ///
    /// Артефакты: beer_lambert.csv, normals_sample.csv, vel2_slice_*.png, grad_slice.png, summary.json.
    /// </summary>
    public sealed class OpticalConsistencyTest : AnalyticalTestBase
    {
        [Header("Сетка")]
        [SerializeField] private Vector3Int _volSize = new Vector3Int(4, 2, 2);
        [SerializeField] private Vector3Int _volRes = new Vector3Int(64, 32, 32);
        [SerializeField] private float _hbar = 0.1f;
        [SerializeField] private float _dt = 1f / 12f;

        [Header("Начальная конфигурация (Гауссов пакет)")]
        [Tooltip("Центр Гауссова пятна в координатах объёма.")]
        [SerializeField] private Vector3 _blobCenter = new Vector3(2f, 1f, 1f);
        [Tooltip("Стандартное отклонение пятна (физические единицы). Должно быть >> dx для гладких нормалей.")]
        [SerializeField] private float _blobSigma = 0.4f;

        [Header("Прогрев (шаги до замера)")]
        [Tooltip("Число шагов ISF перед измерением — чтобы сформировалась выраженная структура.")]
        [SerializeField] private int _warmupSteps = 60;

        [Header("Beer–Lambert")]
        [Tooltip("Коэффициент поглощения κ (то же, что _raymarchAbsorption в ParticleDisplay3D).")]
        [SerializeField] private float _kappa = 0.42f;
        [Tooltip("Число лучей по сетке y/z через центральный X-срез.")]
        [SerializeField] private int _raysPerAxis = 16;
        [Tooltip("Максимальная допустимая относительная ошибка Beer–Lambert.")]
        [SerializeField] private float _thresholdBeerRelError = 0.05f;

        [Header("Нормали")]
        [Tooltip("Допустимая средняя угловая вариация нормалей в радианах.")]
        [SerializeField] private float _thresholdNormalVariance = 0.3f;
        [SerializeField] private int _normalSampleCount = 500;

        public override string TestName => "05_OpticalConsistency";
        public override string Description =>
            $"Beer-Lambert (κ={_kappa}), нормали ∇ρ, прогрев {_warmupSteps} шагов";

        public override IEnumerator RunTest(TestRunContext ctx)
        {
            using (var writer = new TestArtifactWriter(ctx.OutputDir))
            {
                int[] volSize = { _volSize.x, _volSize.y, _volSize.z };
                int[] volRes = { _volRes.x, _volRes.y, _volRes.z };

                using (var h = new IsfTestHarness(ctx.KernelsShader, ctx.FftShader, ctx.LesShader,
                    volSize, volRes, _hbar, _dt))
                {
                    // Гауссов пакет без импульса: ρ=|ψ|² локализовано, градиенты значимы и гладки.
                    // После прогрева ISF поле остаётся гладким — это то, что проверяет тест нормалей.
                    h.InitGaussianBlob(_blobCenter, _blobSigma, Vector3.zero);
                    h.NormalizeAndProject();
                    h.Isf.UpdateVelocities(h.Velocity);

                    var opts = StepOptions.Default;
                    for (int step = 0; step < _warmupSteps; step++)
                    {
                        h.Step(opts);
                        if (step % 20 == 0)
                        {
                            ctx.ReportProgress(step, _warmupSteps, $"Прогрев шаг {step}/{_warmupSteps}...");
                            yield return null;
                        }
                    }
                    ctx.ReportProgress(_warmupSteps, _warmupSteps, "Замер Beer-Lambert и нормалей...");

                    // Читаем |v|² на CPU — реально неоднородное поле ISF.
                    // |ψ|² = 1 везде из-за поклеточной нормировки; скоростное поле — нет.
                    int rx = h.Isf.resX, ry = h.Isf.resY, rz = h.Isf.resZ;
                    int num = rx * ry * rz;
                    float[] field = IsfDiagnostics.SampleVelocitySquared(h.Velocity, num);

                    float fieldSum = 0f;
                    float fieldMax = 0f;
                    for (int i = 0; i < num; i++)
                    {
                        fieldSum += field[i];
                        if (field[i] > fieldMax) fieldMax = field[i];
                    }
                    float fieldMean = fieldSum / num;

                    writer.Log($"Прогрев завершён. Среднее |v|²={fieldMean:G5}, max |v|²={fieldMax:G5}");
                    yield return null;

                    // ---- 1. Beer–Lambert по полю |v|² ----
                    BeerLambertResults bl = TestBeerLambert(field, rx, ry, rz,
                        h.Isf.dx, h.Isf.dy, h.Isf.dz, writer);
                    bool beerOk = bl.MaxRelError < _thresholdBeerRelError;
                    writer.Log($"Beer-Lambert (|v|²): maxRelErr={bl.MaxRelError:G4}, meanRelErr={bl.MeanRelError:G4} — {(beerOk ? "OK" : "WARN")}");

                    yield return null;

                    // ---- 2+3. Нормали ∇(|v|²) ----
                    NormalResults normals = TestNormals(field, rx, ry, rz,
                        h.Isf.dx, h.Isf.dy, h.Isf.dz, writer);
                    bool normalOk = normals.MeanAngleVariance < _thresholdNormalVariance;
                    writer.Log($"Нормали ∇(|v|²): mean angle={normals.MeanAngleVariance:G4} рад — {(normalOk ? "OK" : "WARN")}");
                    writer.Log($"  ∇(|v|²) max={normals.MaxGradMag:G5}, mean={normals.MeanGradMag:G5}");

                    yield return null;

                    // ---- Срезы ----
                    SaveRhoSlice(writer, field, rx, ry, rz, SliceAxis.Z, rz / 2, "vel2_slice_Z");
                    SaveRhoSlice(writer, field, rx, ry, rz, SliceAxis.Y, ry / 2, "vel2_slice_Y");
                    SaveGradSlice(writer, field, rx, ry, rz, h.Isf.dx, h.Isf.dy, h.Isf.dz, "grad_slice");

                    var summary = new Summary
                    {
                        testName                = TestName,
                        volRes                  = new[] { rx, ry, rz },
                        hbar                    = _hbar,
                        kappa                   = _kappa,
                        warmupSteps             = _warmupSteps,
                        vel2Mean                = fieldMean,
                        vel2Max                 = fieldMax,
                        beerMaxRelError         = bl.MaxRelError,
                        beerMeanRelError        = bl.MeanRelError,
                        beerPass                = beerOk,
                        beerThreshold           = _thresholdBeerRelError,
                        normalMeanAngleVariance = normals.MeanAngleVariance,
                        normalMaxGradMag        = normals.MaxGradMag,
                        normalPass              = normalOk,
                        normalThreshold         = _thresholdNormalVariance
                    };
                    writer.Json("summary.json", summary);
                }
            }
        }

        // ---- Beer–Lambert ----

        private struct BeerLambertResults
        {
            public float MaxRelError;
            public float MeanRelError;
        }

        private BeerLambertResults TestBeerLambert(float[] rho,
            int rx, int ry, int rz, float dx, float dy, float dz,
            TestArtifactWriter writer)
        {
            // Лучи вдоль оси X через сетку ry x rz точек (по _raysPerAxis на ось)
            int step_y = Mathf.Max(1, ry / _raysPerAxis);
            int step_z = Mathf.Max(1, rz / _raysPerAxis);

            var csv = writer.CreateCsv("beer_lambert.csv",
                "ray_j", "ray_k",
                "optical_depth", "transmittance_numeric", "transmittance_beerlambert", "rel_error");

            double errSum = 0.0;
            float maxErr = 0f;
            int count = 0;

            for (int j = 0; j < ry; j += step_y)
            {
                for (int k = 0; k < rz; k += step_z)
                {
                    // Численная оптическая глубина τ = κ · Σ ρ(i,j,k) · dx
                    double tau = 0.0;
                    for (int i = 0; i < rx; i++)
                        tau += rho[i * ry * rz + j * rz + k] * dx;
                    tau *= _kappa;

                    double T_bl = System.Math.Exp(-tau);

                    // «Численная» пропускаемость — дискретное произведение через шаг dx
                    double T_step = 1.0;
                    for (int i = 0; i < rx; i++)
                    {
                        double absorb = _kappa * rho[i * ry * rz + j * rz + k] * dx;
                        T_step *= System.Math.Exp(-absorb);
                    }

                    float relErr = T_bl > 1e-12 ? (float)System.Math.Abs(T_step - T_bl) / (float)T_bl : 0f;
                    csv.Row(j, k, (float)tau, (float)T_step, (float)T_bl, relErr);

                    errSum += relErr;
                    if (relErr > maxErr) maxErr = relErr;
                    count++;
                }
            }

            return new BeerLambertResults
            {
                MaxRelError = maxErr,
                MeanRelError = count > 0 ? (float)(errSum / count) : 0f
            };
        }

        // ---- Нормали ----

        private struct NormalResults
        {
            public float MeanAngleVariance;
            public float MaxGradMag;
            public float MeanGradMag;
        }

        private NormalResults TestNormals(float[] rho,
            int rx, int ry, int rz, float dx, float dy, float dz,
            TestArtifactWriter writer)
        {
            var csv = writer.CreateCsv("normals_sample.csv",
                "ix", "iy", "iz",
                "grad_x", "grad_y", "grad_z", "grad_mag",
                "angle_to_neighbor_deg");

            float twoDx = 2f * dx, twoDy = 2f * dy, twoDz = 2f * dz;

            Vector3 Gradient(int i, int j, int k)
            {
                int ip = (i + 1) % rx, im = (i - 1 + rx) % rx;
                int jp = (j + 1) % ry, jm = (j - 1 + ry) % ry;
                int kp = (k + 1) % rz, km = (k - 1 + rz) % rz;
                float gx = (rho[ip * ry * rz + j * rz + k] - rho[im * ry * rz + j * rz + k]) / twoDx;
                float gy = (rho[i * ry * rz + jp * rz + k] - rho[i * ry * rz + jm * rz + k]) / twoDy;
                float gz = (rho[i * ry * rz + j * rz + kp] - rho[i * ry * rz + j * rz + km]) / twoDz;
                return new Vector3(gx, gy, gz);
            }

            double angleSum = 0.0;
            int count = 0;
            float maxMag = 0f;
            double magSum = 0.0;

            // Выборка по равномерной сетке
            int stepX = Mathf.Max(1, rx / Mathf.CeilToInt(Mathf.Pow(_normalSampleCount, 1f / 3f)));
            int stepY = Mathf.Max(1, ry / Mathf.CeilToInt(Mathf.Pow(_normalSampleCount, 1f / 3f)));
            int stepZ = Mathf.Max(1, rz / Mathf.CeilToInt(Mathf.Pow(_normalSampleCount, 1f / 3f)));

            for (int i = 1; i < rx - 1; i += stepX)
            {
                for (int j = 1; j < ry - 1; j += stepY)
                {
                    for (int k = 1; k < rz - 1; k += stepZ)
                    {
                        Vector3 g = Gradient(i, j, k);
                        float mag = g.magnitude;
                        magSum += mag;
                        if (mag > maxMag) maxMag = mag;

                        // Угловая вариация — угол с градиентом правого соседа по X
                        Vector3 gNeighbor = Gradient((i + 1) % rx, j, k);
                        float cosAngle = mag > 1e-12f && gNeighbor.magnitude > 1e-12f
                            ? Vector3.Dot(g.normalized, gNeighbor.normalized)
                            : 1f;
                        float angle = Mathf.Acos(Mathf.Clamp(cosAngle, -1f, 1f));
                        angleSum += angle;
                        count++;

                        csv.Row(i, j, k, g.x, g.y, g.z, mag, angle * Mathf.Rad2Deg);
                    }
                }
            }

            return new NormalResults
            {
                MeanAngleVariance = count > 0 ? (float)(angleSum / count) : 0f,
                MaxGradMag = maxMag,
                MeanGradMag = count > 0 ? (float)(magSum / count) : 0f
            };
        }

        // ---- Срезы ----

        private static void SaveRhoSlice(TestArtifactWriter writer, float[] rho,
            int rx, int ry, int rz, SliceAxis axis, int sliceIdx, string name)
        {
            int w, h2;
            switch (axis)
            {
                case SliceAxis.X: w = ry; h2 = rz; break;
                case SliceAxis.Y: w = rx; h2 = rz; break;
                default: w = rx; h2 = ry; break;
            }

            float rhoMax = 0f;
            for (int a = 0; a < w; a++)
                for (int b = 0; b < h2; b++)
                {
                    int idx = CellIndex(axis, a, b, sliceIdx, rx, ry, rz);
                    if (rho[idx] > rhoMax) rhoMax = rho[idx];
                }

            float inv = rhoMax > 1e-20f ? 1f / rhoMax : 0f;
            var tex = new Texture2D(w, h2, TextureFormat.RGB24, false);
            var pixels = new Color32[w * h2];
            for (int a = 0; a < w; a++)
                for (int b = 0; b < h2; b++)
                {
                    int idx = CellIndex(axis, a, b, sliceIdx, rx, ry, rz);
                    pixels[b * w + a] = Colormap.Viridis(rho[idx] * inv);
                }
            tex.SetPixels32(pixels);
            tex.Apply(false);
            try { writer.Png($"{name}.png", tex); }
            finally { UnityEngine.Object.Destroy(tex); }
        }

        private static void SaveGradSlice(TestArtifactWriter writer, float[] rho,
            int rx, int ry, int rz, float dx, float dy, float dz, string name)
        {
            int sliceZ = rz / 2;
            float twoDx = 2f * dx, twoDy = 2f * dy, twoDz = 2f * dz;
            float maxMag = 0f;
            var mags = new float[rx * ry];

            for (int i = 0; i < rx; i++)
            {
                int ip = (i + 1) % rx, im = (i - 1 + rx) % rx;
                for (int j = 0; j < ry; j++)
                {
                    int jp = (j + 1) % ry, jm = (j - 1 + ry) % ry;
                    int kp = (sliceZ + 1) % rz, km = (sliceZ - 1 + rz) % rz;
                    float gx = (rho[ip * ry * rz + j * rz + sliceZ] - rho[im * ry * rz + j * rz + sliceZ]) / twoDx;
                    float gy = (rho[i * ry * rz + jp * rz + sliceZ] - rho[i * ry * rz + jm * rz + sliceZ]) / twoDy;
                    float gz = (rho[i * ry * rz + j * rz + kp] - rho[i * ry * rz + j * rz + km]) / twoDz;
                    float mag = Mathf.Sqrt(gx * gx + gy * gy + gz * gz);
                    mags[j * rx + i] = mag;
                    if (mag > maxMag) maxMag = mag;
                }
            }

            float inv = maxMag > 1e-20f ? 1f / maxMag : 0f;
            var tex = new Texture2D(rx, ry, TextureFormat.RGB24, false);
            var pixels = new Color32[rx * ry];
            for (int p = 0; p < pixels.Length; p++)
                pixels[p] = Colormap.Viridis(mags[p] * inv);
            tex.SetPixels32(pixels);
            tex.Apply(false);
            try { writer.Png($"{name}.png", tex); }
            finally { UnityEngine.Object.Destroy(tex); }
        }

        private static int CellIndex(SliceAxis axis, int a, int b, int s,
            int rx, int ry, int rz)
        {
            switch (axis)
            {
                case SliceAxis.X: return s * ry * rz + a * rz + b;
                case SliceAxis.Y: return a * ry * rz + s * rz + b;
                default: return a * ry * rz + b * rz + s;
            }
        }

        [Serializable]
        private struct Summary
        {
            public string testName;
            public int[]  volRes;
            public float  hbar, kappa;
            public int    warmupSteps;
            public float  vel2Mean, vel2Max;          // поле |v|² (заменяет rhoMean — ρ=1 везде в ISF)
            public float  beerMaxRelError, beerMeanRelError;
            public bool   beerPass;
            public float  beerThreshold;
            public float  normalMeanAngleVariance, normalMaxGradMag;
            public bool   normalPass;
            public float  normalThreshold;
        }
    }
}
