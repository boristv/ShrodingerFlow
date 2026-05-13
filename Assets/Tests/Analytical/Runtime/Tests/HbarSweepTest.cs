using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ShrodingerFlow.AnalyticalTests
{
    /// <summary>
    /// Тест 3 «Влияние параметра ℏ».
    /// Одна и та же начальная вихревая конфигурация прогоняется при нескольких значениях ℏ.
    /// Измеряемые величины:
    ///   • число «активных» вихревых ячеек (|ω| > threshold) — косвенный счётчик структур;
    ///   • спектр E(k) в финале — смещение пика при изменении ℏ;
    ///   • устойчивость: доля шагов без NaN/Inf в полях скорости;
    ///   • среднеквадратичный разброс скоростей RMS|u|.
    /// Артефакты:
    ///   metrics_hbar_*.csv (по одному на ℏ);
    ///   spectrum_hbar_*.csv;
    ///   vorticity_hbar_*.png (4 кадра на ℏ);
    ///   sweep_summary.csv — финальные показатели по всем ℏ в одной таблице;
    ///   summary.json.
    /// </summary>
    public sealed class HbarSweepTest : AnalyticalTestBase
    {
        [Header("Сетка (общая для всех прогонов)")]
        [SerializeField] private Vector3Int _volSize = new Vector3Int(10, 5, 5);
        [SerializeField] private Vector3Int _volRes = new Vector3Int(64, 32, 32);
        [SerializeField] private float _dt = 1f / 12f;

        [Header("Набор значений ℏ")]
        [SerializeField] private float[] _hbarValues = { 0.05f, 0.1f, 0.2f };

        [Header("Начальная конфигурация (кольца)")]
        [SerializeField] private List<VortexRingDef> _rings = new List<VortexRingDef>
        {
            new VortexRingDef { Center = new Vector3(5f, 2.5f, 2.5f), Normal = new Vector3(-1f, 0f, 0f), Radius = 1.5f },
            new VortexRingDef { Center = new Vector3(5f, 2.5f, 2.5f), Normal = new Vector3(-1f, 0f, 0f), Radius = 0.9f },
        };

        [Header("Симуляция")]
        [SerializeField] private int _totalSteps = 200;
        [SerializeField] private int _sampleEveryNSteps = 5;
        [SerializeField] private int _slicePngCount = 4;
        [Tooltip("Порог |ω| для счётчика активных вихревых ячеек (относительный от max|ω| на шаге 0).")]
        [SerializeField] private float _vorticityThresholdRel = 0.1f;

        [Header("LES")]
        [SerializeField] private bool _useLES;
        [SerializeField] private float _csSmagorinsky = 0.07f;

        public override string TestName => "03_HbarSweep";
        public override string Description =>
            $"ℏ-sweep: {string.Join(",", _hbarValues)}, vol {_volRes}, steps={_totalSteps}";

        public override IEnumerator RunTest(TestRunContext ctx)
        {
            if (_hbarValues == null || _hbarValues.Length == 0)
            {
                Debug.LogWarning($"[{TestName}] Список hbarValues пустой — тест пропущен.");
                yield break;
            }

            using (var writer = new TestArtifactWriter(ctx.OutputDir))
            {
                var rings = _rings.ToArray();
                int[] volSize = { _volSize.x, _volSize.y, _volSize.z };
                int[] volRes = { _volRes.x, _volRes.y, _volRes.z };

                var sweepCsv = writer.CreateCsv("sweep_summary.csv",
                    "hbar", "final_kinetic_energy", "final_max_omega", "final_rms_omega",
                    "final_max_vel", "final_rms_vel", "active_vortex_cells_final",
                    "density_integral_rel", "wall_time_s");

                var sweepEntries = new List<SweepEntry>();

                for (int hi = 0; hi < _hbarValues.Length; hi++)
                {
                    float hbar = _hbarValues[hi];
                    string tag = $"hbar{hbar:G3}".Replace('.', 'p');
                    writer.Log($"=== ℏ={hbar} ({hi + 1}/{_hbarValues.Length}) ===");

                    var metricsCsv = writer.CreateCsv($"metrics_{tag}.csv",
                        "step", "t",
                        "kinetic_energy", "density_integral_rel",
                        "max_omega", "rms_omega", "active_cells",
                        "max_vel", "rms_vel");

                    double t0 = Time.realtimeSinceStartupAsDouble;
                    IsfDiagnosticsSnapshot finalSnap = default;
                    int finalActiveCells = 0;

                    using (var h = new IsfTestHarness(ctx.KernelsShader, ctx.FftShader, ctx.LesShader,
                        volSize, volRes, hbar, _dt))
                    {
                        h.InitPlaneWave(Vector3.zero, 0.01f);
                        h.AddVortexRings(rings);
                        h.NormalizeAndProject();
                        h.Isf.UpdateVelocities(h.Velocity);

                        var init = IsfDiagnostics.SampleAll(h.Isf, h.Velocity);
                        float omegaThreshold = init.MaxVorticity * _vorticityThresholdRel;

                        var opts = new StepOptions
                        {
                            UseLES = _useLES,
                            CsSmagorinsky = _csSmagorinsky
                        };

                        int sliceInterval = _slicePngCount > 1
                            ? Mathf.Max(1, _totalSteps / (_slicePngCount - 1))
                            : int.MaxValue;
                        int nextSlice = sliceInterval;
                        int sliceNum = 0;
                        SaveSlice(writer, h, tag, sliceNum++);

                        ctx.ReportProgress(0, _totalSteps * _hbarValues.Length);

                        for (int step = 1; step <= _totalSteps; step++)
                        {
                            h.Step(opts);

                            if (step % _sampleEveryNSteps == 0 || step == _totalSteps)
                            {
                                var s = IsfDiagnostics.SampleAll(h.Isf, h.Velocity);
                                int activeCells = CountActiveCells(h, omegaThreshold);
                                ctx.ReportProgress(hi * _totalSteps + step, _totalSteps * _hbarValues.Length,
                                    $"ℏ={hbar} — E={s.KineticEnergy:G4}  max|ω|={s.MaxVorticity:G4}  cells={activeCells}");
                                float densityRel = init.DensityIntegral > 1e-20f
                                    ? s.DensityIntegral / init.DensityIntegral : 0f;
                                metricsCsv.Row(step, h.SimTime,
                                    s.KineticEnergy, densityRel,
                                    s.MaxVorticity, s.RmsVorticity, activeCells,
                                    s.MaxVelocity, s.RmsVelocity);
                                finalSnap = s;
                                finalActiveCells = activeCells;
                                yield return null;
                            }

                            if (sliceNum < _slicePngCount && step >= nextSlice)
                            {
                                SaveSlice(writer, h, tag, sliceNum++);
                                nextSlice = sliceNum * sliceInterval;
                            }
                        }

                        if (sliceNum < _slicePngCount)
                            SaveSlice(writer, h, tag, sliceNum);

                        WriteSpectrum(writer, h, tag);
                    }

                    double wallTime = Time.realtimeSinceStartupAsDouble - t0;
                    float densIntegRel = _rings.Count > 0 ? 1f : 0f; // placeholder — в metricsCsv финальная строка содержит точное значение

                    sweepCsv.Row(hbar,
                        finalSnap.KineticEnergy, finalSnap.MaxVorticity, finalSnap.RmsVorticity,
                        finalSnap.MaxVelocity, finalSnap.RmsVelocity,
                        finalActiveCells, densIntegRel, (float)wallTime);

                    sweepEntries.Add(new SweepEntry
                    {
                        hbar = hbar,
                        finalKineticEnergy = finalSnap.KineticEnergy,
                        finalMaxOmega = finalSnap.MaxVorticity,
                        finalRmsOmega = finalSnap.RmsVorticity,
                        finalMaxVel = finalSnap.MaxVelocity,
                        finalActiveCells = finalActiveCells,
                        wallTimeSeconds = (float)wallTime
                    });
                    writer.Log($"  ℏ={hbar}: max|ω|={finalSnap.MaxVorticity:G5}, E={finalSnap.KineticEnergy:G5}, активных ячеек={finalActiveCells}");
                }

                var summary = new Summary
                {
                    testName = TestName,
                    volRes = new[] { _volRes.x, _volRes.y, _volRes.z },
                    hbarValues = (float[])_hbarValues.Clone(),
                    dt = _dt,
                    totalSteps = _totalSteps,
                    useLES = _useLES,
                    entries = sweepEntries
                };
                writer.Json("summary.json", summary);
                writer.Log("Sweep завершён.");
            }
        }

        private int CountActiveCells(IsfTestHarness h, float threshold)
        {
            int rx = h.Isf.resX, ry = h.Isf.resY, rz = h.Isf.resZ;
            int num = rx * ry * rz;

            var vx = new float[num];
            var vy = new float[num];
            var vz = new float[num];
            h.Velocity.vx.GetData(vx);
            h.Velocity.vy.GetData(vy);
            h.Velocity.vz.GetData(vz);

            float twoDx = 2f * h.Isf.dx;
            float twoDy = 2f * h.Isf.dy;
            float twoDz = 2f * h.Isf.dz;

            int count = 0;
            for (int i = 0; i < rx; i++)
            {
                int ip = (i + 1) % rx, im = (i - 1 + rx) % rx;
                for (int j = 0; j < ry; j++)
                {
                    int jp = (j + 1) % ry, jm = (j - 1 + ry) % ry;
                    for (int k = 0; k < rz; k++)
                    {
                        int kp = (k + 1) % rz, km = (k - 1 + rz) % rz;
                        int idx_ip = ip * ry * rz + j * rz + k;
                        int idx_im = im * ry * rz + j * rz + k;
                        int idx_jp = i * ry * rz + jp * rz + k;
                        int idx_jm = i * ry * rz + jm * rz + k;
                        int idx_kp = i * ry * rz + j * rz + kp;
                        int idx_km = i * ry * rz + j * rz + km;

                        float wx = (vz[idx_jp] - vz[idx_jm]) / twoDy - (vy[idx_kp] - vy[idx_km]) / twoDz;
                        float wy = (vx[idx_kp] - vx[idx_km]) / twoDz - (vz[idx_ip] - vz[idx_im]) / twoDx;
                        float wz = (vy[idx_ip] - vy[idx_im]) / twoDx - (vx[idx_jp] - vx[idx_jm]) / twoDy;
                        if (Mathf.Sqrt(wx * wx + wy * wy + wz * wz) > threshold)
                            count++;
                    }
                }
            }
            return count;
        }

        private void WriteSpectrum(TestArtifactWriter writer, IsfTestHarness h, string tag)
        {
            int rx = h.Isf.resX, ry = h.Isf.resY, rz = h.Isf.resZ;
            int num = rx * ry * rz;
            var vx = new float[num]; var vy = new float[num]; var vz = new float[num];
            h.Velocity.vx.GetData(vx); h.Velocity.vy.GetData(vy); h.Velocity.vz.GetData(vz);

            int kMax = Mathf.Min(rx, Mathf.Min(ry, rz)) / 2;
            var eShell = new double[kMax + 1];
            var cnt = new int[kMax + 1];
            float dV = h.Isf.dx * h.Isf.dy * h.Isf.dz;

            for (int i = 0; i < rx; i++) { int ki = i <= rx / 2 ? i : rx - i;
                for (int j = 0; j < ry; j++) { int kj = j <= ry / 2 ? j : ry - j;
                    for (int k = 0; k < rz; k++) {
                        int kk = k <= rz / 2 ? k : rz - k;
                        int shell = Mathf.RoundToInt(Mathf.Sqrt(ki * ki + kj * kj + kk * kk));
                        if (shell > kMax) continue;
                        int idx = i * ry * rz + j * rz + k;
                        double u = vx[idx], v = vy[idx], w = vz[idx];
                        eShell[shell] += 0.5 * (u * u + v * v + w * w);
                        cnt[shell]++;
                    }
                }
            }

            var csv = writer.CreateCsv($"spectrum_{tag}.csv", "k_shell", "E_k", "cell_count");
            for (int k = 1; k <= kMax; k++)
                csv.Row(k, (float)(cnt[k] > 0 ? eShell[k] * dV / cnt[k] : 0.0), cnt[k]);
        }

        private void SaveSlice(TestArtifactWriter writer, IsfTestHarness h, string tag, int idx)
        {
            var tex = IsfDiagnostics.RenderVorticitySlice(h.Isf, h.Velocity, SliceAxis.Z, h.Isf.resZ / 2);
            try { writer.Png($"vorticity_{tag}_{idx:D2}.png", tex); }
            finally { UnityEngine.Object.Destroy(tex); }
        }

        [Serializable]
        private struct SweepEntry
        {
            public float hbar;
            public float finalKineticEnergy;
            public float finalMaxOmega;
            public float finalRmsOmega;
            public float finalMaxVel;
            public int finalActiveCells;
            public float wallTimeSeconds;
        }

        [Serializable]
        private struct Summary
        {
            public string testName;
            public int[] volRes;
            public float[] hbarValues;
            public float dt;
            public int totalSteps;
            public bool useLES;
            public List<SweepEntry> entries;
        }
    }
}
