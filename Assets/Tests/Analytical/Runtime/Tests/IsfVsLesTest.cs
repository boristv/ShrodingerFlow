using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ShrodingerFlow.AnalyticalTests
{
    /// <summary>
    /// Тест 2 «Сравнение ISF и ISF+LES».
    /// Один и тот же вихревой сценарий прогоняется дважды (без LES и с LES).
    /// Измеряемые величины:
    ///   • шум скорости: RMS|u| и max|u| со временем;
    ///   • max|ω|, RMS|ω| — устойчивость мелких структур;
    ///   • спектр кинетической энергии E(k) — пишется в отдельный CSV;
    ///   • визуальный срез |ω| для каждого сценария.
    /// Артефакты:
    ///   metrics_noLES.csv / metrics_LES.csv — поштучные ряды;
    ///   spectrum_noLES.csv / spectrum_LES.csv — радиальный спектр E(k) (финальный кадр);
    ///   summary.json;
    ///   vorticity_slice_*.png для обоих прогонов.
    /// </summary>
    public sealed class IsfVsLesTest : AnalyticalTestBase
    {
        [Header("Сетка")]
        [SerializeField] private Vector3Int _volSize = new Vector3Int(10, 5, 5);
        [SerializeField] private Vector3Int _volRes = new Vector3Int(128, 64, 64);
        [SerializeField] private float _hbar = 0.1f;
        [SerializeField] private float _dt = 1f / 12f;

        [Header("Начальная конфигурация")]
        [SerializeField] private List<VortexRingDef> _rings = new List<VortexRingDef>
        {
            new VortexRingDef { Center = new Vector3(5f, 2.5f, 2.5f), Normal = new Vector3(-1f, 0f, 0f), Radius = 1.5f },
            new VortexRingDef { Center = new Vector3(5f, 2.5f, 2.5f), Normal = new Vector3(-1f, 0f, 0f), Radius = 0.9f },
        };

        [Header("Симуляция")]
        [SerializeField] private int _totalSteps = 300;
        [SerializeField] private int _sampleEveryNSteps = 5;
        [SerializeField] private int _slicePngCount = 3;

        [Header("LES параметры")]
        [SerializeField] private float _csSmagorinsky = 0.07f;
        [SerializeField] private float _kinematicViscosity;

        public override string TestName => "02_IsfVsLes";
        public override string Description =>
            $"ISF vs ISF+LES, vol {_volRes}, hbar={_hbar}, dt={_dt:G3}, steps={_totalSteps}";

        public override IEnumerator RunTest(TestRunContext ctx)
        {
            using (var writer = new TestArtifactWriter(ctx.OutputDir))
            {
                var summary = new Summary
                {
                    testName = TestName,
                    volRes = new[] { _volRes.x, _volRes.y, _volRes.z },
                    hbar = _hbar,
                    dt = _dt,
                    totalSteps = _totalSteps,
                    csSmagorinsky = _csSmagorinsky,
                    kinematicViscosity = _kinematicViscosity,
                };

                int[] volSize = { _volSize.x, _volSize.y, _volSize.z };
                int[] volRes = { _volRes.x, _volRes.y, _volRes.z };
                var rings = _rings.ToArray();

                // ------ прогон без LES ------
                writer.Log("=== Прогон 1/2: без LES ===");
                var metricsNoLes = writer.CreateCsv("metrics_noLES.csv",
                    "step", "t", "kinetic_energy", "max_omega", "rms_omega", "max_vel", "rms_vel");
                var snapshotsNoLes = new List<IsfDiagnosticsSnapshot>();

                IsfDiagnosticsSnapshot finalNoLes;
                yield return RunOnce(ctx, volSize, volRes, rings, false, metricsNoLes, snapshotsNoLes, writer, "noLES");
                finalNoLes = snapshotsNoLes.Count > 0 ? snapshotsNoLes[snapshotsNoLes.Count - 1] : default;

                // ------ прогон с LES ------
                writer.Log("=== Прогон 2/2: с LES ===");
                var metricsLes = writer.CreateCsv("metrics_LES.csv",
                    "step", "t", "kinetic_energy", "max_omega", "rms_omega", "max_vel", "rms_vel");
                var snapshotsLes = new List<IsfDiagnosticsSnapshot>();

                IsfDiagnosticsSnapshot finalLes;
                yield return RunOnce(ctx, volSize, volRes, rings, true, metricsLes, snapshotsLes, writer, "LES");
                finalLes = snapshotsLes.Count > 0 ? snapshotsLes[snapshotsLes.Count - 1] : default;

                // ------ сравнительная таблица ------
                var comparison = writer.CreateCsv("comparison.csv",
                    "metric", "noLES", "LES", "ratio_LES_noLES");
                comparison.Row("kinetic_energy_final", finalNoLes.KineticEnergy, finalLes.KineticEnergy,
                    SafeDiv(finalLes.KineticEnergy, finalNoLes.KineticEnergy));
                comparison.Row("max_omega_final", finalNoLes.MaxVorticity, finalLes.MaxVorticity,
                    SafeDiv(finalLes.MaxVorticity, finalNoLes.MaxVorticity));
                comparison.Row("rms_omega_final", finalNoLes.RmsVorticity, finalLes.RmsVorticity,
                    SafeDiv(finalLes.RmsVorticity, finalNoLes.RmsVorticity));
                comparison.Row("max_vel_final", finalNoLes.MaxVelocity, finalLes.MaxVelocity,
                    SafeDiv(finalLes.MaxVelocity, finalNoLes.MaxVelocity));
                comparison.Row("rms_vel_final", finalNoLes.RmsVelocity, finalLes.RmsVelocity,
                    SafeDiv(finalLes.RmsVelocity, finalNoLes.RmsVelocity));

                summary.noLesFinal = ToSer(finalNoLes);
                summary.lesFinal = ToSer(finalLes);
                writer.Json("summary.json", summary);
                writer.Log("Готово.");
            }
        }

        private IEnumerator RunOnce(TestRunContext ctx, int[] volSize, int[] volRes,
            VortexRingDef[] rings, bool useLES, CsvLog metrics,
            List<IsfDiagnosticsSnapshot> snapshots, TestArtifactWriter writer, string tag)
        {
            using (var h = new IsfTestHarness(ctx.KernelsShader, ctx.FftShader, ctx.LesShader,
                volSize, volRes, _hbar, _dt))
            {
                h.InitPlaneWave(Vector3.zero, 0.01f);
                h.AddVortexRings(rings);
                h.NormalizeAndProject();
                h.Isf.UpdateVelocities(h.Velocity);

                var opts = new StepOptions
                {
                    UseLES = useLES,
                    KinematicViscosity = _kinematicViscosity,
                    CsSmagorinsky = _csSmagorinsky
                };

                // Срезы |ω|: step 0 + равномерно по прогрессу
                int nextSliceStep = _slicePngCount > 1
                    ? Mathf.Max(1, _totalSteps / (_slicePngCount - 1))
                    : int.MaxValue;
                int sliceIdx = 0;
                SaveSlice(writer, h, tag, sliceIdx++);

                ctx.ReportProgress(0, _totalSteps * 2, $"[{tag}] прогрев...");

                for (int step = 1; step <= _totalSteps; step++)
                {
                    h.Step(opts);

                    if (step % _sampleEveryNSteps == 0 || step == _totalSteps)
                    {
                        var s = IsfDiagnostics.SampleAll(h.Isf, h.Velocity);
                        snapshots.Add(s);
                        metrics.Row(step, h.SimTime,
                            s.KineticEnergy, s.MaxVorticity, s.RmsVorticity,
                            s.MaxVelocity, s.RmsVelocity);

                        int globalStep = (useLES ? _totalSteps : 0) + step;
                        ctx.ReportProgress(globalStep, _totalSteps * 2,
                            $"[{tag}] E={s.KineticEnergy:G4}  max|ω|={s.MaxVorticity:G4}  RMS|u|={s.RmsVelocity:G4}");
                        yield return null;
                    }

                    if (sliceIdx < _slicePngCount && step >= nextSliceStep * sliceIdx)
                    {
                        SaveSlice(writer, h, tag, sliceIdx++);
                    }
                }

                if (sliceIdx < _slicePngCount)
                    SaveSlice(writer, h, tag, sliceIdx);

                // Записать спектр E(k) финального кадра
                WriteEnergySpectrum(writer, h, tag);
                writer.Log($"{tag}: final rms|ω|={snapshots[snapshots.Count - 1].RmsVorticity:G5}, max|u|={snapshots[snapshots.Count - 1].MaxVelocity:G5}");
            }
        }

        /// <summary>
        /// Радиальный (shell-averaged) спектр кинетической энергии по компонентам скорости.
        /// E(k) = Σ_{|k_vec|≈k} 0.5·|û(k_vec)|² · (volSize/num)
        /// Спектр пишется как "shell_bin_k, E_k".
        /// </summary>
        private void WriteEnergySpectrum(TestArtifactWriter writer, IsfTestHarness h, string tag)
        {
            int rx = h.Isf.resX, ry = h.Isf.resY, rz = h.Isf.resZ;
            int num = rx * ry * rz;

            var vx = new float[num];
            var vy = new float[num];
            var vz = new float[num];
            h.Velocity.vx.GetData(vx);
            h.Velocity.vy.GetData(vy);
            h.Velocity.vz.GetData(vz);

            int kMax = Mathf.Min(rx, Mathf.Min(ry, rz)) / 2;
            var eShell = new double[kMax + 1];
            var countShell = new int[kMax + 1];

            float dV = h.Isf.dx * h.Isf.dy * h.Isf.dz;

            for (int i = 0; i < rx; i++)
            {
                int ki = i <= rx / 2 ? i : rx - i;
                for (int j = 0; j < ry; j++)
                {
                    int kj = j <= ry / 2 ? j : ry - j;
                    for (int k = 0; k < rz; k++)
                    {
                        int kk = k <= rz / 2 ? k : rz - k;
                        int shell = Mathf.RoundToInt(Mathf.Sqrt(ki * ki + kj * kj + kk * kk));
                        if (shell > kMax) continue;
                        int idx = i * ry * rz + j * rz + k;
                        double u = vx[idx], v = vy[idx], w = vz[idx];
                        eShell[shell] += 0.5 * (u * u + v * v + w * w);
                        countShell[shell]++;
                    }
                }
            }

            // Нормировать на число ячеек и объём — чтобы сравнение ISF/LES было честным
            var csv = writer.CreateCsv($"spectrum_{tag}.csv", "k_shell", "E_k", "cell_count");
            for (int k = 1; k <= kMax; k++)
            {
                double ek = countShell[k] > 0 ? eShell[k] * dV / countShell[k] : 0.0;
                csv.Row(k, (float)ek, countShell[k]);
            }
        }

        private void SaveSlice(TestArtifactWriter writer, IsfTestHarness h, string tag, int index)
        {
            var tex = IsfDiagnostics.RenderVorticitySlice(h.Isf, h.Velocity, SliceAxis.Z, h.Isf.resZ / 2);
            try { writer.Png($"vorticity_{tag}_{index:D2}.png", tex); }
            finally { UnityEngine.Object.Destroy(tex); }
        }

        private static float SafeDiv(float a, float b) =>
            Mathf.Abs(b) < 1e-20f ? 0f : a / b;

        [Serializable]
        private struct SnapshotSer
        {
            public float kineticEnergy, maxOmega, rmsOmega, maxVel, rmsVel;
        }
        private static SnapshotSer ToSer(IsfDiagnosticsSnapshot s) => new SnapshotSer
        {
            kineticEnergy = s.KineticEnergy,
            maxOmega = s.MaxVorticity,
            rmsOmega = s.RmsVorticity,
            maxVel = s.MaxVelocity,
            rmsVel = s.RmsVelocity
        };

        [Serializable]
        private struct Summary
        {
            public string testName;
            public int[] volRes;
            public float hbar, dt, csSmagorinsky, kinematicViscosity;
            public int totalSteps;
            public SnapshotSer noLesFinal;
            public SnapshotSer lesFinal;
        }
    }
}
