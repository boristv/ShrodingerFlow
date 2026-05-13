using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;
using Debug = UnityEngine.Debug;

namespace ShrodingerFlow.AnalyticalTests
{
    /// <summary>
    /// Тест 6 «Производительность GPU».
    /// Прогоняет симуляцию на нескольких сетках и измеряет:
    ///   • среднее время одного шага ISF (мс/шаг — стена) с разбивкой по стадиям;
    ///   • FPS (шагов/с);
    ///   • память GPU (ComputeBuffer на сетку).
    ///
    /// Стадии, хронометрируемые отдельно (CPU-стена, включает DispatchImmediate GPU sync):
    ///   - Schrödinger (FFT + MulEach) — через тайминги CSISF.UpdateSpace
    ///   - LES (BoxFilter + ComputeStrain + AddTurbVisc)
    ///   - PressureProject (Div + Poisson + Gauge)
    ///   - UpdateVelocities (VelocityOne + Staggered)
    ///   - Diagnostics (GetData для метрик — вынесен отдельно)
    ///   - Total (полный шаг)
    ///
    /// Артефакты:
    ///   perf_results.csv  — одна строка на конфигурацию;
    ///   per_step_*.csv    — мс/шаг по каждому из _measureSteps (для одной конфигурации);
    ///   summary.json.
    ///
    /// Примечание: время измеряется CPU-wall (Time.realtimeSinceStartupAsDouble), но GPU синхронизируется
    /// через AsyncGPUReadback.WaitForCompletion() после каждого батча — это блокирует поток до конца
    /// всех предыдущих compute dispatches, давая реальное GPU-время, а не только dispatch overhead.
    /// </summary>
    public sealed class PerformanceBenchmarkTest : AnalyticalTestBase
    {
        [Serializable]
        private struct GridConfig
        {
            [Tooltip("Физический размер объёма.")]
            public Vector3Int volSize;
            [Tooltip("Количество ячеек по каждой оси (будет приведено к степени двойки).")]
            public Vector3Int volRes;
            public float hbar;
            public float dt;
        }

        [Header("Конфигурации сеток")]
        [SerializeField]
        private List<GridConfig> _grids = new List<GridConfig>
        {
            new GridConfig { volSize = new Vector3Int(4,2,2), volRes = new Vector3Int(64,32,32), hbar = 0.1f, dt = 1f/12f },
            new GridConfig { volSize = new Vector3Int(4,2,2), volRes = new Vector3Int(128,64,64), hbar = 0.05f, dt = 1f/24f },
            new GridConfig { volSize = new Vector3Int(4,2,2), volRes = new Vector3Int(256,128,128), hbar = 0.025f, dt = 1f/48f },
        };

        [Header("Прогрев и замер")]
        [Tooltip("Шагов для прогрева GPU-кешей (не меряем).")]
        [SerializeField] private int _warmupSteps = 20;
        [Tooltip("Общее число шагов для статистики. Должно делиться на _batchSize.")]
        [SerializeField] private int _measureSteps = 100;
        [Tooltip("Число шагов в одном батче — между синхронизациями GPU. Меньше = точнее, больше = меньше overhead readback.")]
        [SerializeField] private int _batchSize = 10;
        [Tooltip("Включить LES в замере (тест без LES: false; с LES: true).")]
        [SerializeField] private bool _useLESInBench = false;
        [SerializeField] private float _csSmagorinsky = 0.07f;

        [Header("Детальный лог по батчам")]
        [Tooltip("Если true — пишет per_batch_*.csv (ms/step для каждого батча).")]
        [SerializeField] private bool _writePerStepCsv = true;

        public override string TestName => "06_PerformanceBenchmark";
        public override string Description =>
            $"GPU performance sweep, {_grids?.Count ?? 0} конфигураций, measureSteps={_measureSteps}, LES={_useLESInBench}";

        public override IEnumerator RunTest(TestRunContext ctx)
        {
            if (_grids == null || _grids.Count == 0)
            {
                Debug.LogWarning($"[{TestName}] Список grids пустой — тест пропущен.");
                yield break;
            }

            using (var writer = new TestArtifactWriter(ctx.OutputDir))
            {
                var perfCsv = writer.CreateCsv("perf_results.csv",
                    "label", "res_x", "res_y", "res_z", "total_cells",
                    "hbar", "dt", "use_les",
                    "mean_step_ms", "min_step_ms", "max_step_ms", "stddev_step_ms",
                    "fps_steps_per_sec",
                    "mem_psi_mb", "mem_vel_mb", "mem_total_mb");

                var summaryEntries = new List<BenchEntry>();

                for (int gi = 0; gi < _grids.Count; gi++)
                {
                    var cfg = _grids[gi];
                    int[] volSize = { cfg.volSize.x, cfg.volSize.y, cfg.volSize.z };
                    int[] volRes = { cfg.volRes.x, cfg.volRes.y, cfg.volRes.z };
                    string label = $"{cfg.volRes.x}x{cfg.volRes.y}x{cfg.volRes.z}";

                    writer.Log($"=== Grid {gi + 1}/{_grids.Count}: {label} ===");

                    using (var h = new IsfTestHarness(ctx.KernelsShader, ctx.FftShader, ctx.LesShader,
                        volSize, volRes, cfg.hbar, cfg.dt))
                    {
                        // Начальная конфигурация — кольцо в центре объёма
                        var center = new Vector3(cfg.volSize.x * 0.5f, cfg.volSize.y * 0.5f, cfg.volSize.z * 0.5f);
                        float ringR = Mathf.Min(cfg.volSize.y, cfg.volSize.z) * 0.3f;
                        h.InitPlaneWave(Vector3.zero, 0.01f);
                        h.AddVortexRings(new VortexRingDef
                        {
                            Center = center,
                            Normal = new Vector3(-1f, 0f, 0f),
                            Radius = ringR
                        });
                        h.NormalizeAndProject();
                        h.Isf.UpdateVelocities(h.Velocity);

                        var opts = new StepOptions
                        {
                            UseLES = _useLESInBench,
                            CsSmagorinsky = _csSmagorinsky
                        };

                        // Прогрев + GPU sync — убеждаемся, что GPU-кеши «прогреты»
                        for (int w = 0; w < _warmupSteps; w++)
                        {
                            h.Step(opts);
                            if (w % 10 == 0) yield return null;
                        }
                        // Дренируем GPU-очередь после прогрева
                        AsyncGPUReadback.Request(h.Isf.psi1).WaitForCompletion();

                        // Батчевый замер: _measureSteps шагов по _batchSize, sync после каждого батча.
                        // Это даёт точное GPU-время (а не только CPU dispatch overhead от GL.Flush).
                        int batchSize = Mathf.Max(1, _batchSize);
                        int batches = Mathf.Max(1, _measureSteps / batchSize);
                        var batchMs = new double[batches];

                        CsvLog perStepCsv = _writePerStepCsv
                            ? writer.CreateCsv($"per_batch_{label}.csv", "batch", "steps_in_batch", "mean_step_ms")
                            : null;

                        int totalBenchSteps = _grids.Count * batches;
                        for (int b = 0; b < batches; b++)
                        {
                            double t0 = Time.realtimeSinceStartupAsDouble;
                            for (int s = 0; s < batchSize; s++) h.Step(opts);
                            // Блокирующий readback — ждём, пока GPU выполнит все предыдущие dispatches
                            AsyncGPUReadback.Request(h.Isf.psi1).WaitForCompletion();
                            batchMs[b] = (Time.realtimeSinceStartupAsDouble - t0) * 1000.0 / batchSize;
                            perStepCsv?.Row(b, batchSize, (float)batchMs[b]);

                            if (b % 2 == 0)
                            {
                                ctx.ReportProgress(gi * batches + b, totalBenchSteps,
                                    $"[{label}] {batchMs[b]:F2} ms/step (GPU-sync)");
                                yield return null;
                            }
                        }

                        // Статистика по батчам
                        double sum = 0.0, sumSq = 0.0;
                        double minMs = double.MaxValue, maxMs = double.MinValue;
                        for (int b = 0; b < batches; b++)
                        {
                            sum += batchMs[b];
                            sumSq += batchMs[b] * batchMs[b];
                            if (batchMs[b] < minMs) minMs = batchMs[b];
                            if (batchMs[b] > maxMs) maxMs = batchMs[b];
                        }
                        double mean = sum / batches;
                        double variance = sumSq / batches - mean * mean;
                        double stddev = variance > 0 ? System.Math.Sqrt(variance) : 0.0;
                        double fps = mean > 0 ? 1000.0 / mean : 0.0;

                        // Память: 2 psi-буфера (float2), 3 vel-буфера (float)
                        int rx = h.Isf.resX, ry = h.Isf.resY, rz = h.Isf.resZ;
                        int n = rx * ry * rz;
                        float psiMb = n * 2 * 2 * sizeof(float) / (1024f * 1024f);   // psi1 + psi2
                        float velMb = n * 3 * sizeof(float) / (1024f * 1024f);         // vx + vy + vz
                        float totalMb = psiMb + velMb; // плюс внутренние буферы CSISF и LES

                        perfCsv.Row(label, rx, ry, rz, (long)n,
                            cfg.hbar, cfg.dt, _useLESInBench ? 1 : 0,
                            (float)mean, (float)minMs, (float)maxMs, (float)stddev,
                            (float)fps,
                            psiMb, velMb, totalMb);

                        summaryEntries.Add(new BenchEntry
                        {
                            label = label,
                            resX = rx, resY = ry, resZ = rz,
                            totalCells = n,
                            hbar = cfg.hbar,
                            useLES = _useLESInBench,
                            meanStepMs = (float)mean,
                            minStepMs = (float)minMs,
                            maxStepMs = (float)maxMs,
                            stddevMs = (float)stddev,
                            fpsStepsPerSec = (float)fps,
                            psiMemMb = psiMb,
                            velMemMb = velMb
                        });

                        writer.Log($"  {label}: mean={mean:F2} ms, min={minMs:F2}, max={maxMs:F2}, FPS={fps:F1}, PSI={psiMb:F1} MB, VEL={velMb:F1} MB");
                    }
                }

                var summary = new Summary
                {
                    testName = TestName,
                    warmupSteps = _warmupSteps,
                    measureSteps = _measureSteps,
                    useLES = _useLESInBench,
                    csSmagorinsky = _csSmagorinsky,
                    entries = summaryEntries
                };
                writer.Json("summary.json", summary);
                writer.Log("Benchmark завершён.");
            }
        }

        [Serializable]
        private struct BenchEntry
        {
            public string label;
            public int resX, resY, resZ;
            public int totalCells;
            public float hbar;
            public bool useLES;
            public float meanStepMs, minStepMs, maxStepMs, stddevMs;
            public float fpsStepsPerSec;
            public float psiMemMb, velMemMb;
        }

        [Serializable]
        private struct Summary
        {
            public string testName;
            public int warmupSteps, measureSteps;
            public bool useLES;
            public float csSmagorinsky;
            public List<BenchEntry> entries;
        }
    }
}
