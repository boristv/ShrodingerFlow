using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ShrodingerFlow.AnalyticalTests
{
    /// <summary>
    /// Тест 1 «Сохранение вихревых структур в ISF».
    /// Кладёт одно или несколько вихревых колец, прогоняет N шагов, по пути замеряет:
    ///   • кинетическую энергию E(t);
    ///   • интеграл плотности ∫ρ dV (должен ≈ vol);
    ///   • max|ω| и RMS|ω| (мера резкости вихря и общего «деграда»);
    ///   • энстрофию ∫|ω|² dV;
    ///   • центроид вихревого облака (трекинг кольца).
    /// Дополнительно сохраняет PNG-срезы |ω| на 4 ключевых моментах.
    /// </summary>
    public sealed class VortexRingConservationTest : AnalyticalTestBase
    {
        [Header("Сетка")]
        [SerializeField] private Vector3Int _volSize = new Vector3Int(10, 5, 5);
        [SerializeField] private Vector3Int _volRes = new Vector3Int(128, 64, 64);
        [SerializeField] private float _hbar = 0.1f;
        [SerializeField] private float _dt = 1f / 12f;

        [Header("Фоновый поток (для exp(i k·x))")]
        [SerializeField] private Vector3 _backgroundFlow = new Vector3(-0.2f, 0f, 0f);
        [Tooltip("Если включено — фоновый поток через плоскую волну. Иначе ψ₁=1, ψ₂=psi2Init.")]
        [SerializeField] private bool _useBackgroundPlaneWave = true;
        [SerializeField] private Vector2 _psi2Init = new Vector2(0.01f, 0f);

        [Header("Кольца (1 — одиночный вихрь, 2 — leapfrog)")]
        [SerializeField]
        private List<VortexRingDef> _rings = new List<VortexRingDef>
        {
            new VortexRingDef { Center = new Vector3(5f, 2.5f, 2.5f), Normal = new Vector3(-1f, 0f, 0f), Radius = 1.5f, Thickness = 0f },
            new VortexRingDef { Center = new Vector3(5f, 2.5f, 2.5f), Normal = new Vector3(-1f, 0f, 0f), Radius = 0.9f, Thickness = 0f },
        };

        [Header("Симуляция")]
        [Tooltip("Сколько шагов ISF выполнить.")]
        [SerializeField] private int _totalSteps = 600;
        [Tooltip("Каждые сколько шагов писать строку в metrics.csv.")]
        [SerializeField] private int _sampleEveryNSteps = 5;
        [Tooltip("Сколько срезов |ω| сохранить как PNG (равномерно по прогрессу).")]
        [Range(0, 16)]
        [SerializeField] private int _slicePngCount = 4;

        [Header("Опции шага")]
        [SerializeField] private bool _useLES;
        [SerializeField] private float _kinematicViscosity;
        [SerializeField] private float _csSmagorinsky = 0.07f;

        public override string TestName => "01_VortexRingConservation";
        public override string Description =>
            $"Эволюция {_rings?.Count ?? 0} вихревого(ых) кольца(колец), {_totalSteps} шагов, LES={_useLES}, ν={_kinematicViscosity}.";

        public override IEnumerator RunTest(TestRunContext ctx)
        {
            if (_rings == null || _rings.Count == 0)
            {
                Debug.LogWarning($"[{TestName}] Список колец пустой — тест пропущен.");
                yield break;
            }

            using (var writer = new TestArtifactWriter(ctx.OutputDir))
            {
                int[] volSize = { _volSize.x, _volSize.y, _volSize.z };
                int[] volRes = { _volRes.x, _volRes.y, _volRes.z };

                writer.Log($"vol_size = {_volSize}, vol_res = {_volRes}, hbar = {_hbar}, dt = {_dt}");

                using (var h = new IsfTestHarness(ctx.KernelsShader, ctx.FftShader, ctx.LesShader,
                    volSize, volRes, _hbar, _dt))
                {
                    var ringsArr = _rings.ToArray();
                    if (_useBackgroundPlaneWave)
                        h.InitPlaneWave(_backgroundFlow, _psi2Init.x);
                    else
                        h.InitUniformPsi(_psi2Init);

                    h.AddVortexRings(ringsArr);
                    h.NormalizeAndProject();
                    h.Isf.UpdateVelocities(h.Velocity);

                    var initial = IsfDiagnostics.SampleAll(h.Isf, h.Velocity);
                    writer.Log($"Initial: E={initial.KineticEnergy:G6}, ∫ρdV={initial.DensityIntegral:G6}, max|ω|={initial.MaxVorticity:G6}, Ω={initial.Enstrophy:G6}");

                    var csv = writer.CreateCsv("metrics.csv",
                        "step", "t",
                        "kinetic_energy", "kinetic_energy_rel",
                        "density_integral", "density_integral_rel",
                        "max_omega", "rms_omega", "enstrophy",
                        "max_vel", "rms_vel",
                        "centroid_x", "centroid_y", "centroid_z",
                        "centroid_shift");

                    Vector3 initialCentroid = initial.VorticityCentroid;
                    LogCsvRow(csv, 0, 0f, initial, initial, initialCentroid);
                    SaveSlice(writer, h, 0);

                    var opts = new StepOptions
                    {
                        UseLES = _useLES,
                        KinematicViscosity = _kinematicViscosity,
                        CsSmagorinsky = _csSmagorinsky,
                        GravityPsi2 = null
                    };

                    int nextSliceStep = _slicePngCount > 1
                        ? Mathf.Max(1, _totalSteps / Mathf.Max(1, _slicePngCount - 1))
                        : int.MaxValue;
                    int sliceIndex = 1;

                    ctx.ReportProgress(0, _totalSteps);

                    double simStart = Time.realtimeSinceStartupAsDouble;
                    for (int step = 1; step <= _totalSteps; step++)
                    {
                        h.Step(opts);

                        bool isSampleStep = (step % _sampleEveryNSteps) == 0 || step == _totalSteps;
                        if (isSampleStep)
                        {
                            var s = IsfDiagnostics.SampleAll(h.Isf, h.Velocity);
                            LogCsvRow(csv, step, h.SimTime, s, initial, initialCentroid);

                            string status = $"E={s.KineticEnergy:G4}  max|ω|={s.MaxVorticity:G4}  drift={( s.VorticityCentroid - initialCentroid).magnitude:F3}";
                            ctx.ReportProgress(step, _totalSteps, status);
                            yield return null;
                        }

                        if (sliceIndex < _slicePngCount && step >= nextSliceStep)
                        {
                            SaveSlice(writer, h, sliceIndex);
                            sliceIndex++;
                            nextSliceStep = _slicePngCount > 1
                                ? sliceIndex * (_totalSteps / Mathf.Max(1, _slicePngCount - 1))
                                : int.MaxValue;

                            // Обновить live-preview последним срезом |ω|
                            var preview = IsfDiagnostics.RenderVorticitySlice(h.Isf, h.Velocity,
                                SliceAxis.Z, h.Isf.resZ / 2);
                            ctx.SetPreview(preview);
                            // Не уничтожаем: раннер держит ссылку до следующего preview
                        }
                    }
                    double simWallTime = Time.realtimeSinceStartupAsDouble - simStart;

                    if (sliceIndex < _slicePngCount)
                        SaveSlice(writer, h, sliceIndex);

                    var final = IsfDiagnostics.SampleAll(h.Isf, h.Velocity);

                    var summary = new Summary
                    {
                        testName = TestName,
                        volSize = new[] { _volSize.x, _volSize.y, _volSize.z },
                        volRes = new[] { h.Isf.resX, h.Isf.resY, h.Isf.resZ },
                        hbar = _hbar,
                        dt = _dt,
                        totalSteps = _totalSteps,
                        useLES = _useLES,
                        kinematicViscosity = _kinematicViscosity,
                        csSmagorinsky = _csSmagorinsky,
                        ringCount = ringsArr.Length,
                        wallTimeSeconds = (float)simWallTime,
                        avgStepMillis = (float)(simWallTime * 1000.0 / Mathf.Max(1, _totalSteps)),
                        initial = ToSerializable(initial),
                        final = ToSerializable(final),
                        relEnergy = SafeRel(final.KineticEnergy, initial.KineticEnergy),
                        relEnstrophy = SafeRel(final.Enstrophy, initial.Enstrophy),
                        relMaxOmega = SafeRel(final.MaxVorticity, initial.MaxVorticity),
                        relDensityIntegral = SafeRel(final.DensityIntegral, initial.DensityIntegral),
                        centroidDrift = (final.VorticityCentroid - initialCentroid).magnitude,
                    };
                    writer.Json("summary.json", summary);
                    writer.Log($"Final: E={final.KineticEnergy:G6} (rel {summary.relEnergy:F4}), max|ω|={final.MaxVorticity:G6} (rel {summary.relMaxOmega:F4}), drift={summary.centroidDrift:F4}");
                    writer.Log($"Wall time {simWallTime:F2} c, {summary.avgStepMillis:F2} ms/шаг (включая GetData диагностики).");
                }
            }
        }

        private static void LogCsvRow(CsvLog csv, int step, float t,
            IsfDiagnosticsSnapshot s, IsfDiagnosticsSnapshot initial, Vector3 centroid0)
        {
            csv.Row(
                step, t,
                s.KineticEnergy, SafeRel(s.KineticEnergy, initial.KineticEnergy),
                s.DensityIntegral, SafeRel(s.DensityIntegral, initial.DensityIntegral),
                s.MaxVorticity, s.RmsVorticity, s.Enstrophy,
                s.MaxVelocity, s.RmsVelocity,
                s.VorticityCentroid.x, s.VorticityCentroid.y, s.VorticityCentroid.z,
                (s.VorticityCentroid - centroid0).magnitude);
        }

        private void SaveSlice(TestArtifactWriter writer, IsfTestHarness h, int index)
        {
            int sliceIdx = h.Isf.resZ / 2;
            var tex = IsfDiagnostics.RenderVorticitySlice(h.Isf, h.Velocity, SliceAxis.Z, sliceIdx);
            try
            {
                writer.Png($"vorticity_slice_z_{index:D2}_step{Mathf.RoundToInt(h.SimTime / h.Dt):D5}.png", tex);
            }
            finally
            {
                UnityEngine.Object.Destroy(tex);
            }
        }

        private static float SafeRel(float current, float baseline) =>
            Mathf.Abs(baseline) < 1e-20f ? 0f : current / baseline;

        [Serializable]
        private struct SnapshotSer
        {
            public float kineticEnergy;
            public float densityIntegral;
            public float maxVorticity;
            public float rmsVorticity;
            public float enstrophy;
            public float maxVelocity;
            public float rmsVelocity;
            public Vector3 vorticityCentroid;
        }

        private static SnapshotSer ToSerializable(IsfDiagnosticsSnapshot s) => new SnapshotSer
        {
            kineticEnergy = s.KineticEnergy,
            densityIntegral = s.DensityIntegral,
            maxVorticity = s.MaxVorticity,
            rmsVorticity = s.RmsVorticity,
            enstrophy = s.Enstrophy,
            maxVelocity = s.MaxVelocity,
            rmsVelocity = s.RmsVelocity,
            vorticityCentroid = s.VorticityCentroid
        };

        [Serializable]
        private struct Summary
        {
            public string testName;
            public int[] volSize;
            public int[] volRes;
            public float hbar;
            public float dt;
            public int totalSteps;
            public bool useLES;
            public float kinematicViscosity;
            public float csSmagorinsky;
            public int ringCount;
            public float wallTimeSeconds;
            public float avgStepMillis;
            public SnapshotSer initial;
            public SnapshotSer final;
            public float relEnergy;
            public float relEnstrophy;
            public float relMaxOmega;
            public float relDensityIntegral;
            public float centroidDrift;
        }
    }
}
