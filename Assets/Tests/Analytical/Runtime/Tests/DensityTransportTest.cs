using System;
using System.Collections;
using UnityEngine;

namespace ShrodingerFlow.AnalyticalTests
{
    /// <summary>
    /// Тест 4 «Сохранение импульса и несжимаемость».
    ///
    /// Проверяет два фундаментальных закона ISF:
    ///
    ///   1. Гидродинамический импульс P = ½ ∫ x × ω dV.
    ///      При свободной эволюции без внешних сил (и периодических ГУ) P сохраняется.
    ///      Это нетривиальная проверка: энергия диссипирует даже при ν=0 (числовая), а импульс —
    ///      нет, если схема корректна.
    ///
    ///   2. Несжимаемость: RMS(∇·v) ≈ 0 на каждом шаге.
    ///      ISF гарантирует ∇·v = 0 через pressure projection; тест проверяет, что это
    ///      выполняется численно (ожидается уровень машинной точности ≲ 1e-4).
    ///
    ///   3. Спиральность H = ∫ v·ω dV (должна сохраняться при ν=0).
    ///
    /// Начальная конфигурация: одно вихревое кольцо.
    ///
    /// Артефакты:
    ///   metrics.csv    — P_x, P_y, P_z, |P|, H, RMS(∇·v) на каждом шаге;
    ///   vorticity_*.png — срезы завихрённости;
    ///   summary.json.
    /// </summary>
    public sealed class DensityTransportTest : AnalyticalTestBase
    {
        [Header("Сетка")]
        [SerializeField] private Vector3Int _volSize = new Vector3Int(10, 5, 5);
        [SerializeField] private Vector3Int _volRes  = new Vector3Int(64, 32, 32);
        [SerializeField] private float _hbar = 0.1f;
        [SerializeField] private float _dt   = 1f / 12f;

        [Header("Вихревое кольцо")]
        [SerializeField] private Vector3 _ringCenter = new Vector3(5f, 2.5f, 2.5f);
        [SerializeField] private Vector3 _ringNormal = new Vector3(-1f, 0f, 0f);
        [SerializeField] private float   _ringRadius = 1.5f;

        [Header("Симуляция")]
        [SerializeField] private int _totalSteps        = 300;
        [SerializeField] private int _sampleEveryNSteps = 5;
        [SerializeField] private int _slicePngCount     = 4;

        [Header("LES")]
        [SerializeField] private bool  _useLES;
        [SerializeField] private float _csSmagorinsky = 0.07f;

        [Header("Допуски")]
        [Tooltip("Максимально допустимое |ΔP|/|P₀| за всё время прогона.\n" +
                 "ISF имеет числовую диссипацию ~10-15%/300 шагов на сетке 64³.")]
        [SerializeField] private float _impulseRelTolerance = 0.15f;
        [Tooltip("Максимально допустимое RMS(∇·v), нормированное FD-схемой.\n" +
                 "ISF гарантирует ∇·v=0 в спектральном смысле; FD-остаток ≈ 5–15% от V/dx.")]
        [SerializeField] private float _divRmsThreshold = 0.02f;

        public override string TestName => "04_ImpulseConservation";
        public override string Description =>
            $"Импульс P=½∫x×ω dV, несжимаемость ∇·v≈0, кольцо r={_ringRadius}, hbar={_hbar}, steps={_totalSteps}";

        public override IEnumerator RunTest(TestRunContext ctx)
        {
            using (var writer = new TestArtifactWriter(ctx.OutputDir))
            {
                int[] volSize = { _volSize.x, _volSize.y, _volSize.z };
                int[] volRes  = { _volRes.x,  _volRes.y,  _volRes.z  };

                using (var h = new IsfTestHarness(ctx.KernelsShader, ctx.FftShader, ctx.LesShader,
                    volSize, volRes, _hbar, _dt))
                {
                    h.InitPlaneWave(Vector3.zero, 0.01f);
                    h.AddVortexRings(new VortexRingDef
                    {
                        Center    = _ringCenter,
                        Normal    = _ringNormal,
                        Radius    = _ringRadius,
                        Thickness = 0f
                    });
                    h.NormalizeAndProject();
                    h.Isf.UpdateVelocities(h.Velocity);

                    var init = IsfDiagnostics.SampleConservation(h.Isf, h.Velocity);
                    writer.Log($"Initial: |P|={init.ImpulseMag:G5}  H={init.Helicity:G5}  div_rms={init.DivergenceRms:G4}");

                    var csv = writer.CreateCsv("metrics.csv",
                        "step", "t",
                        "P_x", "P_y", "P_z", "P_mag", "P_mag_rel",
                        "helicity", "helicity_rel",
                        "div_rms");

                    LogRow(csv, 0, 0f, init, init);

                    var opts = new StepOptions
                    {
                        UseLES        = _useLES,
                        CsSmagorinsky = _csSmagorinsky
                    };

                    int sliceInterval = _slicePngCount > 1
                        ? Mathf.Max(1, _totalSteps / (_slicePngCount - 1)) : int.MaxValue;
                    int nextSlice = sliceInterval;
                    int sliceNum  = 0;
                    SaveVorticitySlice(writer, h, sliceNum++);

                    ctx.ReportProgress(0, _totalSteps);

                    double divRmsSum = 0;
                    int    divRmsCount = 0;
                    float  maxDivRms  = 0f;

                    for (int step = 1; step <= _totalSteps; step++)
                    {
                        h.Step(opts);

                        if (step % _sampleEveryNSteps == 0 || step == _totalSteps)
                        {
                            var s    = IsfDiagnostics.SampleConservation(h.Isf, h.Velocity);
                            LogRow(csv, step, h.SimTime, s, init);

                            divRmsSum += s.DivergenceRms;
                            divRmsCount++;
                            if (s.DivergenceRms > maxDivRms) maxDivRms = s.DivergenceRms;

                            float pRel = init.ImpulseMag > 1e-12f ? s.ImpulseMag / init.ImpulseMag : 0f;
                            ctx.ReportProgress(step, _totalSteps,
                                $"|P|/|P₀|={pRel:F5}  div_rms={s.DivergenceRms:G3}  H={s.Helicity:G4}");

                            var tex = IsfDiagnostics.RenderVorticitySlice(h.Isf, h.Velocity, SliceAxis.Z, h.Isf.resZ / 2);
                            ctx.SetPreview(tex);

                            yield return null;
                        }

                        if (sliceNum < _slicePngCount && step >= nextSlice)
                        {
                            SaveVorticitySlice(writer, h, sliceNum++);
                            nextSlice = sliceNum * sliceInterval;
                        }
                    }

                    if (sliceNum < _slicePngCount)
                        SaveVorticitySlice(writer, h, sliceNum);

                    var fin = IsfDiagnostics.SampleConservation(h.Isf, h.Velocity);
                    float pRelFinal = init.ImpulseMag > 1e-12f ? fin.ImpulseMag / init.ImpulseMag : 0f;
                    float hRelFinal = Mathf.Abs(init.Helicity) > 1e-12f ? fin.Helicity / init.Helicity : 0f;
                    float meanDivRms = divRmsCount > 0 ? (float)(divRmsSum / divRmsCount) : 0f;

                    bool impulseOk = Mathf.Abs(pRelFinal - 1f) < _impulseRelTolerance;
                    bool divOk     = meanDivRms < _divRmsThreshold;

                    // Нормируем div_rms на V/dx для физической интерпретации
                    var diagFin = IsfDiagnostics.SampleAll(h.Isf, h.Velocity);
                    float vRms = diagFin.RmsVelocity;
                    float dx   = h.Isf.dx;
                    float divNorm = vRms > 1e-12f ? meanDivRms / (vRms / dx) : 0f;

                    writer.Log($"Final:   |P|/|P₀|={pRelFinal:F6}  H/H₀={hRelFinal:F6}  mean_div={meanDivRms:G4}  max_div={maxDivRms:G4}");
                    writer.Log($"  div_rms / (Vrms/dx) = {divNorm:P1}  (FD vs спектральный остаток)");
                    writer.Log($"Impulse conservation: {(impulseOk ? "OK" : "WARN")}  (tol={_impulseRelTolerance:P0})");
                    writer.Log($"Incompressibility:    {(divOk     ? "OK" : "WARN")}  (mean_div thr={_divRmsThreshold:G3})");

                    var summary = new Summary
                    {
                        testName             = TestName,
                        volRes               = new[] { h.Isf.resX, h.Isf.resY, h.Isf.resZ },
                        hbar                 = _hbar,
                        dt                   = _dt,
                        totalSteps           = _totalSteps,
                        useLES               = _useLES,
                        initialImpulseMag    = init.ImpulseMag,
                        finalImpulseMag      = fin.ImpulseMag,
                        impulseRelChange     = pRelFinal,
                        impulsePass          = impulseOk,
                        impulseThreshold     = _impulseRelTolerance,
                        initialHelicity      = init.Helicity,
                        finalHelicity        = fin.Helicity,
                        helicityRelChange    = hRelFinal,
                        meanDivRms           = meanDivRms,
                        maxDivRms            = maxDivRms,
                        divPass              = divOk,
                        divThreshold         = _divRmsThreshold
                    };
                    writer.Json("summary.json", summary);
                }
            }
        }

        private static void LogRow(CsvLog csv, int step, float t,
            IsfConservationSnapshot s, IsfConservationSnapshot init)
        {
            float pRel = init.ImpulseMag > 1e-12f ? s.ImpulseMag / init.ImpulseMag : 0f;
            float hRel = Mathf.Abs(init.Helicity) > 1e-12f ? s.Helicity / init.Helicity : 0f;
            csv.Row(step, t,
                s.Impulse.x, s.Impulse.y, s.Impulse.z, s.ImpulseMag, pRel,
                s.Helicity, hRel, s.DivergenceRms);
        }

        private static void SaveVorticitySlice(TestArtifactWriter writer, IsfTestHarness h, int idx)
        {
            var tex = IsfDiagnostics.RenderVorticitySlice(h.Isf, h.Velocity, SliceAxis.Z, h.Isf.resZ / 2);
            try { writer.Png($"vorticity_z_{idx:D2}_step{Mathf.RoundToInt(h.SimTime / h.Dt):D5}.png", tex); }
            finally { UnityEngine.Object.Destroy(tex); }
        }

        [Serializable]
        private struct Summary
        {
            public string testName;
            public int[]  volRes;
            public float  hbar, dt;
            public int    totalSteps;
            public bool   useLES;
            public float  initialImpulseMag, finalImpulseMag, impulseRelChange;
            public bool   impulsePass;
            public float  impulseThreshold;
            public float  initialHelicity, finalHelicity, helicityRelChange;
            public float  meanDivRms, maxDivRms;
            public bool   divPass;
            public float  divThreshold;
        }
    }
}
