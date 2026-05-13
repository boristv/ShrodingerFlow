# Аналитические тесты ISF / ShrodingerFlow

Набор «проверяемых численных тестов» (см. план в чате). Не Unity Test Framework — это
обычные `MonoBehaviour`-сценарии, которые запускаются из Play-режима и сохраняют
артефакты (`metrics.csv`, `summary.json`, `vorticity_slice_*.png`, ...) в `AnalyticalTestResults/`.

## Структура

```
Assets/Tests/Analytical/
  Runtime/
    AnalyticalTestRunner.cs       // диспетчер: собирает включённые тесты, выдаёт каталог под артефакты
    AnalyticalTestBase.cs         // абстрактный тест
    TestRunContext.cs             // общие compute-шейдеры + outputDir
    TestArtifactWriter.cs         // CSV / JSON / PNG / лог
    IsfTestHarness.cs             // тонкая обёртка над CSISF/CSVelocity, init ψ + AddCircle
    IsfDiagnostics.cs             // E_kin, |ω|, ∫ρ, центроид, viridis-срезы
    Tests/
      VortexRingConservationTest.cs   // тест 1: сохранение вихревых структур
      IsfVsLesTest.cs                 // тест 2: ISF vs ISF+LES
      HbarSweepTest.cs                // тест 3: влияние ℏ
      DensityTransportTest.cs         // тест 4: перенос α = |ψ|²
      OpticalConsistencyTest.cs       // тест 5: оптическая согласованность
      PerformanceBenchmarkTest.cs     // тест 6: производительность GPU
  Editor/
    AnalyticalTestsSceneSetup.cs  // меню Create/Open Scene и Open Results Folder
  Scenes/
    AnalyticalTests.unity         // одна сцена, к которой прикручены все тесты-компоненты
```

## Как запустить

1. `ShrodingerFlow → Analytical Tests → Create or Update Scene`. Если уже есть — `Open Scene`.
2. Выделите `AnalyticalTests` GameObject. На нём:
   - в `AnalyticalTestRunner` — назначены compute-шейдеры (`SFComputeKernels`, `SFComputeFFT`, `SFComputeLES`), путь `AnalyticalTestResults`, флаги авто-запуска и выхода из Play.
   - один или несколько компонентов `AnalyticalTestBase` (например `VortexRingConservationTest`). У каждого свой флажок «Enabled In Batch» и набор параметров.
3. Нажмите Play.
4. Артефакты появятся в `AnalyticalTestResults/<TestName>/<yyyy-MM-dd_HH-mm-ss>/`. В корне ещё `session_<runId>.json` — общий список тестов и времени.
5. `ShrodingerFlow → Analytical Tests → Open Results Folder` открывает каталог результатов в системном файловом менеджере.

> ⚠ Папка `AnalyticalTestResults/` добавлена в `.gitignore`, поэтому сами артефакты не коммитятся.

## Что записывают тесты по умолчанию

### `metrics.csv`

Одна строка на семпл (`_sampleEveryNSteps`). В Тесте 1 колонки:

| колонка | смысл |
| --- | --- |
| `step`, `t` | номер шага и физическое время `step·dt` |
| `kinetic_energy`, `kinetic_energy_rel` | E(t) = ½∑|u|²·dV и E(t)/E(0) |
| `density_integral`, `density_integral_rel` | ∫ρ dV ≈ vol — проверка консервативности `Normalize` |
| `max_omega`, `rms_omega`, `enstrophy` | резкость и общая интенсивность вихрей |
| `max_vel`, `rms_vel` | контроль «разгона» поля скорости |
| `centroid_*`, `centroid_shift` | центроид вихревого облака (трекинг кольца) |

### `summary.json`

Снимки `initial` / `final`, относительные изменения, `wallTimeSeconds`, `avgStepMillis`,
параметры запуска (`hbar`, `dt`, `volRes`, `useLES`, `kinematicViscosity`...).

### `vorticity_slice_*.png`

Срез |ω| в XY-плоскости через середину объёма (viridis-палитра). По умолчанию 4 кадра
равномерно по прогрессу — хватает для иллюстраций «3–4 кадра эволюции» из плана.

## Артефакты по тестам

| Тест | CSV | PNG | summary.json |
| --- | --- | --- | --- |
| 1 VortexRingConservation | `metrics.csv` — E(t), ∫ρ, \|ω\|, центроид | `vorticity_slice_z_NN.png` × 4 | params + initial/final |
| 2 IsfVsLes | `metrics_noLES.csv`, `metrics_LES.csv`, `comparison.csv`, `spectrum_*.csv` | `vorticity_noLES_*.png` + `vorticity_LES_*.png` | params + finals |
| 3 HbarSweep | `metrics_hbar_*.csv`, `spectrum_hbar_*.csv`, `sweep_summary.csv` | `vorticity_hbar_*_NN.png` | sweep по ℏ |
| 4 DensityTransport | `metrics.csv` — ∫ρ, центр масс, ширина, расхождение с вихрем | `density_slice_NN.png` × 4 | mass conservation |
| 5 OpticalConsistency | `beer_lambert.csv` — τ/T по лучам, `normals_sample.csv` | `density_slice_Z.png`, `density_slice_Y.png`, `normal_slice.png` | pass/fail + пороги |
| 6 PerformanceBenchmark | `perf_results.csv` — ms/step, FPS, mem, `per_step_*.csv` | — | entries по сетке |

## Расширение

Чтобы добавить новый тест:

1. Создайте `Assets/Tests/Analytical/Runtime/Tests/MyTest.cs`, унаследуйтесь от `AnalyticalTestBase`.
2. Реализуйте `TestName`, `Description` и корутину `RunTest(TestRunContext ctx)`.
3. Внутри: создайте `IsfTestHarness`, прогоните шаги, пишите через `TestArtifactWriter`.
4. Добавьте компонент на `AnalyticalTests` GameObject в сцене (или обновите Editor-меню,
   чтобы оно добавляло компонент при создании сцены).

Конвенция приватных полей — `_camelCase` (как в остальном проекте).
