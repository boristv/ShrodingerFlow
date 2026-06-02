using UnityEngine;

/// <summary>
/// Пресеты для «Apply Scenario Defaults»: правьте .asset в инспекторе (как «снимок» с референсных сцен).
/// Jet — JetExampleCS.unity; Sphere — ObstacleExample (SFObstacle); Cylinder — OtherObstacleExample (OtherSFObstacle);
/// TwoSpheres — дефолты SFUnifiedCS / UnifiedCS.unity;
/// Leapfrog — кольца и симуляция из UnifiedCS.unity, объём 10×5×5 (кольца не влезают в 4×2×2).
/// </summary>
[CreateAssetMenu(fileName = "SFUnifiedScenarioPresets", menuName = "ShrodingerFlow/SF Unified Scenario Presets")]
public class SFUnifiedScenarioPresets : ScriptableObject
{
    public const string DefaultAssetPath =
        "Assets/Scenes/JetExample/ComputeShader/SFUnifiedScenarioPresets.asset";

    [Header("Jet — как example_jet.hip (SFJetCS / Apply Scenario Defaults)")]
    public SFUnifiedJetPreset jet;

    [Header("SphereObstacle — как ObstacleExample.unity (SFObstacle)")]
    public SFUnifiedSphereObstaclePreset sphereObstacle;

    [Header("CylinderObstacle — OtherObstacleExample.unity (OtherSFObstacle)")]
    public SFUnifiedCylinderObstaclePreset cylinderObstacle;

    [Header("TwoSpheres — SFUnifiedCS по умолчанию / UnifiedCS.unity")]
    public SFUnifiedTwoSpheresPreset twoSpheres;

    [Header("LeapfrogRings — UnifiedCS.unity (кольца + spawn), vol 10×5×5")]
    public SFUnifiedLeapfrogRingsPreset leapfrogRings;

    [Header("SmokeMaze2D — 2D-лабиринт: источник слева, 3 перегородки, вытяжка справа")]
    public SFUnifiedSmokeMaze2DPreset smokeMaze2D;

    /// <summary>Встроенные значения (копия референсных сцен на момент добавления).</summary>
    public static SFUnifiedScenarioPresets CreateBuiltIn()
    {
        var s = CreateInstance<SFUnifiedScenarioPresets>();
        s.FillBuiltIn();
        return s;
    }

    public void FillBuiltIn()
    {
        // Как example_jet.hip: hbar, dt, samplediv 128 на домене 4×2×2 → 128×64×64, без LES.
        jet = new SFUnifiedJetPreset
        {
            vol_size = new[] { 4, 2, 2 },
            vol_res = new[] { 128, 64, 64 },
            hbar = 0.02f,
            dt = 1f / 48f,
            velocity = new Vector3(1f, 0f, 0f),
            nozzleCen = new Vector3(0.3f, 0.9656632f, 1.0659939f),
            nozzleLen = 0.5f,
            nozzleRad = 0.3f,
            nParticles = 50,
            particleSize = 0.05f,
            stepsPerFrame = 1,
            useLES = false
        };

        // Как ObstacleExample (SFObstacle): vol/hbar/dt/фон/сопло/боксы.
        // В CUDA cuFFT resX=192 допустим; здесь radix-2 FFT — только степени двойки (192 → симуляция ломалась).
        sphereObstacle = new SFUnifiedSphereObstaclePreset
        {
            vol_size = new[] { 6, 2, 2 },
            vol_res = new[] { 128, 64, 64 },
            hbar = 0.02f,
            dt = 1f / 48f,
            velocity = new Vector3(1f, 0f, 0f),
            obstaclePos1 = new Vector3(1.5f, 1f, 1f),
            obstacleRadius1 = 0.5f,
            nozzleCen = new Vector3(0.3f, 0.966f, 1.066f),
            nozzleLen = 0.5f,
            nozzleRad = 0.4f,
            boxSpawnX = new Vector2(3f, 3f),
            boxSpawnY = new Vector2(0.1f, 1.9f),
            boxSpawnZ = new Vector2(0.5f, 1.9f),
            nParticles = 100,
            particleSize = 0.05f,
            stepsPerFrame = 3,
            useLES = true
        };

        // OtherObstacleExample (OtherSFObstacle): те же vol/hbar/dt/фон/сопло/n, боксы как в инспекторе (Y/Z не как у сферы).
        // В сцене resX=192; radix-2 FFT — 128×64×64.
        cylinderObstacle = new SFUnifiedCylinderObstaclePreset
        {
            vol_size = new[] { 6, 2, 2 },
            vol_res = new[] { 128, 64, 64 },
            hbar = 0.02f,
            dt = 1f / 48f,
            velocity = new Vector3(1f, 0f, 0f),
            obstaclePos1 = new Vector3(1.5f, 1f, 1f),
            obstacleRadius1 = 0.5f,
            nozzleCen = new Vector3(0.3f, 0.966f, 1.066f),
            boxSpawnX = new Vector2(3f, 3f),
            boxSpawnY = new Vector2(0.6f, 1.4f),
            boxSpawnZ = new Vector2(0.1f, 1.9f),
            nParticles = 100,
            particleSize = 0.05f,
            stepsPerFrame = 3,
            useLES = true
        };

        // Как сериализованные поля по умолчанию в SFUnifiedCS (сценарий TwoSpheres).
        twoSpheres = new SFUnifiedTwoSpheresPreset
        {
            vol_size = new[] { 4, 2, 2 },
            vol_res = new[] { 64, 32, 32 },
            hbar = 0.1f,
            dt = 1f / 12f,
            velocity = new Vector3(-0.2f, 0f, 0f),
            obstaclePos1 = new Vector3(1.5f, 1f, 1f),
            obstacleRadius1 = 0.5f,
            obstaclePos2 = new Vector3(2.5f, 1f, 1f),
            obstacleRadius2 = 0.5f,
            boxSpawnX = new Vector2(0.3f, 0.3f),
            boxSpawnY = new Vector2(0.5f, 1.5f),
            boxSpawnZ = new Vector2(0.5f, 1.5f),
            nParticles = 50,
            particleSize = 0.1f,
            stepsPerFrame = 3,
            useLES = false
        };

        // UnifiedCS.unity: кольца, hbar/dt/velocity/_nParticles; vol 10×5×5 — иначе радиусы колец не помещаются в домен.
        leapfrogRings = new SFUnifiedLeapfrogRingsPreset
        {
            vol_size = new[] { 10, 5, 5 },
            vol_res = new[] { 128, 64, 64 },
            hbar = 0.1f,
            dt = 1f / 12f,
            velocity = new Vector3(-0.2f, 0f, 0f),
            ring1Radius = 1.5f,
            ring2Radius = 0.9f,
            ring1Normal = new Vector3(-1f, 0f, 0f),
            ring2Normal = new Vector3(-1f, 0f, 0f),
            boxSpawnX = new Vector2(3f, 7f),
            boxSpawnY = new Vector2(0.5f, 4.5f),
            boxSpawnZ = new Vector2(0.5f, 4.5f),
            nParticles = 100000,
            particleSize = 0.1f,
            stepsPerFrame = 3,
            useLES = false
        };

        // 2D-лабиринт 5×1×3 (вид сверху XZ): источник слева, 3 перегородки со смещёнными проходами, вытяжка справа сверху.
        smokeMaze2D = new SFUnifiedSmokeMaze2DPreset
        {
            vol_size = new[] { 5, 1, 3 },
            vol_res = new[] { 128, 32, 64 },
            hbar = 0.05f,
            dt = 1f / 24f,
            wallThickness = 0.12f,
            wallMargin = 0.1f,
            wall1X = 1.0f,
            wall1GapZ = new Vector2(1.0f, 2.0f),
            wall2X = 2.5f,
            wall2SolidZ = new Vector2(0.7f, 2.3f),
            wall3X = 3.8f,
            wall3SolidMaxZ = 2.0f,
            sourceCenter = new Vector3(0.25f, 0.5f, 1.5f),
            sourceHalf = new Vector3(0.12f, 0.45f, 0.45f),
            emitVelocity = new Vector3(0.25f, 0f, 0f),
            ventCenter = new Vector3(4.75f, 0.5f, 2.35f),
            ventHalf = new Vector3(0.12f, 0.45f, 0.35f),
            ventSuction = Vector3.zero,
            particleDispersion = 0f,
            dispersionWallBoost = 0f,
            ventDrift = 0f,
            wallDeflect = 0f,
            pushSearchCells = 36,
            kinematicViscosity = 0.0004f,
            nParticles = 60,
            particleSize = 0.06f,
            stepsPerFrame = 2,
            useLES = true
        };
    }

    public void ApplyTo(SFUnifiedCS target)
    {
        ApplyTo(target, target.CurrentScenario);
    }

    public void ApplyTo(SFUnifiedCS target, SFUnifiedCS.ScenarioType scenario)
    {
        switch (scenario)
        {
            case SFUnifiedCS.ScenarioType.Jet:
                target.ApplyJetPreset(jet);
                break;
            case SFUnifiedCS.ScenarioType.SphereObstacle:
                target.ApplySphereObstaclePreset(sphereObstacle);
                break;
            case SFUnifiedCS.ScenarioType.CylinderObstacle:
                target.ApplyCylinderObstaclePreset(cylinderObstacle);
                break;
            case SFUnifiedCS.ScenarioType.TwoSpheres:
                target.ApplyTwoSpheresPreset(twoSpheres);
                break;
            case SFUnifiedCS.ScenarioType.LeapfrogRings:
                target.ApplyLeapfrogRingsPreset(leapfrogRings);
                break;
            case SFUnifiedCS.ScenarioType.Cigarette:
                target.ApplyCigaretteHipDefaults();
                break;
            case SFUnifiedCS.ScenarioType.InkCollision:
                target.ApplyInkCollisionHipDefaults();
                break;
            case SFUnifiedCS.ScenarioType.ObliqueRingCollision:
                target.ApplyObliqueRingCollisionHipDefaults();
                break;
            case SFUnifiedCS.ScenarioType.SmokeMaze2D:
                target.ApplySmokeMaze2DPreset(smokeMaze2D);
                break;
        }
    }
}
