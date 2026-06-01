using UnityEngine;
using ComputeShaderSF;
using ShrodingerFlow.Particles;
#if UNITY_EDITOR
using UnityEditor;
#endif

// «Apply Scenario Defaults» берёт числа из SFUnifiedScenarioPresets (ассет или встроенная копия),
// а не из switch в этом файле. Референс: JetExampleCS.unity, UnifiedCS.unity — см. SFUnifiedScenarioPresets.cs.

public class SFUnifiedCS : SFBase, ISimulationParticleSizeSource, IRaymarchDensitySource
{
    public enum ScenarioType
    {
        Jet,
        SphereObstacle,
        CylinderObstacle,
        TwoSpheres,
        LeapfrogRings,
        /// <summary>example_cigarette.hip: фон U, гравитация на ψ₂, heat по ψ₁ в сфере, граница jet ↑.</summary>
        Cigarette,
        /// <summary>example_ink_collision.hip / карточка Ink drop: два шара ±1, скорости ∓1 по X, границы ψ на каждом шаге.</summary>
        InkCollision,
        /// <summary>Вид сверху (XZ): кольцо слева, нормаль +X; кольцо с большой Z, нормаль −Z; плоскости перпендикулярны, встреча в центре объёма.</summary>
        ObliqueRingCollision,
        /// <summary>2D-лабиринт (вид сверху XZ): источник слева, перегородки с разными проходами, вытяжка справа сверху.</summary>
        SmokeMaze2D
    }

    [Header("Compute Shaders")]
    [SerializeField] private ComputeShader _kernelsShader;
    [SerializeField] private ComputeShader _fftShader;
    [SerializeField] private ComputeShader _particlesShader;
    [SerializeField] private ComputeShader _lesShader;

    [Header("Пресеты (контекстное меню Apply Scenario Defaults)")]
    [Tooltip("Если задано — используются эти значения. Иначе в Editor подставляется SFUnifiedScenarioPresets.asset, в рантайме/без ассета — встроенная копия (CreateBuiltIn).")]
    [SerializeField] private SFUnifiedScenarioPresets _scenarioPresetsOverride;

    [Header("Тип сценария")]
    [SerializeField] private ScenarioType _scenario = ScenarioType.TwoSpheres;

    [Header("Базовые параметры ISF")]
    [SerializeField] private int[] vol_size = { 4, 2, 2 };
    [SerializeField] private int[] vol_res = { 64, 32, 32 };
    [SerializeField] private float hbar = 0.1f;
    [SerializeField] private float dt = 1f / 12f;

    [Header("Скорость потока / фона")]
    [SerializeField] private Vector3 _velocity = new Vector3(-0.2f, 0f, 0f);

    [Header("Jet / Nozzle (Jet, SphereObstacle, CylinderObstacle)")]
    [SerializeField] private Vector3 _nozzleCen = new Vector3(0.3f, 0.9656632f, 1.0659939f);
    [SerializeField] private float _nozzleLen = 0.5f;
    [SerializeField] private float _nozzleRad = 0.3f;

    [Header("Obstacle 1 (SphereObstacle, CylinderObstacle, TwoSpheres)")]
    [SerializeField] private Vector3 _obstaclePos1 = new Vector3(1.5f, 1f, 1f);
    [SerializeField] private float _obstacleRadius1 = 0.5f;

    [Header("Obstacle 2 (TwoSpheres)")]
    [SerializeField] private Vector3 _obstaclePos2 = new Vector3(2.5f, 1f, 1f);
    [SerializeField] private float _obstacleRadius2 = 0.5f;

    [Header("Leapfrog Rings")]
    [SerializeField] private float _ring1Radius = 1.5f;
    [SerializeField] private float _ring2Radius = 0.9f;
    [SerializeField] private Vector3 _ring1Normal = new Vector3(-1, 0, 0);
    [SerializeField] private Vector3 _ring2Normal = new Vector3(-1, 0, 0);

    [Header("Oblique rings — вид сверху XZ: слева → +X, сверху (+Z) → −Z, центр объёма")]
    [Tooltip("AddCircle: ось вдоль нормали; перенос вихря — против +n, поэтому заданы −X и +Z, чтобы полёт был к центру (+X и −Z).")]
    [SerializeField] private Vector3 _obliqueRing1Center = new Vector3(1.4f, 2.5f, 2.5f);
    [SerializeField] private Vector3 _obliqueRing1Normal = new Vector3(-1f, 0f, 0f);
    [SerializeField] private Vector3 _obliqueRing2Center = new Vector3(2.5f, 2.5f, 3.6f);
    [SerializeField] private Vector3 _obliqueRing2Normal = new Vector3(0f, 0f, 1f);
    [Tooltip("Радиус трубки в AddCircle (расстояние от оси center+n), не «большой» радиус тора.")]
    [SerializeField] private float _obliqueRingRadius = 0.6f;
    [SerializeField] private float _obliquePsi2Re = 0.05f;

    [Header("Smoke Maze 2D — лабиринт (вид сверху XZ)")]
    [Tooltip("Толщина перегородок и внешних стен.")]
    [SerializeField] private float _mazeWallThickness = 0.12f;
    [Tooltip("Отступ внешних стен от границы домена.")]
    [SerializeField] private float _mazeWallMargin = 0.1f;
    [Tooltip("Перегородка 1: X; проём по Z между gapZ.x и gapZ.y (середина открыта).")]
    [SerializeField] private float _mazeWall1X = 1.0f;
    [SerializeField] private Vector2 _mazeWall1GapZ = new Vector2(1.0f, 2.0f);
    [Tooltip("Перегородка 2: X; сплошная по Z между solidZ.x и solidZ.y (верх/низ открыты).")]
    [SerializeField] private float _mazeWall2X = 2.5f;
    [SerializeField] private Vector2 _mazeWall2SolidZ = new Vector2(0.7f, 2.3f);
    [Tooltip("Перегородка 3: X; сплошная снизу до solidMaxZ (верх открыт).")]
    [SerializeField] private float _mazeWall3X = 3.8f;
    [SerializeField] private float _mazeWall3SolidMaxZ = 2.0f;
    [Tooltip("Источник дыма (оранжевый на схеме): центр и полуразмер AABB.")]
    [SerializeField] private Vector3 _mazeSourceCenter = new Vector3(0.25f, 0.5f, 1.5f);
    [SerializeField] private Vector3 _mazeSourceHalf = new Vector3(0.12f, 0.45f, 0.25f);
    [Tooltip("Скорость выброса из источника: умеренный поток вправо. Не 0 — иначе дым «застывает» у источника.")]
    [SerializeField] private Vector3 _mazeEmitVelocity = new Vector3(0.25f, 0f, 0f);
    [Tooltip("Вытяжка (зелёная): центр и полуразмер AABB.")]
    [SerializeField] private Vector3 _mazeVentCenter = new Vector3(4.75f, 0.5f, 2.35f);
    [SerializeField] private Vector3 _mazeVentHalf = new Vector3(0.12f, 0.45f, 0.35f);
    [Tooltip("Опциональный фазовый драйвер в зоне вытяжки. Обычно 0: тяга задается мягким drift, чтобы не создавать обратный поток к источнику.")]
    [SerializeField] private Vector3 _mazeVentSuction = Vector3.zero;
    [Tooltip("Турбулентная дисперсия трассеров в плоскости XZ (доли ячейки): рассеивание для поиска смещённых проходов.")]
    [SerializeField, Range(0f, 1f)] private float _mazeParticleDispersion = 0.18f;
    [Tooltip("Усиление дисперсии у стен (×) — помогает огибать препятствия.")]
    [SerializeField, Range(0f, 5f)] private float _mazeDispersionWallBoost = 1.2f;
    [Tooltip("Скорость мягкого продольного потока трассеров по лабиринту; у стен добавляется curl.")]
    [SerializeField] private float _mazeVentDrift = 0.12f;
    [Tooltip("Шаг обхода стены вверх/вниз к проёму (доли ячейки).")]
    [SerializeField, Range(0f, 2f)] private float _mazeWallDeflect = 0.55f;
    [Tooltip("Радиус поиска свободной ячейки при столкновении (в ячейках сетки; для проходов нужно ≥30).")]
    [SerializeField, Range(4, 48)] private int _mazePushSearchCells = 36;

    [Header("Cigarette — example_cigarette.hip")]
    [Tooltip("Фоновый поток U для начальной плоской волны (k = U/hbar).")]
    [SerializeField] private Vector3 _cigaretteBackgroundU = new Vector3(0.1f, 0f, 0f);
    [Tooltip("Jet для penalization на маске (в hip: 0, 1, 0).")]
    [SerializeField] private Vector3 _cigaretteJet = new Vector3(0f, 1f, 0f);
    [Tooltip("g: фаза на ψ₂, dot(g,P)*dt/hbar.")]
    [SerializeField] private Vector3 _cigaretteGravity = new Vector3(0f, 1f, 0f);
    [Tooltip("Центр сферы в координатах объёма [0..vol] (hip nozzle → Unity).")]
    [SerializeField] private Vector3 _cigaretteHeatSphereCen = new Vector3(1f, 0.5f, 1.6127148f);
    [SerializeField] private float _cigaretteHeatSphereRad = 0.2f;

    [Header("Частицы")]
    [SerializeField] private int _nParticles = 50;
    [SerializeField] private float _particleSize = 0.1f;

    [Header("Начальное расположение частиц (Box-спавн)")]
    [SerializeField] private Vector2 _boxSpawnX = new Vector2(0.3f, 0.3f);
    [SerializeField] private Vector2 _boxSpawnY = new Vector2(0.5f, 1.5f);
    [SerializeField] private Vector2 _boxSpawnZ = new Vector2(0.5f, 1.5f);

    [Header("ISF — гравитация на ψ₂ (сила через фазу g·x, не через v)")]
    [Tooltip("Для сценариев кроме Cigarette: тот же шаг GravityPsi2, что в example_cigarette (после первой нормировки).")]
    [SerializeField] private bool _applyPsi2Gravity;
    [Tooltip("Вектор g в фазе (g·P)·dt/ℏ на ψ₂; масштаб как у Cigarette — подбором.")]
    [SerializeField] private Vector3 _psi2Gravity = new Vector3(0f, 1f, 0f);

    [Header("Управление")]
    [SerializeField] private bool _useLES = false;
    [Tooltip("Кинематическая ν: лапласиан поля скорости после извлечения из ψ; при LES добавляется к ν_t Smagorinsky (единицы как dx²/шаг).")]
    [SerializeField] private float _kinematicViscosity;
    [SerializeField] private bool _paused;
    [SerializeField, Range(1, 20)] private int _stepsPerFrame = 3;

    private CSISF _isf;
    private CSParticles _particles;
    private CSVelocity _vel;
    private ComputeBuffer _maskBuf1, _maskBuf2;
    private ComputeBuffer _mazeWallMaskBuf, _mazeSourceMaskBuf, _mazeVentMaskBuf;
    private Vector3 _mazeVentMin, _mazeVentMax;

    private ParticleGpuBuffers _particleBuffers;
    private ParticleDisplay3D _particleDisplay;
    private Vector3[] _renderPos;
    private Vector3[] _renderVel;
    private float[] _pxArr, _pyArr, _pzArr;
    private Vector3[] _prevPos;
    /// <summary>Сглаженная скорость для цвета в шейдере (те же индексы, что у частиц GPU).</summary>
    private Vector3[] _displayVelSmooth;
    private int _particlesCount;

    private const float DisplayVelocityBlend = 0.32f;

    private float _kvecX, _kvecY, _kvecZ;
    private float _omega;
    private bool _initialized;
    private int _compactCounter;
    private bool _spawnEachStep;
    private bool _boundaryEachStep;

    /// <summary>Последнее значение <see cref="_particleSize"/>, с которым синхронизировали <see cref="ParticleDisplay3D.scale"/>.</summary>
    private float _particleSizeSyncedForDisplay = float.NaN;

    /// <summary>Текущий выбранный сценарий (для применения пресетов из SFUnifiedScenarioPresets).</summary>
    public ScenarioType CurrentScenario => _scenario;

    public float SimulationParticleSize => _particleSize;

    #region Lifecycle

    private void Start()
    {
        _particleBuffers = GetComponent<ParticleGpuBuffers>();
        _particleDisplay = GetComponent<ParticleDisplay3D>();

        _isf = new CSISF();
        _isf.Init(_kernelsShader, _fftShader, _lesShader, vol_size, vol_res, hbar, dt);

        bool oneTimeParticles = _scenario == ScenarioType.LeapfrogRings
                             || _scenario == ScenarioType.ObliqueRingCollision
                             || _scenario == ScenarioType.TwoSpheres
                             || _scenario == ScenarioType.InkCollision;
        int maxParticles = oneTimeParticles
            ? _nParticles
            : _nParticles * (_scenario == ScenarioType.SmokeMaze2D ? 4000 : 1000);

        _particles = new CSParticles();
        _particles.Init(_particlesShader, maxParticles, _isf);

        _vel = new CSVelocity(_isf.resX, _isf.resY, _isf.resZ);

        _particleBuffers?.EnsureCapacity(maxParticles);
        _renderPos = new Vector3[maxParticles];
        _renderVel = new Vector3[maxParticles];
        _pxArr = new float[maxParticles];
        _pyArr = new float[maxParticles];
        _pzArr = new float[maxParticles];
        _prevPos = new Vector3[maxParticles];
        _displayVelSmooth = new Vector3[maxParticles];

        InitScenario();
        _initialized = true;
        if (float.IsNaN(_particleSizeSyncedForDisplay))
            _particleSizeSyncedForDisplay = _particleSize;
        SyncParticleDisplayScaleFromSimulation();
    }

    private void OnValidate()
    {
        if (float.IsNaN(_particleSizeSyncedForDisplay))
            _particleSizeSyncedForDisplay = _particleSize;

        if (_particleDisplay == null)
            _particleDisplay = GetComponent<ParticleDisplay3D>();

        if (_particleDisplay != null && _particleDisplay.AutomaticSimulationScale)
            SyncParticleDisplayScaleFromSimulation();
        else if (!Mathf.Approximately(_particleSize, _particleSizeSyncedForDisplay))
            _particleSizeSyncedForDisplay = _particleSize;
    }

    private void SyncParticleDisplayScaleFromSimulation()
    {
        if (_particleDisplay == null)
            _particleDisplay = GetComponent<ParticleDisplay3D>();
        if (_particleDisplay == null || !_particleDisplay.AutomaticSimulationScale)
            return;
        _particleDisplay.ApplyAutomaticScaleFromSimulation(_particleSize);
        _particleSizeSyncedForDisplay = _particleSize;
    }

    private void Update()
    {
        if (!_initialized) return;

        if (!_paused)
        {
            for (int step = 0; step < _stepsPerFrame; step++)
            {
                iterator++;
                SimulationStep();
            }

            if (_spawnEachStep
                && (_scenario == ScenarioType.Jet || _scenario == ScenarioType.Cigarette
                    || _scenario == ScenarioType.SmokeMaze2D))
            {
                _compactCounter++;
                if (_compactCounter >= 60)
                {
                    _compactCounter = 0;
                    _particles.CompactParticles(_pxArr, _pyArr, _pzArr,
                        vol_size[0], vol_size[1], vol_size[2], _prevPos, _displayVelSmooth);
                    _particlesCount = _particles.Size;
                }
            }
        }

        UpdateParticleSystem();
    }

    private void OnDestroy()
    {
        _maskBuf1?.Release();
        _maskBuf2?.Release();
        _mazeWallMaskBuf?.Release();
        _mazeSourceMaskBuf?.Release();
        _mazeVentMaskBuf?.Release();
        _vel?.Dispose();
        _particles?.Dispose();
        _isf?.Dispose();
    }

    #endregion

    #region Scenario Init

    private void InitScenario()
    {
        switch (_scenario)
        {
            case ScenarioType.Jet:
                InitPsiUniform();
                _maskBuf1 = BuildJetMask();
                ComputeKvecAndOmega(_velocity);
                RunInitBoundary(_maskBuf1, _kvecX, _kvecY, _kvecZ, 0f, 10);
                _spawnEachStep = true;
                _boundaryEachStep = true;
                break;

            case ScenarioType.SphereObstacle:
                InitPsiWithPhase();
                _maskBuf1 = BuildSphereMask(_obstaclePos1, _obstacleRadius1);
                RunInitBoundary(_maskBuf1, 0f, 0f, 0f, 0f, 10);
                _kvecX = _kvecY = _kvecZ = 0f;
                _omega = 0f;
                _spawnEachStep = true;
                _boundaryEachStep = true;
                break;

            case ScenarioType.CylinderObstacle:
                InitPsiWithPhase();
                _maskBuf1 = BuildCylinderMask(_obstaclePos1, _obstacleRadius1);
                RunInitBoundary(_maskBuf1, 0f, 0f, 0f, 0f, 10);
                _kvecX = _kvecY = _kvecZ = 0f;
                _omega = 0f;
                _spawnEachStep = true;
                _boundaryEachStep = true;
                break;

            case ScenarioType.TwoSpheres:
                InitPsiUniform();
                _maskBuf1 = BuildSphereMask(_obstaclePos1, _obstacleRadius1);
                _maskBuf2 = BuildSphereMask(_obstaclePos2, _obstacleRadius2);
                float kx = _velocity.x / hbar;
                float ky = _velocity.y / hbar;
                float kz = _velocity.z / hbar;
                RunInitTwoSpheres(kx, ky, kz, 10);
                SpawnParticlesInSpheres();
                _spawnEachStep = false;
                _boundaryEachStep = false;
                break;

            case ScenarioType.LeapfrogRings:
                InitPsiWithPhaseAndCircles();
                _isf.Normalize();
                _isf.PressureProject();
                SpawnParticlesInBox(_nParticles);
                _particlesCount = _particles.Size;
                _spawnEachStep = false;
                _boundaryEachStep = false;
                break;

            case ScenarioType.ObliqueRingCollision:
                InitPsiObliqueRingsHip();
                _isf.Normalize();
                _isf.PressureProject();
                SpawnParticlesOnObliqueRings(_nParticles);
                _particlesCount = _particles.Size;
                _spawnEachStep = false;
                _boundaryEachStep = false;
                break;

            case ScenarioType.Cigarette:
                InitPsiPlaneWave(_cigaretteBackgroundU);
                _maskBuf1 = BuildSphereMask(_cigaretteHeatSphereCen, _cigaretteHeatSphereRad);
                ComputeKvecAndOmega(_cigaretteJet);
                RunInitBoundary(_maskBuf1, _kvecX, _kvecY, _kvecZ, 0f, 10);
                _spawnEachStep = true;
                _boundaryEachStep = true;
                break;

            case ScenarioType.InkCollision:
                InitPsiUniform();
                _maskBuf1 = BuildSphereMask(_obstaclePos1, _obstacleRadius1);
                _maskBuf2 = BuildSphereMask(_obstaclePos2, _obstacleRadius2);
                float ikx = _velocity.x / hbar;
                float iky = _velocity.y / hbar;
                float ikz = _velocity.z / hbar;
                RunInitTwoSpheres(ikx, iky, ikz, 10);
                ComputeKvecAndOmega(_velocity);
                SpawnParticlesInSpheres();
                _spawnEachStep = false;
                _boundaryEachStep = true;
                break;

            case ScenarioType.SmokeMaze2D:
                InitPsiPlaneWave(new Vector3(0.06f, 0f, 0f));
                _mazeWallMaskBuf = BuildMazeWallMask();
                _mazeSourceMaskBuf = BuildBoxMask(_mazeSourceCenter, _mazeSourceHalf);
                _mazeVentMaskBuf = BuildBoxMask(_mazeVentCenter, _mazeVentHalf);
                _mazeVentMin = _mazeVentCenter - _mazeVentHalf;
                _mazeVentMax = _mazeVentCenter + _mazeVentHalf;
                RunInitBoundary(_mazeWallMaskBuf, 0f, 0f, 0f, 0f, 10);
                _spawnEachStep = true;
                _boundaryEachStep = false;
                break;
        }
    }

    #endregion

    #region Psi Initialization

    private void InitPsiUniform()
    {
        int num = _isf.num;
        var tmp1 = new Vector2[num];
        var tmp2 = new Vector2[num];
        for (int i = 0; i < num; i++)
        {
            tmp1[i] = new Vector2(1f, 0f);
            tmp2[i] = new Vector2(0.01f, 0f);
        }
        _isf.psi1.SetData(tmp1);
        _isf.psi2.SetData(tmp2);
        _isf.Normalize();
    }

    private void InitPsiWithPhase()
    {
        InitPsiPlaneWave(_velocity);
    }

    /// <summary>Set_background_flow из example_cigarette.hip: ψ₁=exp(i k·P), ψ₂=0.01·exp(i k·P), k=U/hbar.</summary>
    private void InitPsiPlaneWave(Vector3 backgroundU)
    {
        int num = _isf.num;
        float kx = backgroundU.x / hbar;
        float ky = backgroundU.y / hbar;
        float kz = backgroundU.z / hbar;

        var tmp1 = new Vector2[num];
        var tmp2 = new Vector2[num];
        for (int i = 0; i < num; i++)
        {
            float phase = kx * _isf.pxCPU[i] + ky * _isf.pyCPU[i] + kz * _isf.pzCPU[i];
            float c = Mathf.Cos(phase), s = Mathf.Sin(phase);
            tmp1[i] = new Vector2(c, s);
            tmp2[i] = new Vector2(c * 0.01f, s * 0.01f);
        }
        _isf.psi1.SetData(tmp1);
        _isf.psi2.SetData(tmp2);
        _isf.Normalize();
    }

    private void InitPsiWithPhaseAndCircles()
    {
        int num = _isf.num;
        float kx = _velocity.x / hbar;
        float ky = _velocity.y / hbar;
        float kz = _velocity.z / hbar;

        var tmp1 = new Vector2[num];
        var tmp2 = new Vector2[num];
        for (int i = 0; i < num; i++)
        {
            float phase = kx * _isf.pxCPU[i] + ky * _isf.pyCPU[i] + kz * _isf.pzCPU[i];
            float c = Mathf.Cos(phase), s = Mathf.Sin(phase);
            tmp1[i] = new Vector2(c, s);
            tmp2[i] = new Vector2(c * 0.01f, s * 0.01f);
        }

        float d = _isf.dx * 5f;
        Vector3 center = new Vector3(vol_size[0] / 2f, vol_size[1] / 2f, vol_size[2] / 2f);

        AddCircle(tmp1, center, _ring1Normal, _ring1Radius, d);
        AddCircle(tmp1, center, _ring2Normal, _ring2Radius, d);

        _isf.psi1.SetData(tmp1);
        _isf.psi2.SetData(tmp2);
    }

    /// <summary>Set_Constant_One + Add_Vortex_Ring1/2 + Add_small_second_component из example_oblique_collision.hip (без фоновой exp(ik·P)).</summary>
    private void InitPsiObliqueRingsHip()
    {
        int num = _isf.num;
        var tmp1 = new Vector2[num];
        var tmp2 = new Vector2[num];
        for (int i = 0; i < num; i++)
        {
            tmp1[i] = new Vector2(1f, 0f);
            tmp2[i] = new Vector2(_obliquePsi2Re, 0f);
        }

        float d = _isf.dx * 5f;
        AddCircle(tmp1, _obliqueRing1Center, _obliqueRing1Normal, _obliqueRingRadius, d);
        AddCircle(tmp1, _obliqueRing2Center, _obliqueRing2Normal, _obliqueRingRadius, d);

        _isf.psi1.SetData(tmp1);
        _isf.psi2.SetData(tmp2);
    }

    private void AddCircle(Vector2[] psi, Vector3 center, Vector3 normal, float radius, float d)
    {
        normal = normal.normalized;
        int rx = _isf.resX, ry = _isf.resY, rz = _isf.resZ;

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
                    if (rPerp2 < radius * radius)
                    {
                        if (z > 0f && z <= d / 2f)
                            alpha = -Mathf.PI * (2f * z / d - 1f);
                        else if (z <= 0f && z >= -d / 2f)
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

    #endregion

    #region Mask Building

    private ComputeBuffer BuildJetMask()
    {
        int num = _isf.num;
        var mask = new int[num];
        for (int i = 0; i < num; i++)
        {
            float px = _isf.pxCPU[i];
            float py = _isf.pyCPU[i];
            float pz = _isf.pzCPU[i];
            bool inJet = Mathf.Abs(px - _nozzleCen.x) <= _nozzleLen / 2f
                && (py - _nozzleCen.y) * (py - _nozzleCen.y)
                 + (pz - _nozzleCen.z) * (pz - _nozzleCen.z) <= _nozzleRad * _nozzleRad;
            mask[i] = inJet ? 1 : 0;
        }
        var buf = new ComputeBuffer(num, sizeof(int));
        buf.SetData(mask);
        return buf;
    }

    private ComputeBuffer BuildSphereMask(Vector3 pos, float radius)
    {
        int num = _isf.num;
        var mask = new int[num];
        float r2 = radius * radius;
        for (int i = 0; i < num; i++)
        {
            float dx = _isf.pxCPU[i] - pos.x;
            float dy = _isf.pyCPU[i] - pos.y;
            float dz = _isf.pzCPU[i] - pos.z;
            mask[i] = (dx * dx + dy * dy + dz * dz <= r2) ? 1 : 0;
        }
        var buf = new ComputeBuffer(num, sizeof(int));
        buf.SetData(mask);
        return buf;
    }

    private ComputeBuffer BuildCylinderMask(Vector3 pos, float radius)
    {
        int num = _isf.num;
        var mask = new int[num];
        float r2 = radius * radius;
        for (int i = 0; i < num; i++)
        {
            float dx = _isf.pxCPU[i] - pos.x;
            float dy = _isf.pyCPU[i] - pos.y;
            mask[i] = (dx * dx + dy * dy <= r2) ? 1 : 0;
        }
        var buf = new ComputeBuffer(num, sizeof(int));
        buf.SetData(mask);
        return buf;
    }

    private static bool IsInsideAabb(float px, float py, float pz, Vector3 center, Vector3 half)
    {
        return px >= center.x - half.x && px <= center.x + half.x
            && py >= center.y - half.y && py <= center.y + half.y
            && pz >= center.z - half.z && pz <= center.z + half.z;
    }

    private ComputeBuffer BuildBoxMask(Vector3 center, Vector3 half)
    {
        int num = _isf.num;
        var mask = new int[num];
        for (int i = 0; i < num; i++)
            mask[i] = IsInsideAabb(_isf.pxCPU[i], _isf.pyCPU[i], _isf.pzCPU[i], center, half) ? 1 : 0;
        var buf = new ComputeBuffer(num, sizeof(int));
        buf.SetData(mask);
        return buf;
    }

    private bool IsMazeOuterWall(float px, float pz)
    {
        float m = _mazeWallMargin;
        float wx = vol_size[0], wz = vol_size[2];
        return px < m || px > wx - m || pz < m || pz > wz - m;
    }

    private bool IsMazeSlatSolid(float px, float pz)
    {
        float t = _mazeWallThickness * 0.5f;

        if (Mathf.Abs(px - _mazeWall1X) <= t)
            return pz < _mazeWall1GapZ.x || pz > _mazeWall1GapZ.y;

        if (Mathf.Abs(px - _mazeWall2X) <= t)
            return pz >= _mazeWall2SolidZ.x && pz <= _mazeWall2SolidZ.y;

        if (Mathf.Abs(px - _mazeWall3X) <= t)
            return pz < _mazeWall3SolidMaxZ;

        return false;
    }

    /// <summary>Внешние стены + 3 перегородки с разными проходами; без источника и вытяжки.</summary>
    private ComputeBuffer BuildMazeWallMask()
    {
        int num = _isf.num;
        var mask = new int[num];
        float hx = _isf.dx * 0.5f;
        float hz = _isf.dz * 0.5f;
        for (int i = 0; i < num; i++)
        {
            float px = _isf.pxCPU[i], py = _isf.pyCPU[i], pz = _isf.pzCPU[i];
            if (IsInsideAabb(px, py, pz, _mazeSourceCenter, _mazeSourceHalf)
                || IsInsideAabb(px, py, pz, _mazeVentCenter, _mazeVentHalf))
            {
                mask[i] = 0;
                continue;
            }
            // Углы ячейки по XZ — чтобы тонкие стены не «пропускали» поток между центрами соседей.
            bool solid = IsMazeSolidAt(px, py, pz)
                || IsMazeSolidAt(px - hx, py, pz - hz)
                || IsMazeSolidAt(px + hx, py, pz - hz)
                || IsMazeSolidAt(px - hx, py, pz + hz)
                || IsMazeSolidAt(px + hx, py, pz + hz);
            mask[i] = solid ? 1 : 0;
        }
        var buf = new ComputeBuffer(num, sizeof(int));
        buf.SetData(mask);
        return buf;
    }

    private bool IsMazeSolidAt(float px, float py, float pz)
    {
        if (IsInsideAabb(px, py, pz, _mazeSourceCenter, _mazeSourceHalf)) return false;
        if (IsInsideAabb(px, py, pz, _mazeVentCenter, _mazeVentHalf)) return false;
        return IsMazeOuterWall(px, pz) || IsMazeSlatSolid(px, pz);
    }

    private void ComputeKvecAndOmega(Vector3 vel)
    {
        _kvecX = vel.x / hbar;
        _kvecY = vel.y / hbar;
        _kvecZ = vel.z / hbar;
        _omega = vel.sqrMagnitude / (2f * hbar);
    }

    #endregion

    #region Boundary Initialization

    private void RunInitBoundary(ComputeBuffer mask, float kx, float ky, float kz,
        float phaseOffset, int iterations)
    {
        for (int iter = 0; iter < iterations; iter++)
        {
            _isf.ApplyJetBoundary(mask, kx, ky, kz, phaseOffset);
            _isf.PressureProject();
        }
    }

    private void RunInitTwoSpheres(float kx, float ky, float kz, int iterations)
    {
        for (int iter = 0; iter < iterations; iter++)
        {
            _isf.ApplyJetBoundary(_maskBuf1, kx, ky, kz, 0f);
            _isf.ApplyJetBoundary(_maskBuf2, -kx, -ky, -kz, 0f);
            _isf.PressureProject();
        }
    }

    #endregion

    #region Simulation

    private void SimulationStep()
    {
        if (_scenario == ScenarioType.SmokeMaze2D)
        {
            SimulationStepSmokeMaze();
            return;
        }

        _isf.kinematicViscosity = _kinematicViscosity;
        if (_scenario == ScenarioType.Cigarette)
            _isf.UpdateCigaretteSpace(_useLES, _cigaretteGravity, _maskBuf1);
        else
            _isf.UpdateSpace(_useLES,
                _applyPsi2Gravity ? _psi2Gravity : (Vector3?)null);

        if (_boundaryEachStep)
        {
            if (_scenario == ScenarioType.InkCollision)
            {
                float phaseOffset = -_omega * dt * iterator;
                _isf.ApplyJetBoundary(_maskBuf1, _kvecX, _kvecY, _kvecZ, phaseOffset);
                _isf.ApplyJetBoundary(_maskBuf2, -_kvecX, -_kvecY, -_kvecZ, phaseOffset);
                _isf.PressureProject();
            }
            else
            {
                float phaseOffset = (_scenario == ScenarioType.Jet
                                    || _scenario == ScenarioType.Cigarette)
                    ? -_omega * dt * iterator
                    : 0f;
                _isf.ApplyJetBoundary(_maskBuf1, _kvecX, _kvecY, _kvecZ, phaseOffset);
                _isf.PressureProject();
            }
        }

        if (_spawnEachStep)
            SpawnNozzleParticles();

        _isf.UpdateVelocities(_vel);
        _particles.CalculateMovement(_vel);

        if (_scenario != ScenarioType.Jet && _scenario != ScenarioType.Cigarette
            && _scenario != ScenarioType.ObliqueRingCollision)
            _particles.WrapPositions(vol_size[0], vol_size[1], vol_size[2]);
    }

    private void SimulationStepSmokeMaze()
    {
        _isf.kinematicViscosity = _kinematicViscosity;
        _isf.UpdateSpace(_useLES, null);

        float invH = 1f / hbar;
        _isf.ApplyJetBoundary(_mazeWallMaskBuf, 0f, 0f, 0f, 0f);
        _isf.PressureProject();
        _isf.ApplyJetBoundary(_mazeWallMaskBuf, 0f, 0f, 0f, 0f);
        _isf.ApplyJetBoundary(_mazeSourceMaskBuf,
            _mazeEmitVelocity.x * invH, _mazeEmitVelocity.y * invH, _mazeEmitVelocity.z * invH, 0f);
        if (_mazeVentSuction.sqrMagnitude > 1e-8f)
        {
            _isf.ApplyJetBoundary(_mazeVentMaskBuf,
                _mazeVentSuction.x * invH, _mazeVentSuction.y * invH, _mazeVentSuction.z * invH, 0f);
        }
        _isf.PressureProject();

        SpawnMazeSourceParticles();
        _isf.UpdateVelocities(_vel);
        _particles.CalculateMovement(_vel, clampSampling: true);
        _particles.DriftInSourceBox(_mazeSourceCenter, _mazeSourceHalf, _mazeEmitVelocity, dt);
        _particles.DriftTowardVent(_mazeVentCenter, _mazeVentDrift, dt, _mazeWallMaskBuf);
        float mcell = Mathf.Min(_isf.dx, _isf.dz);
        if (_mazeParticleDispersion > 1e-4f)
        {
            _particles.DisperseMaze(_mazeWallMaskBuf,
                _mazeParticleDispersion * mcell, _mazeDispersionWallBoost, iterator);
        }
        if (_mazeWallDeflect > 1e-4f)
        {
            _particles.DeflectAtWall(_mazeWallMaskBuf,
                _mazeWallDeflect * mcell, _mazePushSearchCells);
        }
        for (int pass = 0; pass < 2; pass++)
        {
            _particles.PushOutOfSolidMask(_mazeWallMaskBuf, _mazePushSearchCells,
                _mazeVentCenter, 0.45f);
        }
        _particles.ClampPositionsToVolume(vol_size[0], vol_size[1], vol_size[2]);
        CompactMazeParticles();
    }

    private void SpawnMazeSourceParticles()
    {
        var xx = new float[_nParticles];
        var yy = new float[_nParticles];
        var zz = new float[_nParticles];
        for (int i = 0; i < _nParticles; i++)
        {
            xx[i] = Random.Range(_mazeSourceCenter.x - _mazeSourceHalf.x, _mazeSourceCenter.x + _mazeSourceHalf.x);
            yy[i] = Random.Range(_mazeSourceCenter.y - _mazeSourceHalf.y, _mazeSourceCenter.y + _mazeSourceHalf.y);
            zz[i] = Random.Range(_mazeSourceCenter.z - _mazeSourceHalf.z, _mazeSourceCenter.z + _mazeSourceHalf.z);
        }
        _particles.AddParticles(xx, yy, zz, _nParticles);
        _particlesCount = _particles.Size;
    }

    /// <summary>Удалить частицы вне домена и в зоне вытяжки (имитация удаления дыма).</summary>
    private void CompactMazeParticles()
    {
        if (_particlesCount == 0) return;

        _particles.ReadPositions(_pxArr, _pyArr, _pzArr);
        for (int i = 0; i < _particlesCount; i++)
        {
            float px = _pxArr[i], py = _pyArr[i], pz = _pzArr[i];
            if (px >= _mazeVentMin.x && px <= _mazeVentMax.x
                && py >= _mazeVentMin.y && py <= _mazeVentMax.y
                && pz >= _mazeVentMin.z && pz <= _mazeVentMax.z)
            {
                _pxArr[i] = _pyArr[i] = _pzArr[i] = -1f;
            }
        }
        _particles.WritePositions(_pxArr, _pyArr, _pzArr, _particlesCount);
        _particles.CompactParticles(_pxArr, _pyArr, _pzArr,
            vol_size[0], vol_size[1], vol_size[2], _prevPos, _displayVelSmooth);
        _particlesCount = _particles.Size;
    }

    #endregion

    #region Particle Spawning

    private void SpawnNozzleParticles()
    {
        var xx = new float[_nParticles];
        var yy = new float[_nParticles];
        var zz = new float[_nParticles];

        if (_scenario == ScenarioType.CylinderObstacle)
        {
            for (int i = 0; i < _nParticles; i++)
            {
                xx[i] = _nozzleCen.x;
                yy[i] = Random.Range(_boxSpawnY.x, _boxSpawnY.y);
                zz[i] = Random.Range(_boxSpawnZ.x, _boxSpawnZ.y);
            }
        }
        else if (_scenario == ScenarioType.Cigarette)
        {
            float r = 0.9f * _cigaretteHeatSphereRad;
            for (int i = 0; i < _nParticles; i++)
            {
                Vector3 d = Random.onUnitSphere;
                xx[i] = _cigaretteHeatSphereCen.x + r * d.x;
                yy[i] = _cigaretteHeatSphereCen.y + r * d.y;
                zz[i] = _cigaretteHeatSphereCen.z + r * d.z;
            }
        }
        else
        {
            for (int i = 0; i < _nParticles; i++)
            {
                float rt = Random.value * 2f * Mathf.PI;
                xx[i] = _nozzleCen.x;
                yy[i] = _nozzleCen.y + 0.9f * _nozzleRad * Mathf.Cos(rt);
                zz[i] = _nozzleCen.z + 0.9f * _nozzleRad * Mathf.Sin(rt);
            }
        }

        bool ring = _scenario != ScenarioType.Jet && _scenario != ScenarioType.Cigarette;
        _particles.AddParticles(xx, yy, zz, _nParticles, ring);
        _particlesCount = _particles.Size;
    }

    private void SpawnParticlesInBox(int count)
    {
        var xx = new float[count];
        var yy = new float[count];
        var zz = new float[count];
        for (int i = 0; i < count; i++)
        {
            xx[i] = Random.Range(_boxSpawnX.x, _boxSpawnX.y);
            yy[i] = Random.Range(_boxSpawnY.x, _boxSpawnY.y);
            zz[i] = Random.Range(_boxSpawnZ.x, _boxSpawnZ.y);
        }
        _particles.AddParticles(xx, yy, zz, count);
        _particlesCount = _particles.Size;
    }

    /// <summary>Частицы у двух вихрей как у <see cref="AddCircle"/>: цилиндрическая трубка вдоль нормали (радиус = tube, длина ≈ <see cref="InitPsiObliqueRingsHip"/> d = 5·dx), а не окружность в одной плоскости.</summary>
    private void SpawnParticlesOnObliqueRings(int count)
    {
        if (count <= 0) return;
        var xx = new float[count];
        var yy = new float[count];
        var zz = new float[count];
        int half = count / 2;
        float d = _isf.dx * 5f;
        float halfAxial = d * 0.48f;
        float tubeR = _obliqueRingRadius;
        float jitter = _isf.dx * 0.6f;
        FillAddCircleVortexTube(xx, yy, zz, 0, half,
            _obliqueRing1Center, _obliqueRing1Normal, tubeR, halfAxial, jitter);
        FillAddCircleVortexTube(xx, yy, zz, half, count,
            _obliqueRing2Center, _obliqueRing2Normal, tubeR, halfAxial, jitter);
        _particles.AddParticles(xx, yy, zz, count);
        _particlesCount = _particles.Size;
    }

    /// <summary>Соответствует ядру AddCircle: ось center + s·n, |s|≤halfAxial; сечение — круг радиуса tubeRadius в плоскости ⊥ n.</summary>
    private static void FillAddCircleVortexTube(float[] xx, float[] yy, float[] zz,
        int from, int to, Vector3 center, Vector3 normal, float tubeRadius, float halfAxial, float jitter)
    {
        Vector3 n = normal.normalized;
        Vector3 aux = Mathf.Abs(n.y) < 0.99f ? Vector3.up : Vector3.right;
        Vector3 e1 = Vector3.Cross(aux, n);
        if (e1.sqrMagnitude < 1e-8f)
            e1 = Vector3.Cross(Vector3.forward, n);
        e1.Normalize();
        Vector3 e2 = Vector3.Cross(n, e1);
        for (int i = from; i < to; i++)
        {
            float s = Random.Range(-halfAxial, halfAxial);
            float t = Random.Range(0f, 2f * Mathf.PI);
            Vector3 p = center + s * n + tubeRadius * (Mathf.Cos(t) * e1 + Mathf.Sin(t) * e2);
            if (jitter > 0f)
                p += Random.insideUnitSphere * jitter;
            xx[i] = p.x;
            yy[i] = p.y;
            zz[i] = p.z;
        }
    }

    private void SpawnParticlesInSpheres()
    {
        int half = _nParticles / 2;
        var xx = new float[_nParticles];
        var yy = new float[_nParticles];
        var zz = new float[_nParticles];

        FillSphereParticles(xx, yy, zz, 0, half,
            _obstaclePos2, new Vector3(_obstacleRadius2, _obstacleRadius2, _obstacleRadius2));
        FillSphereParticles(xx, yy, zz, half, _nParticles,
            _obstaclePos1, new Vector3(_obstacleRadius1, _obstacleRadius1, _obstacleRadius1));

        _particles.AddParticles(xx, yy, zz, _nParticles);
        _particlesCount = _particles.Size;
    }

    private static void FillSphereParticles(float[] xx, float[] yy, float[] zz,
        int from, int to, Vector3 center, Vector3 size)
    {
        float r2 = size.x * size.x;
        for (int i = from; i < to; i++)
        {
            float px, py, pz;
            do
            {
                px = Random.Range(center.x - size.x, center.x + size.x);
                py = Random.Range(center.y - size.y, center.y + size.y);
                pz = Random.Range(center.z - size.z, center.z + size.z);
            } while ((px - center.x) * (px - center.x)
                   + (py - center.y) * (py - center.y)
                   + (pz - center.z) * (pz - center.z) > r2);

            xx[i] = px;
            yy[i] = py;
            zz[i] = pz;
        }
    }

    #endregion

    #region Particle Rendering

    private void UpdateParticleSystem()
    {
        if (_particlesCount == 0 || _particleBuffers == null) return;

        _particles.ReadPositions(_pxArr, _pyArr, _pzArr);

        var offset = transform.position;

        bool cull = _scenario == ScenarioType.Jet || _scenario == ScenarioType.Cigarette
            || _scenario == ScenarioType.ObliqueRingCollision
            || _scenario == ScenarioType.SmokeMaze2D;
        float maxX = vol_size[0], maxY = vol_size[1], maxZ = vol_size[2];
        float velThreshold = maxX * maxX + maxY * maxY + maxZ * maxZ;
        int visible = 0;

        for (int i = 0; i < _particlesCount; i++)
        {
            float px = _pxArr[i], py = _pyArr[i], pz = _pzArr[i];

            var pos = new Vector3(px, py, pz) + offset;
            var lastPos = _prevPos[i];
            _prevPos[i] = pos;

            var vel = pos - lastPos;
            if (vel.sqrMagnitude > velThreshold)
                vel = Vector3.zero;

            _displayVelSmooth[i] = Vector3.Lerp(_displayVelSmooth[i], vel, Mathf.Clamp01(DisplayVelocityBlend));

            if (cull && (px < 0f || px > maxX || py < 0f || py > maxY || pz < 0f || pz > maxZ))
                continue;

            _renderPos[visible] = pos;
            _renderVel[visible] = _displayVelSmooth[i];
            visible++;
        }

        _particleBuffers.Upload(_renderPos, _renderVel, visible);
    }

    #endregion

    #region Presets

    private SFUnifiedScenarioPresets ResolveScenarioPresets()
    {
        if (_scenarioPresetsOverride != null)
            return _scenarioPresetsOverride;
#if UNITY_EDITOR
        var fromProject = AssetDatabase.LoadAssetAtPath<SFUnifiedScenarioPresets>(
            SFUnifiedScenarioPresets.DefaultAssetPath);
        if (fromProject != null)
            return fromProject;
#endif
        return SFUnifiedScenarioPresets.CreateBuiltIn();
    }

    [ContextMenu("Apply Scenario Defaults")]
    private void ApplyScenarioDefaults()
    {
        var src = ResolveScenarioPresets();
        src.ApplyTo(this);
        Debug.Log($"[SFUnifiedCS] Applied scenario defaults for {_scenario} (presets: {(src == _scenarioPresetsOverride ? "override field" : "asset / built-in")})");
    }

    public void ApplyJetPreset(SFUnifiedJetPreset p)
    {
        vol_size = (int[])p.vol_size?.Clone() ?? new[] { 4, 2, 2 };
        vol_res = (int[])p.vol_res?.Clone() ?? new[] { 128, 64, 64 };
        hbar = p.hbar;
        dt = p.dt;
        _velocity = p.velocity;
        _nozzleCen = p.nozzleCen;
        _nozzleLen = p.nozzleLen;
        _nozzleRad = p.nozzleRad;
        _nParticles = p.nParticles;
        _particleSize = p.particleSize;
        SyncParticleDisplayScaleFromSimulation();
        _stepsPerFrame = p.stepsPerFrame;
        _useLES = p.useLES;
    }

    public void ApplySphereObstaclePreset(SFUnifiedSphereObstaclePreset p)
    {
        vol_size = (int[])p.vol_size?.Clone() ?? new[] { 4, 2, 2 };
        vol_res = (int[])p.vol_res?.Clone() ?? new[] { 64, 32, 32 };
        hbar = p.hbar;
        dt = p.dt;
        _velocity = p.velocity;
        _obstaclePos1 = p.obstaclePos1;
        _obstacleRadius1 = p.obstacleRadius1;
        _nozzleCen = p.nozzleCen;
        _nozzleLen = p.nozzleLen;
        _nozzleRad = p.nozzleRad;
        _boxSpawnX = p.boxSpawnX;
        _boxSpawnY = p.boxSpawnY;
        _boxSpawnZ = p.boxSpawnZ;
        _nParticles = p.nParticles;
        _particleSize = p.particleSize;
        SyncParticleDisplayScaleFromSimulation();
        _stepsPerFrame = p.stepsPerFrame;
        _useLES = p.useLES;
    }

    public void ApplyCylinderObstaclePreset(SFUnifiedCylinderObstaclePreset p)
    {
        vol_size = (int[])p.vol_size?.Clone() ?? new[] { 4, 2, 2 };
        vol_res = (int[])p.vol_res?.Clone() ?? new[] { 64, 32, 32 };
        hbar = p.hbar;
        dt = p.dt;
        _velocity = p.velocity;
        _obstaclePos1 = p.obstaclePos1;
        _obstacleRadius1 = p.obstacleRadius1;
        _nozzleCen = p.nozzleCen;
        _boxSpawnX = p.boxSpawnX;
        _boxSpawnY = p.boxSpawnY;
        _boxSpawnZ = p.boxSpawnZ;
        _nParticles = p.nParticles;
        _particleSize = p.particleSize;
        SyncParticleDisplayScaleFromSimulation();
        _stepsPerFrame = p.stepsPerFrame;
        _useLES = p.useLES;
    }

    public void ApplyTwoSpheresPreset(SFUnifiedTwoSpheresPreset p)
    {
        vol_size = (int[])p.vol_size?.Clone() ?? new[] { 4, 2, 2 };
        vol_res = (int[])p.vol_res?.Clone() ?? new[] { 64, 32, 32 };
        hbar = p.hbar;
        dt = p.dt;
        _velocity = p.velocity;
        _obstaclePos1 = p.obstaclePos1;
        _obstacleRadius1 = p.obstacleRadius1;
        _obstaclePos2 = p.obstaclePos2;
        _obstacleRadius2 = p.obstacleRadius2;
        _boxSpawnX = p.boxSpawnX;
        _boxSpawnY = p.boxSpawnY;
        _boxSpawnZ = p.boxSpawnZ;
        _nParticles = p.nParticles;
        _particleSize = p.particleSize;
        SyncParticleDisplayScaleFromSimulation();
        _stepsPerFrame = p.stepsPerFrame;
        _useLES = p.useLES;
    }

    public void ApplyLeapfrogRingsPreset(SFUnifiedLeapfrogRingsPreset p)
    {
        vol_size = (int[])p.vol_size?.Clone() ?? new[] { 10, 5, 5 };
        vol_res = (int[])p.vol_res?.Clone() ?? new[] { 128, 64, 64 };
        hbar = p.hbar;
        dt = p.dt;
        _velocity = p.velocity;
        _ring1Radius = p.ring1Radius;
        _ring2Radius = p.ring2Radius;
        _ring1Normal = p.ring1Normal;
        _ring2Normal = p.ring2Normal;
        _boxSpawnX = p.boxSpawnX;
        _boxSpawnY = p.boxSpawnY;
        _boxSpawnZ = p.boxSpawnZ;
        _nParticles = p.nParticles;
        _particleSize = p.particleSize;
        SyncParticleDisplayScaleFromSimulation();
        _stepsPerFrame = p.stepsPerFrame;
        _useLES = p.useLES;
    }

    /// <summary>Параметры как в example_cigarette.hip. Задайте до входа в Play (инициализация CSISF в Start).</summary>
    [ContextMenu("Apply Cigarette (hip) defaults")]
    public void ApplyCigaretteHipDefaults()
    {
        _scenario = ScenarioType.Cigarette;
        vol_size = new[] { 3, 6, 3 };
        vol_res = new[] { 64, 128, 64 };
        hbar = 0.03f;
        dt = 1f / 48f;
        _useLES = false;
        _stepsPerFrame = 1;
        _cigaretteBackgroundU = new Vector3(0.1f, 0f, 0f);
        _cigaretteJet = new Vector3(0f, 1f, 0f);
        _cigaretteGravity = new Vector3(0f, 1f, 0f);
        _cigaretteHeatSphereCen = new Vector3(1f, 0.5f, 1.6127148f);
        _cigaretteHeatSphereRad = 0.2f;
        _nParticles = 50;
        _particleSize = 0.05f;
        SyncParticleDisplayScaleFromSimulation();
    }

    /// <summary>Карточка Ink drop + example_ink_collision.hip (домен 0…4, центры сфер в hip −1 и +1 по X → Unity 1 и 3 при Y,Z=2).</summary>
    [ContextMenu("Apply Ink collision (hip card) defaults")]
    public void ApplyInkCollisionHipDefaults()
    {
        _scenario = ScenarioType.InkCollision;
        vol_size = new[] { 4, 4, 4 };
        vol_res = new[] { 128, 128, 128 };
        hbar = 0.02f;
        dt = 1f / 48f;
        _velocity = new Vector3(1f, 0f, 0f);
        _obstaclePos1 = new Vector3(1f, 2f, 2f);
        _obstaclePos2 = new Vector3(3f, 2f, 2f);
        _obstacleRadius1 = 0.45f;
        _obstacleRadius2 = 0.45f;
        _useLES = false;
        _stepsPerFrame = 1;
        _nParticles = 100000;
        _particleSize = 0.1f;
        SyncParticleDisplayScaleFromSimulation();
    }

    /// <summary>Карточка 5³, hbar 0.05, 64³, dt 1/24, r=0.6. Схема сверху: слева ось +X, сверху по Z ось −Z, центры C−approach·n.</summary>
    [ContextMenu("Apply Oblique ring collision (hip card) defaults")]
    public void ApplyObliqueRingCollisionHipDefaults()
    {
        _scenario = ScenarioType.ObliqueRingCollision;
        vol_size = new[] { 5, 5, 5 };
        vol_res = new[] { 64, 64, 64 };
        hbar = 0.05f;
        dt = 1f / 24f;
        _velocity = Vector3.zero;
        float m = Mathf.Min(vol_size[0], vol_size[1], vol_size[2]);
        float approach = m * 0.22f;
        var C = new Vector3(vol_size[0] * 0.5f, vol_size[1] * 0.5f, vol_size[2] * 0.5f);
        var nFromLeft = new Vector3(-1f, 0f, 0f);
        var nFromHighZ = new Vector3(0f, 0f, 1f);
        _obliqueRing1Center = C - approach * new Vector3(1f, 0f, 0f);
        _obliqueRing1Normal = nFromLeft;
        _obliqueRing2Center = C + approach * new Vector3(0f, 0f, 1f);
        _obliqueRing2Normal = nFromHighZ;
        _obliqueRingRadius = 0.6f;
        _obliquePsi2Re = 0.05f;
        _useLES = false;
        _stepsPerFrame = 1;
        _nParticles = 100000;
        _particleSize = 0.08f;
        SyncParticleDisplayScaleFromSimulation();
    }

    public void ApplySmokeMaze2DPreset(SFUnifiedSmokeMaze2DPreset p)
    {
        vol_size = (int[])p.vol_size?.Clone() ?? new[] { 5, 1, 3 };
        vol_res = (int[])p.vol_res?.Clone() ?? new[] { 128, 32, 96 };
        hbar = p.hbar;
        dt = p.dt;
        _mazeWallThickness = p.wallThickness;
        _mazeWallMargin = p.wallMargin;
        _mazeWall1X = p.wall1X;
        _mazeWall1GapZ = p.wall1GapZ;
        _mazeWall2X = p.wall2X;
        _mazeWall2SolidZ = p.wall2SolidZ;
        _mazeWall3X = p.wall3X;
        _mazeWall3SolidMaxZ = p.wall3SolidMaxZ;
        _mazeSourceCenter = p.sourceCenter;
        _mazeSourceHalf = p.sourceHalf;
        _mazeEmitVelocity = p.emitVelocity;
        _mazeVentCenter = p.ventCenter;
        _mazeVentHalf = p.ventHalf;
        _mazeVentSuction = p.ventSuction;
        _mazeParticleDispersion = p.particleDispersion;
        _mazeDispersionWallBoost = p.dispersionWallBoost;
        _mazeVentDrift = p.ventDrift;
        _mazeWallDeflect = p.wallDeflect;
        _mazePushSearchCells = p.pushSearchCells;
        _kinematicViscosity = p.kinematicViscosity;
        _nParticles = p.nParticles;
        _particleSize = p.particleSize;
        SyncParticleDisplayScaleFromSimulation();
        _stepsPerFrame = p.stepsPerFrame;
        _useLES = p.useLES;
    }

    [ContextMenu("Apply Smoke Maze 2D defaults")]
    public void ApplySmokeMaze2DDefaults()
    {
        _scenario = ScenarioType.SmokeMaze2D;
        ApplySmokeMaze2DPreset(ResolveScenarioPresets().smokeMaze2D);
        Debug.Log("[SFUnifiedCS] SmokeMaze2D defaults applied.");
    }

    #endregion

    #region Gizmos

    [ExecuteAlways]
    private void OnDrawGizmos()
    {
        var volV3 = new Vector3(vol_size[0], vol_size[1], vol_size[2]);
        Gizmos.color = Color.yellow;
        Gizmos.DrawWireCube(transform.position + volV3 / 2f, volV3);

        switch (_scenario)
        {
            case ScenarioType.SphereObstacle:
                Gizmos.color = new Color(1, 0, 0, 0.5f);
                Gizmos.DrawWireSphere(transform.position + _obstaclePos1, _obstacleRadius1);
                break;
            case ScenarioType.CylinderObstacle:
                Gizmos.color = new Color(1, 0, 0, 0.5f);
                Gizmos.DrawWireSphere(transform.position + _obstaclePos1, _obstacleRadius1);
                break;
            case ScenarioType.TwoSpheres:
                Gizmos.color = new Color(1, 0, 0, 0.5f);
                Gizmos.DrawWireSphere(transform.position + _obstaclePos1, _obstacleRadius1);
                Gizmos.color = new Color(0, 0, 1, 0.5f);
                Gizmos.DrawWireSphere(transform.position + _obstaclePos2, _obstacleRadius2);
                break;
            case ScenarioType.InkCollision:
                Gizmos.color = new Color(0.9f, 0.2f, 0.2f, 0.55f);
                Gizmos.DrawWireSphere(transform.position + _obstaclePos1, _obstacleRadius1);
                Gizmos.color = new Color(0.3f, 0.6f, 1f, 0.55f);
                Gizmos.DrawWireSphere(transform.position + _obstaclePos2, _obstacleRadius2);
                break;
            case ScenarioType.ObliqueRingCollision:
                Gizmos.color = new Color(0.95f, 0.25f, 0.2f, 0.65f);
                Gizmos.DrawWireSphere(transform.position + _obliqueRing1Center, _obliqueRingRadius);
                Gizmos.color = new Color(0.2f, 0.45f, 1f, 0.65f);
                Gizmos.DrawWireSphere(transform.position + _obliqueRing2Center, _obliqueRingRadius);
                break;
            case ScenarioType.Cigarette:
                Gizmos.color = new Color(0.2f, 0.8f, 0.3f, 0.6f);
                Gizmos.DrawWireSphere(transform.position + _cigaretteHeatSphereCen,
                    _cigaretteHeatSphereRad);
                break;
            case ScenarioType.SmokeMaze2D:
                DrawSmokeMazeGizmos();
                break;
        }
    }

    private void DrawSmokeMazeGizmos()
    {
        var o = transform.position;
        float wy = vol_size[1];
        float t = _mazeWallThickness;

        Gizmos.color = new Color(1f, 0.55f, 0.1f, 0.85f);
        Gizmos.DrawWireCube(o + _mazeSourceCenter, _mazeSourceHalf * 2f);

        Gizmos.color = new Color(0.2f, 0.85f, 0.35f, 0.85f);
        Gizmos.DrawWireCube(o + _mazeVentCenter, _mazeVentHalf * 2f);

        Gizmos.color = new Color(0.25f, 0.45f, 1f, 0.75f);
        DrawMazeSlatSegment(_mazeWall1X, 0f, _mazeWall1GapZ.x, wy, t);
        DrawMazeSlatSegment(_mazeWall1X, _mazeWall1GapZ.y, vol_size[2], wy, t);
        DrawMazeSlatSegment(_mazeWall2X, _mazeWall2SolidZ.x, _mazeWall2SolidZ.y, wy, t);
        DrawMazeSlatSegment(_mazeWall3X, 0f, _mazeWall3SolidMaxZ, wy, t);

        float m = _mazeWallMargin;
        float wx = vol_size[0], wz = vol_size[2];
        DrawMazeSlatSegment(m, 0f, wz, wy, m);
        DrawMazeSlatSegment(wx - m, 0f, wz, wy, m);
        Gizmos.DrawWireCube(o + new Vector3(wx * 0.5f, wy * 0.5f, m * 0.5f), new Vector3(wx, wy, m));
        Gizmos.DrawWireCube(o + new Vector3(wx * 0.5f, wy * 0.5f, wz - m * 0.5f), new Vector3(wx, wy, m));
    }

    private void DrawMazeSlatSegment(float centerX, float z0, float z1, float heightY, float thicknessX)
    {
        if (z1 <= z0) return;
        var o = transform.position;
        float cz = (z0 + z1) * 0.5f;
        var center = new Vector3(centerX, heightY * 0.5f, cz);
        var size = new Vector3(thicknessX, heightY, z1 - z0);
        Gizmos.DrawWireCube(o + center, size);
    }

    #endregion

    #region IRaymarchDensitySource

    public bool TryGetPsiVolume(out ComputeBuffer psi1, out ComputeBuffer psi2,
        out Vector3 volumeMinWorld, out Vector3 volumeSizeWorld,
        out int resX, out int resY, out int resZ)
    {
        psi1 = null;
        psi2 = null;
        volumeMinWorld = default;
        volumeSizeWorld = default;
        resX = resY = resZ = 0;

        if (!_initialized || _isf == null)
            return false;

        psi1 = _isf.psi1;
        psi2 = _isf.psi2;
        volumeMinWorld = transform.position;
        volumeSizeWorld = new Vector3(vol_size[0], vol_size[1], vol_size[2]);
        resX = _isf.resX;
        resY = _isf.resY;
        resZ = _isf.resZ;
        return true;
    }

    #endregion
}
