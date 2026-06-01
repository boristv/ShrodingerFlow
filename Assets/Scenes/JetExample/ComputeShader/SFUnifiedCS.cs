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
        /// <summary>Домен = внутренность ёмкости; стенки — penalization в полосе у границ; начальное состояние — прямой блок жидкости. Гравитация и вязкость — через существующие поля.</summary>
        RectangularContainer,
        /// <summary>Цифровой двойник распространения дыма: помещение из 2 комнат с проёмом (стены — penalization), источник дыма (инжекция χ + плавучесть), вытяжка (сток χ + подсос + удаление трассеров).</summary>
        RoomSmoke
    }

    [Header("Compute Shaders")]
    [SerializeField] private ComputeShader _kernelsShader;
    [Tooltip("ISF со стенками сетки, χ и AddGravityToVelocity. Только RectangularContainer; для Jet/колец не используется. В Editor подставляется из папки, если пусто.")]
    [SerializeField] private ComputeShader _containerIsfShader;
    [SerializeField] private ComputeShader _fftShader;
    [SerializeField] private ComputeShader _particlesShader;
    [Tooltip("Отдельный compute с ConstrainParticlesToLiquidChi. Для ёмкости+χ в Editor подставляется из папки, если поле пустое; для билда перетащите SFComputeParticlesChiConstrain.")]
    [SerializeField] private ComputeShader _particlesChiConstrainShader;
    [Tooltip("Интерполяция u без wrap + кламп трассеров (SFComputeParticlesWall). Только RectangularContainer; в Editor подставляется из папки, если пусто.")]
    [SerializeField] private ComputeShader _particlesWallShader;
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

    [Header("Прямоугольная ёмкость (RectangularContainer)")]
    [Tooltip("Минимальный угол блока начальной жидкости в координатах объёма [0..vol_size]. Нужны «гравитация на ψ₂» и клип интерполяции скорости (включено в коде для этого сценария).")]
    [SerializeField] private Vector3 _containerFluidMin = new Vector3(0.35f, 2.2f, 0.4f);
    [Tooltip("Максимальный угол блока начальной жидкости.")]
    [SerializeField] private Vector3 _containerFluidMax = new Vector3(2.65f, 3.55f, 2.6f);
    [Tooltip("Толщина слоя ячеек-«стенок» у границ домена (penalization), в тех же единицах, что vol_size.")]
    [SerializeField] private float _containerWallThickness = 0.1f;
    [Tooltip("Ёмкость: разброс трассы после клампа, в долях min(Δx,Δy,Δz). Большие значения дают заметный джиттер у границы домена.")]
    [SerializeField] private float _containerParticleTracerJitter = 0.08f;
    [Tooltip("Поле χ (газ/жидкость): отдельная от |ψ| транспортировка и «вакуум» ψ в газе; только RectangularContainer.")]
    [SerializeField] private bool _useLiquidChiField;
    [Tooltip("Ячейка — жидкость, если χ ≥ порога; иначе после каждой нормировки/фазы ψ сбрасывается к вакууму.")]
    [SerializeField, Range(0f, 1f)] private float _liquidChiThreshold = 0.5f;
    [Tooltip("Подтяжка трассеров к χ: 0 = только поле скорости (стабильнее). Включайте >0 если нужны маркеры строго в жидкости.")]
    [SerializeField, Range(0f, 1f)] private float _liquidParticleConstrainStrength = 0f;
    [Tooltip("Мягкая зона у порога χ: пока сэмпл χ ≥ (порог − margin), подтяжка не включается — меньше скачков на границе жидкости.")]
    [SerializeField, Range(0.02f, 0.25f)] private float _liquidParticleChiSoftMargin = 0.1f;

    [Header("Помещение с дымом (RoomSmoke)")]
    [Tooltip("Толщина внешних стен помещения (penalization), в единицах vol_size.")]
    [SerializeField] private float _roomWallThickness = 0.12f;
    [Tooltip("Координата X внутренней перегородки между двумя комнатами.")]
    [SerializeField] private float _roomPartitionX = 2f;
    [Tooltip("Полутолщина перегородки по X.")]
    [SerializeField] private float _roomPartitionThickness = 0.12f;
    [Tooltip("Центр дверного проёма по Z.")]
    [SerializeField] private float _roomDoorCenterZ = 2f;
    [Tooltip("Ширина дверного проёма по Z.")]
    [SerializeField] private float _roomDoorWidth = 1f;
    [Tooltip("Высота дверного проёма от пола (Y).")]
    [SerializeField] private float _roomDoorHeight = 1.5f;

    [Tooltip("Центр зоны источника дыма (комната A), координаты объёма.")]
    [SerializeField] private Vector3 _roomSourceCenter = new Vector3(0.9f, 0.45f, 2f);
    [Tooltip("Полуразмер зоны источника дыма (AABB).")]
    [SerializeField] private Vector3 _roomSourceHalf = new Vector3(0.22f, 0.22f, 0.22f);
    [Tooltip("Начальная скорость выброса дыма из источника (горячий выброс вверх). k = v/ℏ.")]
    [SerializeField] private Vector3 _roomEmitVelocity = new Vector3(0f, 0.6f, 0f);
    [Tooltip("Значение концентрации χ, нагнетаемое в источнике каждый шаг.")]
    [SerializeField, Range(0f, 1f)] private float _roomChiInject = 1f;

    [Tooltip("Центр зоны вытяжки (комната B, у потолка), координаты объёма.")]
    [SerializeField] private Vector3 _roomVentCenter = new Vector3(3.1f, 2.5f, 2f);
    [Tooltip("Полуразмер зоны вытяжки (AABB).")]
    [SerializeField] private Vector3 _roomVentHalf = new Vector3(0.4f, 0.22f, 0.5f);
    [Tooltip("Скорость подсоса вытяжки (втягивающее граничное условие на u). k = v/ℏ.")]
    [SerializeField] private Vector3 _roomVentSuction = new Vector3(0f, 1.2f, 0f);
    [Tooltip("Доля дыма, остающаяся в зоне вытяжки за шаг (сток χ): меньше = сильнее вытягивает.")]
    [SerializeField, Range(0f, 1f)] private float _roomVentDecay = 0.8f;

    [Tooltip("Коэффициент плавучести β: u += β·χ·dir·dt (горячий дым легче воздуха → подъём).")]
    [SerializeField] private float _roomBuoyancyBeta = 6f;
    [Tooltip("Направление плавучести (обычно +Y).")]
    [SerializeField] private Vector3 _roomBuoyancyDir = new Vector3(0f, 1f, 0f);
    [Tooltip("Скорость всплытия дыма относительно воздуха (drift-flux): несёт χ и трассеры вверх сквозь спокойный воздух. Разрывает «блокировку» маски скорости.")]
    [SerializeField] private float _roomSmokeRiseSpeed = 0.6f;
    [Tooltip("Турбулентная диффузия дыма χ за шаг (alpha = D·dt/h²). Расширяет султан и даёт растекание под потолком/через проём. Стабильно до ~0.16.")]
    [SerializeField, Range(0f, 0.16f)] private float _roomSmokeDiffusion = 0.08f;
    [Tooltip("Турбулентная дисперсия трассеров внутри дыма (в долях ячейки): случайное блуждание, расширяющее столб в клубящееся облако.")]
    [SerializeField, Range(0f, 1f)] private float _roomTracerDispersion = 0.15f;
    [Tooltip("Амплитуда curl-noise турбулентности (вихревое поле скорости внутри дыма): клубление и боковое вовлечение, столб перестаёт быть прямым.")]
    [SerializeField] private float _roomTurbAmplitude = 0.6f;
    [Tooltip("Пространственная частота curl-noise (1/world): крупнее → крупные вихри, мельче → мелкая турбулентность.")]
    [SerializeField] private float _roomTurbScale = 1.6f;
    private float _curlTime;

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
    [Tooltip("Для RectangularContainer при включении гравитации фаза добавляется и на ψ₁ (иначе доминирующий |ψ₁| почти не даёт скорости). Остальные сценарии — только ψ₂, как в hip.")]
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
    private ComputeBuffer _roomWallMaskBuf, _roomSourceMaskBuf, _roomVentMaskBuf;
    private Vector3 _roomVentMin, _roomVentMax;

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

    /// <summary>Сценарии, которым нужен расширенный ISF (стенки сетки, χ, среднее u) и стеночный шейдер частиц.</summary>
    private bool UsesContainerExt =>
        _scenario == ScenarioType.RectangularContainer || _scenario == ScenarioType.RoomSmoke;

#if UNITY_EDITOR
    private const string ParticlesChiConstrainAssetPath =
        "Assets/Scenes/JetExample/ComputeShader/SFComputeParticlesChiConstrain.compute";
    private const string ContainerIsfAssetPath =
        "Assets/Scenes/JetExample/ComputeShader/SFComputeKernelsContainer.compute";
    private const string ParticlesWallAssetPath =
        "Assets/Scenes/JetExample/ComputeShader/SFComputeParticlesWall.compute";
#endif

    /// <summary>Подтяжка трассеров к χ: отдельный compute-asset, не смешиваем с основным шейдером частиц (Metal/CB).</summary>
    private ComputeShader ResolveParticlesChiConstrainShader()
    {
        if (_particlesChiConstrainShader != null)
            return _particlesChiConstrainShader;
#if UNITY_EDITOR
        if ((_scenario == ScenarioType.RectangularContainer && _useLiquidChiField)
            || _scenario == ScenarioType.RoomSmoke)
            return AssetDatabase.LoadAssetAtPath<ComputeShader>(ParticlesChiConstrainAssetPath);
#endif
        return null;
    }

    private ComputeShader ResolveContainerIsfShader()
    {
        if (_containerIsfShader != null)
            return _containerIsfShader;
#if UNITY_EDITOR
        if (UsesContainerExt)
            return AssetDatabase.LoadAssetAtPath<ComputeShader>(ContainerIsfAssetPath);
#endif
        return null;
    }

    private ComputeShader ResolveParticlesWallShader()
    {
        if (_particlesWallShader != null)
            return _particlesWallShader;
#if UNITY_EDITOR
        if (UsesContainerExt)
            return AssetDatabase.LoadAssetAtPath<ComputeShader>(ParticlesWallAssetPath);
#endif
        return null;
    }

    #region Lifecycle

    private void Start()
    {
        _particleBuffers = GetComponent<ParticleGpuBuffers>();
        _particleDisplay = GetComponent<ParticleDisplay3D>();

        _isf = new CSISF();
        bool roomSmoke = _scenario == ScenarioType.RoomSmoke;
        ComputeShader containerIsf = UsesContainerExt ? ResolveContainerIsfShader() : null;
        _isf.Init(_kernelsShader, _fftShader, _lesShader, vol_size, vol_res, hbar, dt, containerIsf);
        _isf.clampGridBorders = UsesContainerExt;
        _isf.useLiquidChiField =
            (_scenario == ScenarioType.RectangularContainer && _useLiquidChiField) || roomSmoke;
        _isf.liquidChiThreshold = _liquidChiThreshold;
        // Дым = χ-«фаза» (как жидкость в контейнере): воздух (χ<порога) → ψ вакуум и u=0 (спокоен, без FFT-шума),
        // дым (χ≥порога) когерентно поднимается плавучестью. Без этого подавления резкие границы разносятся FFT в шум.
        _isf.chiAffectsPsiVacuum = true;
        _isf.maskVelocityWithChi = true;

        bool oneTimeParticles = _scenario == ScenarioType.LeapfrogRings
                             || _scenario == ScenarioType.ObliqueRingCollision
                             || _scenario == ScenarioType.TwoSpheres
                             || _scenario == ScenarioType.InkCollision
                             || _scenario == ScenarioType.RectangularContainer;
        int maxParticles = oneTimeParticles ? _nParticles : _nParticles * 1000;

        _particles = new CSParticles();
        _particles.Init(_particlesShader, maxParticles, _isf, ResolveParticlesChiConstrainShader(),
            ResolveParticlesWallShader());

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

#if UNITY_EDITOR
        if (UsesContainerExt && _containerIsfShader == null)
            _containerIsfShader = AssetDatabase.LoadAssetAtPath<ComputeShader>(ContainerIsfAssetPath);
        if (UsesContainerExt && _particlesWallShader == null)
            _particlesWallShader = AssetDatabase.LoadAssetAtPath<ComputeShader>(ParticlesWallAssetPath);
#endif

        if (_initialized && _isf != null)
        {
            bool roomSmoke = _scenario == ScenarioType.RoomSmoke;
            _isf.useLiquidChiField =
                (_scenario == ScenarioType.RectangularContainer && _useLiquidChiField) || roomSmoke;
            _isf.liquidChiThreshold = _liquidChiThreshold;
            _isf.chiAffectsPsiVacuum = true;
            _isf.maskVelocityWithChi = true;
        }
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

            // Подтяжку к χ по χ делаем раз за кадр: при stepsPerFrame>1 вызов после каждого подшага
            // суммировался и давал сильный дребезг трассеров у границы жидкости.
            if (_scenario == ScenarioType.RectangularContainer && _isf.useLiquidChiField)
                ApplyLiquidChiParticleConstrain();

            if (_spawnEachStep
                && (_scenario == ScenarioType.Jet || _scenario == ScenarioType.Cigarette
                    || _scenario == ScenarioType.RoomSmoke))
            {
                // RoomSmoke компактим чаще: трассеры гибнут в вытяжке и должны быстро освобождать место.
                int compactEvery = _scenario == ScenarioType.RoomSmoke ? 12 : 60;
                _compactCounter++;
                if (_compactCounter >= compactEvery)
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
        _roomWallMaskBuf?.Release();
        _roomSourceMaskBuf?.Release();
        _roomVentMaskBuf?.Release();
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

            case ScenarioType.RectangularContainer:
                InitPsiFluidBlock();
                _maskBuf1 = BuildContainerWallMask();
                _kvecX = _kvecY = _kvecZ = 0f;
                _omega = 0f;
                RunInitBoundary(_maskBuf1, 0f, 0f, 0f, 0f, 10);
                if (_isf.useLiquidChiField)
                    UploadInitialLiquidChiForContainer();
                SpawnParticlesInFluidBlock(_nParticles);
                _particlesCount = _particles.Size;
                _spawnEachStep = false;
                _boundaryEachStep = true;
                break;

            case ScenarioType.RoomSmoke:
                InitPsiUniform();
                _roomWallMaskBuf = BuildRoomLayoutMask();
                _roomSourceMaskBuf = BuildBoxMask(_roomSourceCenter, _roomSourceHalf);
                _roomVentMaskBuf = BuildBoxMask(_roomVentCenter, _roomVentHalf);
                _roomVentMin = _roomVentCenter - _roomVentHalf;
                _roomVentMax = _roomVentCenter + _roomVentHalf;
                UploadInitialChiZero();
                _kvecX = _kvecY = _kvecZ = 0f;
                _omega = 0f;
                RunInitBoundary(_roomWallMaskBuf, 0f, 0f, 0f, 0f, 10);
                _particlesCount = 0;
                _spawnEachStep = true;
                _boundaryEachStep = true;
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

    private void GetFluidBlockBounds(out Vector3 min, out Vector3 max)
    {
        min = Vector3.Min(_containerFluidMin, _containerFluidMax);
        max = Vector3.Max(_containerFluidMin, _containerFluidMax);
        var vmax = new Vector3(vol_size[0], vol_size[1], vol_size[2]);
        min = Vector3.Max(min, Vector3.zero);
        max = Vector3.Min(max, vmax);
    }

    /// <summary>Отрицательно внутри AABB; снаружи — расстояние до ближайшей точки на поверхности.</summary>
    private static float BoxSignedDistanceToFluidAabb(float px, float py, float pz, Vector3 bmin, Vector3 bmax)
    {
        float cx = Mathf.Clamp(px, bmin.x, bmax.x);
        float cy = Mathf.Clamp(py, bmin.y, bmax.y);
        float cz = Mathf.Clamp(pz, bmin.z, bmax.z);
        float dx = px - cx, dy = py - cy, dz = pz - cz;
        float dOut = Mathf.Sqrt(dx * dx + dy * dy + dz * dz);
        if (dOut > 1e-8f)
            return dOut;
        float dInX = Mathf.Min(px - bmin.x, bmax.x - px);
        float dInY = Mathf.Min(py - bmin.y, bmax.y - py);
        float dInZ = Mathf.Min(pz - bmin.z, bmax.z - pz);
        return -Mathf.Min(dInX, Mathf.Min(dInY, dInZ));
    }

    /// <summary>
    /// Неподвижная жидкость в прямоугольнике (геометрия — для частиц/χ); ψ задаётся согласованно по всему домену.
    /// Normalize делается по ячейке — модули сравниваются только внутри (ψ₁,ψ₂); снаружи нужна та же пропорция,
    /// что внутри блока, иначе скачок составляющих даёт ложную скорость на границе (в т.ч. «улетает в сторону») без всякого χ.
    /// </summary>
    private void InitPsiFluidBlock()
    {
        GetFluidBlockBounds(out Vector3 fmin, out Vector3 fmax);
        int num = _isf.num;
        var in1 = new Vector2(1f, 0f);
        var in2 = new Vector2(0.01f, 0f);
        const float outsideScale = 1e-6f;
        var tmp1 = new Vector2[num];
        var tmp2 = new Vector2[num];
        float feather = 3f * Mathf.Min(_isf.dx, Mathf.Min(_isf.dy, _isf.dz));
        feather = Mathf.Max(feather, 1e-6f);
        for (int i = 0; i < num; i++)
        {
            float px = _isf.pxCPU[i], py = _isf.pyCPU[i], pz = _isf.pzCPU[i];
            float sd = BoxSignedDistanceToFluidAabb(px, py, pz, fmin, fmax);
            float win = sd <= 0f ? 1f : Mathf.Exp(-sd / feather);
            tmp1[i] = Vector2.Lerp(in1 * outsideScale, in1, win);
            tmp2[i] = Vector2.Lerp(in2 * outsideScale, in2, win);
        }
        _isf.psi1.SetData(tmp1);
        _isf.psi2.SetData(tmp2);
        _isf.Normalize();
    }

    /// <summary>χ=1 в том же AABB, что начальный блок жидкости; χ=0 в газе (несжимаемая «пустота» в терминах ψ задаётся яхром в CSISF).</summary>
    private void UploadInitialLiquidChiForContainer()
    {
        GetFluidBlockBounds(out Vector3 fmin, out Vector3 fmax);
        int n = _isf.num;
        var chi = new float[n];
        for (int i = 0; i < n; i++)
        {
            float px = _isf.pxCPU[i], py = _isf.pyCPU[i], pz = _isf.pzCPU[i];
            bool inside = px >= fmin.x && px <= fmax.x
                       && py >= fmin.y && py <= fmax.y
                       && pz >= fmin.z && pz <= fmax.z;
            chi[i] = inside ? 1f : 0f;
        }
        _isf.UploadLiquidChi(chi);
    }

    /// <summary>RoomSmoke: чистый воздух — χ=0 по всему объёму (дым появляется только из источника).</summary>
    private void UploadInitialChiZero()
    {
        var chi = new float[_isf.num];
        _isf.UploadLiquidChi(chi);
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

    /// <summary>Маска стенок ёмкости: ячейки в полосе у границы домена [0, vol_size] (как у твёрдого препятствия, k=0).</summary>
    private ComputeBuffer BuildContainerWallMask()
    {
        int num = _isf.num;
        var mask = new int[num];
        float t = Mathf.Max(0f, _containerWallThickness);
        float wx = vol_size[0], wy = vol_size[1], wz = vol_size[2];
        for (int i = 0; i < num; i++)
        {
            float px = _isf.pxCPU[i], py = _isf.pyCPU[i], pz = _isf.pzCPU[i];
            bool wall = px < t || px > wx - t
                     || py < t || py > wy - t
                     || pz < t || pz > wz - t;
            mask[i] = wall ? 1 : 0;
        }
        var buf = new ComputeBuffer(num, sizeof(int));
        buf.SetData(mask);
        return buf;
    }

    /// <summary>Маска AABB (1 внутри коробки center±half) в координатах объёма.</summary>
    private ComputeBuffer BuildBoxMask(Vector3 center, Vector3 half)
    {
        int num = _isf.num;
        var mask = new int[num];
        Vector3 bmin = center - half;
        Vector3 bmax = center + half;
        for (int i = 0; i < num; i++)
        {
            float px = _isf.pxCPU[i], py = _isf.pyCPU[i], pz = _isf.pzCPU[i];
            bool inside = px >= bmin.x && px <= bmax.x
                       && py >= bmin.y && py <= bmax.y
                       && pz >= bmin.z && pz <= bmax.z;
            mask[i] = inside ? 1 : 0;
        }
        var buf = new ComputeBuffer(num, sizeof(int));
        buf.SetData(mask);
        return buf;
    }

    /// <summary>
    /// Маска твёрдого тела помещения: внешние стены (полоса у границ домена) + внутренняя перегородка по X
    /// с дверным проёмом (z вокруг центра, y ниже высоты двери). Из солида вычитаются зоны источника и вытяжки,
    /// чтобы они оставались открытыми.
    /// </summary>
    private ComputeBuffer BuildRoomLayoutMask()
    {
        int num = _isf.num;
        var mask = new int[num];
        float t = Mathf.Max(0f, _roomWallThickness);
        float wx = vol_size[0], wy = vol_size[1], wz = vol_size[2];
        float halfPart = Mathf.Max(0f, _roomPartitionThickness);
        float doorZmin = _roomDoorCenterZ - _roomDoorWidth * 0.5f;
        float doorZmax = _roomDoorCenterZ + _roomDoorWidth * 0.5f;

        Vector3 srcMin = _roomSourceCenter - _roomSourceHalf;
        Vector3 srcMax = _roomSourceCenter + _roomSourceHalf;
        Vector3 ventMin = _roomVentCenter - _roomVentHalf;
        Vector3 ventMax = _roomVentCenter + _roomVentHalf;

        for (int i = 0; i < num; i++)
        {
            float px = _isf.pxCPU[i], py = _isf.pyCPU[i], pz = _isf.pzCPU[i];

            bool border = px < t || px > wx - t
                       || py < t || py > wy - t
                       || pz < t || pz > wz - t;

            bool partition = Mathf.Abs(px - _roomPartitionX) <= halfPart;
            bool door = pz >= doorZmin && pz <= doorZmax && py <= _roomDoorHeight;
            bool solid = border || (partition && !door);

            if (solid)
            {
                bool inSource = px >= srcMin.x && px <= srcMax.x
                             && py >= srcMin.y && py <= srcMax.y
                             && pz >= srcMin.z && pz <= srcMax.z;
                bool inVent = px >= ventMin.x && px <= ventMax.x
                           && py >= ventMin.y && py <= ventMax.y
                           && pz >= ventMin.z && pz <= ventMax.z;
                if (inSource || inVent)
                    solid = false;
            }

            mask[i] = solid ? 1 : 0;
        }
        var buf = new ComputeBuffer(num, sizeof(int));
        buf.SetData(mask);
        return buf;
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
        if (_scenario == ScenarioType.RoomSmoke)
        {
            SimulationStepRoomSmoke();
            return;
        }

        _isf.kinematicViscosity = _kinematicViscosity;
        _isf.clampGridBorders = _scenario == ScenarioType.RectangularContainer;
        _isf.liquidChiThreshold = _liquidChiThreshold;
        bool isContainer = _scenario == ScenarioType.RectangularContainer;
        if (_scenario == ScenarioType.Cigarette)
            _isf.UpdateCigaretteSpace(_useLES, _cigaretteGravity, _maskBuf1);
        else
        {
            // Для контейнера фазовую гравитацию на ψ не используем — сила задаётся в ApplyVelocityGravity
            // (иначе накопление фазы → «перемотка» atan2 в VelocityOne → скачки).
            Vector3? phaseGrav = (_applyPsi2Gravity && !isContainer) ? _psi2Gravity : (Vector3?)null;
            _isf.UpdateSpace(_useLES, phaseGrav, false);
        }

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
        // Gravity ДО адвекции χ: поле скорости должно нести гравитацию,
        // иначе χ advect-ится по нулевому полю и блок жидкости остаётся на месте.
        if (isContainer && _applyPsi2Gravity && _psi2Gravity.sqrMagnitude > 1e-14f)
            _isf.ApplyVelocityGravity(_psi2Gravity, _vel);
        if (_isf.useLiquidChiField)
        {
            _isf.AdvectLiquidChi(_vel);
            _isf.ApplyLiquidChiVelocityMask(_vel);
        }
        _particles.CalculateMovement(_vel, isContainer);

        if (_scenario == ScenarioType.RectangularContainer)
        {
            _particles.ClampPositionsToVolume(vol_size[0], vol_size[1], vol_size[2],
                _containerParticleTracerJitter, iterator);
        }
        else if (_scenario != ScenarioType.Jet && _scenario != ScenarioType.Cigarette
            && _scenario != ScenarioType.ObliqueRingCollision)
            _particles.WrapPositions(vol_size[0], vol_size[1], vol_size[2]);
    }

    private void ApplyLiquidChiParticleConstrain()
    {
        float mcell = Mathf.Min(_isf.dx, Mathf.Min(_isf.dy, _isf.dz));
        Vector3 gdir = Vector3.zero;
        float gPerCell = 0f;
        if (_applyPsi2Gravity && _psi2Gravity.sqrMagnitude > 1e-14f)
        {
            gdir = _psi2Gravity.normalized;
            gPerCell = 0.45f / Mathf.Max(mcell, 1e-6f);
        }

        _particles.ConstrainToLiquidChi(_isf.LiquidChiBuffer, _liquidChiThreshold,
            _liquidParticleConstrainStrength, gdir, gPerCell, _liquidParticleChiSoftMargin);
    }

    /// <summary>
    /// Шаг цифрового двойника дыма: ISF → стены/перегородка (solid) + источник (jet ↑) + вытяжка (подсос) →
    /// инжекция χ и спавн трассеров → плавучесть(χ) → адвекция χ и сток в вытяжке → адвекция трассеров,
    /// выталкивание из стен, кламп к объёму, удаление в вытяжке.
    /// </summary>
    private void SimulationStepRoomSmoke()
    {
        _isf.kinematicViscosity = _kinematicViscosity;
        _isf.clampGridBorders = true;

        _isf.UpdateSpace(_useLES, null, false);

        // Граничные условия на скорость: стены k=0 (непротекание), вытяжка — подсос.
        // Источник НЕ задаём струёй (давало «лазерный» столб) — дым поднимается плавучестью + дрейфом.
        float invH = 1f / hbar;
        _isf.ApplyJetBoundary(_roomWallMaskBuf, 0f, 0f, 0f, 0f);
        _isf.ApplyJetBoundary(_roomVentMaskBuf,
            _roomVentSuction.x * invH, _roomVentSuction.y * invH, _roomVentSuction.z * invH, 0f);
        _isf.PressureProject();

        SpawnRoomSourceParticles();
        _isf.InjectChi(_roomSourceMaskBuf, _roomChiInject);

        _isf.UpdateVelocities(_vel);
        if (_roomBuoyancyBeta > 1e-6f)
            _isf.ApplyBuoyancy(_roomBuoyancyDir, _roomBuoyancyBeta, _vel);

        // Когерентная вихревая турбулентность (curl-noise) в дыму: клубление + боковое вовлечение,
        // султан расширяется и перестаёт быть прямым. Анимируется во времени (эволюция вихрей).
        if (_roomTurbAmplitude > 1e-5f)
        {
            _curlTime += dt;
            float thrLo = Mathf.Max(0f, _liquidChiThreshold - 0.25f);
            Vector3 tOff = new Vector3(_curlTime * 0.6f, _curlTime * 0.9f, _curlTime * 0.45f);
            _isf.AddCurlTurbulence(_vel, _roomTurbAmplitude, _roomTurbScale, tOff,
                thrLo, _liquidChiThreshold);
        }

        // Drift-flux: дым всплывает относительно воздуха — несём χ вверх даже там, где поле u замаскировано.
        Vector3 riseDrift = _roomBuoyancyDir.normalized * _roomSmokeRiseSpeed;
        _isf.AdvectLiquidChi(_vel, riseDrift);
        // Турбулентная диффузия χ: расширяет султан вбок и растекает дым под потолком/через проём
        // (маска скорости держит воздух спокойным, поэтому без диффузии дым шёл бы тонким столбом).
        _isf.DiffuseChi(_roomSmokeDiffusion);
        _isf.VentChiSink(_roomVentMaskBuf, _roomVentDecay);
        // Воздух (низкий χ) — неподвижен: гасит FFT-шум и удерживает трассеры в дыму.
        _isf.ApplyLiquidChiVelocityMask(_vel);

        _particles.CalculateMovement(_vel, true);
        // Тот же дрейф для трассеров: поднимаются вместе с дымом сквозь спокойный воздух.
        float thr = _liquidChiThreshold;
        _particles.AddBuoyantDrift(_isf.LiquidChiBuffer, riseDrift, dt,
            Mathf.Max(0f, thr - 0.2f), thr);
        // Турбулентная дисперсия: случайное блуждание трассеров внутри дыма — столб → клубящееся облако.
        if (_roomTracerDispersion > 1e-4f)
        {
            float mcell = Mathf.Min(_isf.dx, Mathf.Min(_isf.dy, _isf.dz));
            _particles.DisperseInChi(_isf.LiquidChiBuffer,
                _roomTracerDispersion * mcell, iterator, Mathf.Max(0.05f, thr - 0.35f));
        }
        _particles.PushOutOfSolid(_roomWallMaskBuf, 4);
        _particles.ClampPositionsToVolume(vol_size[0], vol_size[1], vol_size[2],
            _containerParticleTracerJitter, iterator);
        _particles.KillInVent(_roomVentMin, _roomVentMax);
        // Трассеры вне дыма (рассеялись) — убрать, чтобы не копился статичный «замёрзший» шар.
        _particles.KillLowChi(_isf.LiquidChiBuffer, Mathf.Max(0.05f, thr - 0.35f));
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

    /// <summary>Непрерывная эмиссия трассеров дыма в зоне источника (кольцевой буфер: старые перезаписываются).</summary>
    private void SpawnRoomSourceParticles()
    {
        int count = _nParticles;
        if (count <= 0) return;
        Vector3 bmin = _roomSourceCenter - _roomSourceHalf;
        Vector3 bmax = _roomSourceCenter + _roomSourceHalf;
        var xx = new float[count];
        var yy = new float[count];
        var zz = new float[count];
        for (int i = 0; i < count; i++)
        {
            xx[i] = Random.Range(bmin.x, bmax.x);
            yy[i] = Random.Range(bmin.y, bmax.y);
            zz[i] = Random.Range(bmin.z, bmax.z);
        }
        _particles.AddParticles(xx, yy, zz, count, ring: true);
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

    private void SpawnParticlesInFluidBlock(int count)
    {
        if (count <= 0) return;
        GetFluidBlockBounds(out Vector3 fmin, out Vector3 fmax);
        if (fmax.x <= fmin.x || fmax.y <= fmin.y || fmax.z <= fmin.z)
        {
            Debug.LogWarning("[SFUnifiedCS] RectangularContainer: блок жидкости вырожден, частицы не заспавнены.");
            return;
        }
        var xx = new float[count];
        var yy = new float[count];
        var zz = new float[count];
        for (int i = 0; i < count; i++)
        {
            xx[i] = Random.Range(fmin.x, fmax.x);
            yy[i] = Random.Range(fmin.y, fmax.y);
            zz[i] = Random.Range(fmin.z, fmax.z);
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
            || _scenario == ScenarioType.RoomSmoke;
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

    public void ApplyRectangularContainerPreset(SFUnifiedRectangularContainerPreset p)
    {
        _scenario = ScenarioType.RectangularContainer;
        vol_size = (int[])p.vol_size?.Clone() ?? new[] { 3, 4, 3 };
        vol_res = (int[])p.vol_res?.Clone() ?? new[] { 64, 64, 64 };
        hbar = p.hbar;
        dt = p.dt;
        _containerFluidMin = p.fluidMin;
        _containerFluidMax = p.fluidMax;
        _containerWallThickness = p.wallThickness;
        _applyPsi2Gravity = p.applyPsi2Gravity;
        _psi2Gravity = p.psi2Gravity;
        _kinematicViscosity = p.kinematicViscosity;
        _nParticles = p.nParticles;
        _particleSize = p.particleSize;
        SyncParticleDisplayScaleFromSimulation();
        _stepsPerFrame = p.stepsPerFrame;
        _useLES = p.useLES;
        _useLiquidChiField = p.useLiquidChiField;
        _liquidChiThreshold = p.liquidChiThreshold;
    }

    /// <summary>Параметры из SFUnifiedScenarioPresets.rectangularContainer (или встроенные при отсутствии ассета).</summary>
    [ContextMenu("Apply Rectangular container defaults")]
    public void ApplyRectangularContainerDefaults()
    {
        ApplyRectangularContainerPreset(ResolveScenarioPresets().rectangularContainer);
        Debug.Log("[SFUnifiedCS] RectangularContainer defaults applied.");
    }

    public void ApplyRoomSmokePreset(SFUnifiedRoomSmokePreset p)
    {
        _scenario = ScenarioType.RoomSmoke;
        vol_size = (int[])p.vol_size?.Clone() ?? new[] { 4, 3, 4 };
        vol_res = (int[])p.vol_res?.Clone() ?? new[] { 64, 64, 64 };
        hbar = p.hbar;
        dt = p.dt;
        _roomWallThickness = p.wallThickness;
        _roomPartitionX = p.partitionX;
        _roomPartitionThickness = p.partitionThickness;
        _roomDoorCenterZ = p.doorCenterZ;
        _roomDoorWidth = p.doorWidth;
        _roomDoorHeight = p.doorHeight;
        _roomSourceCenter = p.sourceCenter;
        _roomSourceHalf = p.sourceHalf;
        _roomEmitVelocity = p.emitVelocity;
        _roomChiInject = p.chiInject;
        _roomVentCenter = p.ventCenter;
        _roomVentHalf = p.ventHalf;
        _roomVentSuction = p.ventSuction;
        _roomVentDecay = p.ventDecay;
        _roomBuoyancyBeta = p.buoyancyBeta;
        _roomBuoyancyDir = p.buoyancyDir;
        _roomSmokeRiseSpeed = p.smokeRiseSpeed;
        _roomSmokeDiffusion = p.smokeDiffusion;
        _roomTracerDispersion = p.tracerDispersion;
        _roomTurbAmplitude = p.turbAmplitude;
        _roomTurbScale = p.turbScale;
        _kinematicViscosity = p.kinematicViscosity;
        _nParticles = p.nParticles;
        _particleSize = p.particleSize;
        SyncParticleDisplayScaleFromSimulation();
        _stepsPerFrame = p.stepsPerFrame;
        _useLES = p.useLES;
    }

    /// <summary>Параметры из SFUnifiedScenarioPresets.roomSmoke (или встроенные при отсутствии ассета).</summary>
    [ContextMenu("Apply Room smoke defaults")]
    public void ApplyRoomSmokeDefaults()
    {
        ApplyRoomSmokePreset(ResolveScenarioPresets().roomSmoke);
        Debug.Log("[SFUnifiedCS] RoomSmoke defaults applied.");
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
            case ScenarioType.RectangularContainer:
                GetFluidBlockBounds(out Vector3 fmin, out Vector3 fmax);
                Vector3 fc = (fmin + fmax) * 0.5f;
                Vector3 fsize = fmax - fmin;
                Gizmos.color = new Color(0.25f, 0.55f, 1f, 0.75f);
                Gizmos.DrawWireCube(transform.position + fc, fsize);
                break;
            case ScenarioType.RoomSmoke:
                DrawRoomSmokeGizmos();
                break;
        }
    }

    private void DrawRoomSmokeGizmos()
    {
        var pos = transform.position;
        // Перегородка с проёмом.
        Gizmos.color = new Color(0.7f, 0.7f, 0.7f, 0.5f);
        var partCenter = new Vector3(_roomPartitionX, vol_size[1] * 0.5f, vol_size[2] * 0.5f);
        var partSize = new Vector3(_roomPartitionThickness * 2f, vol_size[1], vol_size[2]);
        Gizmos.DrawWireCube(pos + partCenter, partSize);
        // Дверной проём.
        Gizmos.color = new Color(0.95f, 0.85f, 0.2f, 0.8f);
        var doorCenter = new Vector3(_roomPartitionX, _roomDoorHeight * 0.5f, _roomDoorCenterZ);
        var doorSize = new Vector3(_roomPartitionThickness * 2f, _roomDoorHeight, _roomDoorWidth);
        Gizmos.DrawWireCube(pos + doorCenter, doorSize);
        // Источник дыма.
        Gizmos.color = new Color(0.2f, 0.85f, 0.3f, 0.8f);
        Gizmos.DrawWireCube(pos + _roomSourceCenter, _roomSourceHalf * 2f);
        // Вытяжка.
        Gizmos.color = new Color(0.3f, 0.6f, 1f, 0.85f);
        Gizmos.DrawWireCube(pos + _roomVentCenter, _roomVentHalf * 2f);
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
