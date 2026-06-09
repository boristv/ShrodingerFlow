using System;
using System.IO;
using UnityEngine;
using ComputeShaderSF;
using ShrodingerFlow.Particles;

/// <summary>
/// Гибридная волново-плотностная модель ISF (диссертация Тюлькина): ψ — вихревое ядро, α — поле состояния
/// (плотность дыма / интерфейс / оптика), переносимое скоростью ũ из ψ+LES. v1: ИЗОТЕРМИЧЕСКАЯ 3D-комната.
/// Шаг (гл. 5.1): ISF(ψ)→u→LES→проекция→advect(α)→границы. Силы вводятся через состояние (фазовый импульс на входе).
/// Единые маски (гл. 4.4.9): стены/источник/сток видимы и ψ-границе, и переносу α, и оптике.
/// </summary>
/// <summary>Осесимметричный бокс (центр+полуразмер) для данных плана: перегородка, колонна, вытяжка.</summary>
[System.Serializable]
public struct HybridBox
{
    public Vector3 center;
    public Vector3 half;
}

public class SFHybridCS : SFBase, IRaymarchScalarFieldSource, IRaymarchSolidMaskSource
{
    [Header("Compute Shaders")]
    [SerializeField] private ComputeShader _kernelsShader;
    [SerializeField] private ComputeShader _fftShader;
    [SerializeField] private ComputeShader _particlesShader;
    [SerializeField] private ComputeShader _lesShader;
    [Tooltip("SFHybridScalar.compute — перенос плотностно-фазового поля α.")]
    [SerializeField] private ComputeShader _scalarShader;

    [Header("Базовые параметры ISF")]
    [SerializeField] private int[] vol_size = { 4, 3, 4 };
    [SerializeField] private int[] vol_res = { 64, 32, 64 };
    [SerializeField] private float hbar = 0.05f;
    [SerializeField] private float dt = 1f / 24f;

    [Header("Комната (вид: короб vol_size; стены по периметру с отступом margin)")]
    [SerializeField] private float _wallMargin = 0.12f;
    [Tooltip("Вход дыма (струя): центр и полуразмер AABB на/у стены, + скорость впуска.")]
    [SerializeField] private Vector3 _inletCenter = new Vector3(0.25f, 0.8f, 2.0f);
    [SerializeField] private Vector3 _inletHalf = new Vector3(0.14f, 0.35f, 0.35f);
    [SerializeField] private Vector3 _inletVelocity = new Vector3(0.7f, 0f, 0f);
    [Tooltip("Вытяжка: центр и полуразмер AABB; α откачивается, фаза задаёт отток.")]
    [SerializeField] private Vector3 _ventCenter = new Vector3(3.75f, 2.2f, 2.0f);
    [SerializeField] private Vector3 _ventHalf = new Vector3(0.14f, 0.35f, 0.45f);
    [SerializeField] private Vector3 _ventVelocity = new Vector3(0.5f, 0f, 0f);
    [Tooltip("Вытяжка включена: тяга в зоне вытяжки + откачка α + удаление трассеров. Выкл — вытяжка не работает (дым только накапливается/уходит пассивно). Можно щёлкать в Play.")]
    [SerializeField] private bool _ventEnabled = true;
    [Tooltip("Препятствие-мебель (короб): центр и полуразмер. Нулевой размер — выкл.")]
    [SerializeField] private Vector3 _obstacleCenter = new Vector3(2.0f, 0.9f, 2.0f);
    [SerializeField] private Vector3 _obstacleHalf = new Vector3(0.4f, 0.9f, 0.4f);
    [Tooltip("План этажа: перегородки/колонны (твёрдые боксы). Проёмы = промежутки между боксами. Пусто = простая комната.")]
    public HybridBox[] wallSegments;
    [Tooltip("Доп. вытяжки (помимо основной): откачка α + проём в стене. Пусто = только основная.")]
    public HybridBox[] extraVents;
    [Tooltip("Рисовать внешние стены в раймарче. Обычно выкл: они закрывают обзор. Границей потока остаются в любом случае.")]
    [SerializeField] private bool _renderOuterWalls;
    [Tooltip("Меши перегородок/колонн (полупрозрачные) — чтобы стены было видно в режимах Billboard/Shaded (в раймарче их затирает полноэкранный проход, там работают объёмные стены).")]
    [SerializeField] private bool _showWallMeshes = true;
    [SerializeField] private Color _wallMeshColor = new Color(0.55f, 0.6f, 0.7f, 0.22f);
    [Tooltip("Открыть грани ±X (направление продувки): сквозной поток вместо запечатанной коробки. В закрытом объёме ISF поток застаивается (ветер гасится давлением). Вкл = проветриваемая комната/аэродинам. труба.")]
    [SerializeField] private bool _openFlowFaces = true;

    [Header("Поток (впуск, непрерывная подкачка)")]
    [Tooltip("ВКЛ — впуск = всё входное сечение у −X грани (Дирихле-сквозняк, добивает до дальнего конца, авто-масштаб с размером комнаты). ВЫКЛ — локализованная струя-зона ниже (3D-плюм, но бьёт лишь часть комнаты).")]
    [SerializeField] private bool _inflowFullFace;
    [Tooltip("Глубина входного сечения по X (для полного впуска).")]
    [SerializeField] private float _inflowDepth = 0.5f;
    [Tooltip("Зона впуска (струя): когда полный впуск выключен. Смести центр по Y для асимметрии/вертикали.")]
    [SerializeField] private Vector3 _inflowCenter = new Vector3(0.2f, 1.0f, 2.0f);
    [SerializeField] private Vector3 _inflowHalf = new Vector3(0.35f, 0.6f, 0.8f);
    [SerializeField, Range(1, 8)] private int _boundaryIters = 4;

    [Header("Плотностно-фазовое поле α")]
    [Tooltip("Целевое α в источнике (концентрация дыма на входе).")]
    [SerializeField, Range(0f, 1f)] private float _alphaSourceValue = 1f;
    [Tooltip("Множитель α в вытяжке за шаг (<1 — откачка).")]
    [SerializeField, Range(0f, 1f)] private float _alphaSinkFactor = 0.6f;
    [Tooltip("Регуляризующая диффузия D_α (мягкость интерфейса, заполнение застойных зон).")]
    [SerializeField] private float _alphaDiffusion = 0.0015f;
    [SerializeField, Range(0, 4)] private int _alphaDiffuseIters = 1;

    [Header("Тепло (v3): температура, остывание, плавучесть от T")]
    [Tooltip("Сила плавучести ∝ T (вверх, +Y), через фазу ψ. Горячее у очага всплывает; остывший дым нейтрален → растекается/оседает. 0 — выкл.")]
    [SerializeField] private float _buoyancy = 3f;
    [Tooltip("Диффузия температуры (тепловая) — сглаживает поле T.")]
    [SerializeField] private float _thermalDiffusion = 0.003f;
    [Tooltip("Скорость остывания (1/с): дым теряет тепло по мере распространения. Больше → быстрее теряет подъём (потолочный слой оседает раньше). 0 — дым вечно горячий/всплывает.")]
    [SerializeField] private float _cooling = 0.4f;

    [Header("Турбулентность")]
    [Tooltip("Vorticity confinement: возвращает мелкие завихрения (клубящийся факел), размытые численной диффузией. 0 — выкл. Действует на транспортную скорость, ψ не трогает.")]
    [SerializeField] private float _vorticityConfinement = 3f;
    [Tooltip("MacCormack-перенос α/T (низкодиффузионный): держит концентрацию и вихри → клубление и вовлечение (entrainment) видны, дым 'пухнет' при подъёме. Выкл — обычный полулагранж (размытее).")]
    [SerializeField] private bool _macCormack = true;

    [Header("Трассеры — опциональный A/B-режим (для Billboard/Shaded)")]
    [Tooltip("Считать частицы-трассеры параллельно α. Несутся ТОЙ ЖЕ скоростью ISF+LES. Не влияют на α-поле. Видны в режимах Billboard/Shaded.")]
    [SerializeField] private bool _spawnTracers = true;
    [Tooltip("Сколько трассеров рождать за шаг. МЕНЬШЕ → дольше живёт каждый (медленнее оборот при том же бюджете) → шлейф добивает дальше.")]
    [SerializeField, Range(0, 1000)] private int _tracerPerStep = 30;
    [Tooltip("Бюджет популяции: держим столько одновременно. Уходят в вытяжке/за границей; при переполнении перерабатываются САМЫЕ СТАРЫЕ. Срока по таймеру нет. Эффективная «жизнь» ≈ бюджет / спавн-за-шаг — поэтому большой бюджет + малый спавн = долгий шлейф.")]
    [SerializeField, Range(2000, 400000)] private int _tracerPopulationCap = 120000;
    [Tooltip("Диффузия трассеров D (0 — чистая адвекция по характеристикам: контраст с α).")]
    [SerializeField] private float _tracerDiffusion;

    [Header("Управление")]
    [SerializeField] private bool _useLES = true;
    [SerializeField] private float _kinematicViscosity = 0.00005f;
    [SerializeField] private bool _paused;
    [SerializeField, Range(1, 20)] private int _stepsPerFrame = 2;

    [Header("Диагностика")]
    [SerializeField] private bool _debugMetrics = true;
    [SerializeField, Range(30, 600)] private int _debugEverySteps = 120;

    [Header("Проверяемые метрики (раздел диссертации)")]
    [Tooltip("Тест сохранения объёма α: засеять гауссову каплю, выключить источник/сток/тепло/плавучесть — поле только переносится. ∫α должно сохраняться; в журнал печатается дрейф % (консервативность схемы переноса).")]
    [SerializeField] private bool _conservationTest;
    [Tooltip("Сбрасывать снимки поля скорости в файлы VelocityDumps/ — для расчёта энергетического спектра E(k) при разных ℏ скриптом spectrum.py.")]
    [SerializeField] private bool _dumpVelocity;
    [SerializeField, Range(60, 3000)] private int _dumpEverySteps = 300;

    private CSISF _isf;
    private CSVelocity _vel;
    private CSScalarField _alpha;   // дым (плотность/оптика)
    private CSScalarField _temp;    // температура T (нагрев у очага, остывание) → плавучесть
    private CSVorticityConfine _vc; // vorticity confinement на транспортной скорости
    private ComputeBuffer _solidMask, _sourceMask, _sinkMask;
    private ComputeBuffer _inflowMask;      // широкое −X сечение: задаёт скорость впуска (Дирихле)
    private ComputeBuffer _visualSolidMask; // только то, что рисуем (мебель; внешние стены опц.)
    private ComputeBuffer _buoyB;           // потенциал плавучести (вертикальный интеграл α)
    private bool _initialized;
    // Кэш конфигурации масок — для перестроения в рантайме при изменении.
    private bool _lastOpenFlowFaces, _lastInflowFullFace;
    private float _lastInflowDepth;
    private double _alphaBaseline = -1; // базовый ∫α для теста сохранения объёма

    // Параметры струи входа (бегущая фаза, как в сценарии Jet).
    private float _kInX, _kInY, _kInZ, _omega;

    // Трассеры (A/B), полностью независимы от α.
    private CSParticles _particles;
    private ParticleGpuBuffers _particleBuffers;
    private Vector3[] _renderPos, _renderVel, _prevPos, _displayVelSmooth;
    private float[] _pxArr, _pyArr, _pzArr;
    private int _particlesCount, _maxParticles;
    private const float DisplayVelBlend = 0.32f;

    private void Start()
    {
        _isf = new CSISF();
        _isf.Init(_kernelsShader, _fftShader, _lesShader, vol_size, vol_res, hbar, dt);
        _vel = new CSVelocity(_isf.resX, _isf.resY, _isf.resZ);
        _alpha = new CSScalarField(_scalarShader, _isf.resX, _isf.resY, _isf.resZ, _isf.dx, _isf.dy, _isf.dz);
        _temp = new CSScalarField(_scalarShader, _isf.resX, _isf.resY, _isf.resZ, _isf.dx, _isf.dy, _isf.dz);
        _vc = new CSVorticityConfine(_scalarShader, _isf.resX, _isf.resY, _isf.resZ, _isf.dx, _isf.dy, _isf.dz);
        _buoyB = new ComputeBuffer(_isf.num, sizeof(float));

        _kInX = _inletVelocity.x / hbar;
        _kInY = _inletVelocity.y / hbar;
        _kInZ = _inletVelocity.z / hbar;
        _omega = _inletVelocity.sqrMagnitude / (2f * hbar);

        BuildMasks();
        _lastOpenFlowFaces = _openFlowFaces;
        _lastInflowFullFace = _inflowFullFace;
        _lastInflowDepth = _inflowDepth;
        InitPsiPlaneWave(_inletVelocity * 0.2f);
        RunInitBoundary(8);

        InitTracers();
        BuildWallMeshes();
        if (_conservationTest) SeedAlphaBlob();

        _initialized = true;
    }

    /// <summary>Гауссова капля α в центре домена — начальное условие для теста сохранения объёма (чистый перенос).</summary>
    private void SeedAlphaBlob()
    {
        int n = _isf.num;
        var b = new float[n];
        float cx = vol_size[0] * 0.5f, cy = vol_size[1] * 0.5f, cz = vol_size[2] * 0.5f;
        float sig = Mathf.Min(vol_size[0], Mathf.Min(vol_size[1], vol_size[2])) * 0.12f;
        float inv2s2 = 1f / (2f * sig * sig);
        for (int i = 0; i < n; i++)
        {
            float ex = _isf.pxCPU[i] - cx, ey = _isf.pyCPU[i] - cy, ez = _isf.pzCPU[i] - cz;
            b[i] = Mathf.Exp(-(ex * ex + ey * ey + ez * ez) * inv2s2);
        }
        _alpha.Alpha.SetData(b);
    }

    /// <summary>Полупрозрачные боксы перегородок/колонн как реальная геометрия — видны в Billboard/Shaded.</summary>
    private void BuildWallMeshes()
    {
        if (!_showWallMeshes) return;

        var temp = GameObject.CreatePrimitive(PrimitiveType.Cube);
        Mesh cube = temp.GetComponent<MeshFilter>().sharedMesh;
        Destroy(temp);

        var sh = Shader.Find("Universal Render Pipeline/Unlit");
        if (sh == null) sh = Shader.Find("Sprites/Default");
        var mat = new Material(sh) { color = _wallMeshColor };
        mat.SetColor("_BaseColor", _wallMeshColor);
        mat.SetFloat("_Surface", 1f);   // transparent
        mat.SetFloat("_ZWrite", 0f);
        mat.SetInt("_SrcBlend", (int)UnityEngine.Rendering.BlendMode.SrcAlpha);
        mat.SetInt("_DstBlend", (int)UnityEngine.Rendering.BlendMode.OneMinusSrcAlpha);
        mat.EnableKeyword("_SURFACE_TYPE_TRANSPARENT");
        mat.renderQueue = 3000;

        var parent = new GameObject("WallMeshes").transform;
        parent.SetParent(transform, false);

        if (wallSegments != null)
            foreach (var w in wallSegments) SpawnWallMesh(parent, cube, mat, w.center, w.half);
        SpawnWallMesh(parent, cube, mat, _obstacleCenter, _obstacleHalf);
    }

    private static void SpawnWallMesh(Transform parent, Mesh cube, Material mat, Vector3 c, Vector3 h)
    {
        if (h.sqrMagnitude < 1e-6f) return;
        var g = new GameObject("wall");
        g.transform.SetParent(parent, false);
        g.transform.localPosition = c;
        g.transform.localScale = h * 2f;
        g.AddComponent<MeshFilter>().sharedMesh = cube;
        g.AddComponent<MeshRenderer>().sharedMaterial = mat;
    }

    private void InitTracers()
    {
        _particleBuffers = GetComponent<ParticleGpuBuffers>();
        if (!_spawnTracers || _particlesShader == null) return;

        // Буфер = бюджет популяции + небольшой запас на приток до переработки.
        _maxParticles = Mathf.Clamp(_tracerPopulationCap + _tracerPerStep * 8, 256, 600000);
        _particles = new CSParticles();
        _particles.Init(_particlesShader, _maxParticles, _isf);
        _particleBuffers?.EnsureCapacity(_maxParticles);
        _renderPos = new Vector3[_maxParticles];
        _renderVel = new Vector3[_maxParticles];
        _prevPos = new Vector3[_maxParticles];
        _displayVelSmooth = new Vector3[_maxParticles];
        _pxArr = new float[_maxParticles];
        _pyArr = new float[_maxParticles];
        _pzArr = new float[_maxParticles];
    }

    private void Update()
    {
        if (!_initialized) return;
        RebuildMasksIfNeeded();
        if (_paused) return;
        for (int s = 0; s < _stepsPerFrame; s++)
        {
            iterator++;
            Step();
        }
        UpdateTracerDisplay();
        if (_dumpVelocity && _dumpEverySteps > 0 && iterator % _dumpEverySteps == 0)
            DumpVelocity();
    }

    private void OnDestroy()
    {
        ReleaseMasks();
        _buoyB?.Release();
        _particles?.Dispose();
        _alpha?.Dispose();
        _temp?.Dispose();
        _vc?.Dispose();
        _vel?.Dispose();
        _isf?.Dispose();
    }

    #region Step

    private void Step()
    {
        _isf.kinematicViscosity = _kinematicViscosity;
        _isf.UpdateSpace(_useLES, null);

        // Плавучесть (v3): фаза ψ из вертикального интеграла ТЕМПЕРАТУРЫ ⇒ подъём ∝ T. До проекции.
        if (_buoyancy != 0f && !_conservationTest)
        {
            _temp.ComputeBuoyancyPotential(_buoyB);
            _isf.ApplyPhaseField(_buoyB, _buoyancy * dt / hbar);
        }

        // Непрерывная струя: бегущая фаза -ω·t на входе (постоянная подкачка импульса, как сопло в Jet).
        float invH = 1f / hbar;
        float jetPhase = -_omega * dt * iterator;
        for (int b = 0; b < _boundaryIters; b++)
        {
            _isf.ApplyJetBoundary(_solidMask, 0f, 0f, 0f, 0f);
            _isf.ApplyJetBoundary(_inflowMask, _kInX, _kInY, _kInZ, jetPhase);
            if (_ventEnabled)
                _isf.ApplyJetBoundary(_sinkMask,
                    _ventVelocity.x * invH, _ventVelocity.y * invH, _ventVelocity.z * invH, 0f);
            _isf.PressureProject();
        }

        // Восстановление стабилизированной скорости и перенос полей той же скоростью (гл. 5.1.4–5.1.5).
        _isf.UpdateVelocities(_vel);
        // Vorticity confinement на транспортной ũ (клубящийся факел; ψ не трогаем).
        if (_vorticityConfinement > 0f)
            _vc.Apply(_vel, _vorticityConfinement, dt);
        // Температура: нагрев у очага (T=1), перенос, диффузия, остывание (decay). Сток вытяжки не охлаждает (factor=1).
        if (!_conservationTest)
            _temp.Step(_vel, _solidMask, _sourceMask, _sinkMask,
                dt, _thermalDiffusion, 1f, 1f, 1, _cooling * dt, _macCormack);
        // В тесте сохранения: источник/сток/затухание выключены — чистый перенос (проверка консервативности схемы).
        float aSource = _conservationTest ? 0f : _alphaSourceValue;
        float aSink = _conservationTest ? 1f : (_ventEnabled ? _alphaSinkFactor : 1f);
        _alpha.Step(_vel, _solidMask, _sourceMask, _sinkMask,
            dt, _alphaDiffusion, aSource, aSink, _alphaDiffuseIters, 0f, _macCormack);

        StepTracers();
        LogMetrics();
    }

    // Трассеры несутся ТОЙ ЖЕ скоростью _vel; на α не влияют (отдельная подсистема для A/B и Billboard/Shaded).
    private void StepTracers()
    {
        if (!_spawnTracers || _particles == null) return;
        SpawnTracers();
        _particles.CalculateMovement(_vel, clampSampling: true);
        if (_tracerDiffusion > 0f)
            _particles.DiffuseTracers(Mathf.Sqrt(2f * _tracerDiffusion * dt), iterator);
        _particles.PushOutOfSolidMask(_solidMask, 12, _ventCenter, 0f);
        _particles.ClampPositionsToVolume(vol_size[0], vol_size[1], vol_size[2]);
        CompactTracers();
    }

    private void SpawnTracers()
    {
        int n = _tracerPerStep;
        if (n <= 0) return;
        var xx = new float[n];
        var yy = new float[n];
        var zz = new float[n];
        for (int i = 0; i < n; i++)
        {
            xx[i] = Random.Range(_inletCenter.x - _inletHalf.x, _inletCenter.x + _inletHalf.x);
            yy[i] = Random.Range(_inletCenter.y - _inletHalf.y, _inletCenter.y + _inletHalf.y);
            zz[i] = Random.Range(_inletCenter.z - _inletHalf.z, _inletCenter.z + _inletHalf.z);
        }
        _particles.AddParticles(xx, yy, zz, n);
        _particlesCount = _particles.Size;
    }

    private void CompactTracers()
    {
        if (_particlesCount == 0) return;
        _particles.ReadPositions(_pxArr, _pyArr, _pzArr);
        Vector3 vmin = _ventCenter - _ventHalf, vmax = _ventCenter + _ventHalf;
        // Массив в порядке возраста (фронт — старейшие). При переполнении бюджета перерабатываем старейших.
        int over = _particlesCount - _tracerPopulationCap;
        for (int i = 0; i < _particlesCount; i++)
        {
            bool recycleOldest = over > 0 && i < over;
            bool inVent = _ventEnabled
                       && _pxArr[i] >= vmin.x && _pxArr[i] <= vmax.x
                       && _pyArr[i] >= vmin.y && _pyArr[i] <= vmax.y
                       && _pzArr[i] >= vmin.z && _pzArr[i] <= vmax.z;
            if (recycleOldest || inVent)
                _pxArr[i] = _pyArr[i] = _pzArr[i] = -1f;
        }
        _particles.WritePositions(_pxArr, _pyArr, _pzArr, _particlesCount);
        _particles.CompactParticles(_pxArr, _pyArr, _pzArr,
            vol_size[0], vol_size[1], vol_size[2], _prevPos, _displayVelSmooth);
        _particlesCount = _particles.Size;
    }

    private void UpdateTracerDisplay()
    {
        if (!_spawnTracers || _particles == null || _particlesCount == 0 || _particleBuffers == null) return;
        _particles.ReadPositions(_pxArr, _pyArr, _pzArr);
        var offset = transform.position;
        float maxX = vol_size[0], maxY = vol_size[1], maxZ = vol_size[2];
        float velThr = maxX * maxX + maxY * maxY + maxZ * maxZ;
        int visible = 0;
        for (int i = 0; i < _particlesCount; i++)
        {
            float px = _pxArr[i], py = _pyArr[i], pz = _pzArr[i];
            var pos = new Vector3(px, py, pz) + offset;
            var last = _prevPos[i];
            _prevPos[i] = pos;
            var vel = pos - last;
            if (vel.sqrMagnitude > velThr) vel = Vector3.zero;
            _displayVelSmooth[i] = Vector3.Lerp(_displayVelSmooth[i], vel, DisplayVelBlend);
            if (px < 0f || px > maxX || py < 0f || py > maxY || pz < 0f || pz > maxZ)
                continue;
            _renderPos[visible] = pos;
            _renderVel[visible] = _displayVelSmooth[i];
            visible++;
        }
        _particleBuffers.Upload(_renderPos, _renderVel, visible);
    }

    private void RunInitBoundary(int iterations)
    {
        for (int i = 0; i < iterations; i++)
        {
            _isf.ApplyJetBoundary(_solidMask, 0f, 0f, 0f, 0f);
            _isf.ApplyJetBoundary(_inflowMask, _kInX, _kInY, _kInZ, 0f);
            _isf.PressureProject();
        }
    }

    #endregion

    #region Init / Masks

    private void InitPsiPlaneWave(Vector3 backgroundU)
    {
        int num = _isf.num;
        float kx = backgroundU.x / hbar, ky = backgroundU.y / hbar, kz = backgroundU.z / hbar;
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

    private static bool InAabb(float px, float py, float pz, Vector3 c, Vector3 h)
    {
        return px >= c.x - h.x && px <= c.x + h.x
            && py >= c.y - h.y && py <= c.y + h.y
            && pz >= c.z - h.z && pz <= c.z + h.z;
    }

    private bool IsOuterWall(float px, float py, float pz)
    {
        float m = _wallMargin;
        // Пол/потолок (Y) и боковые стены (Z) — всегда. Грани ±X (продувка) — опционально открыты.
        bool yz = py < m || py > vol_size[1] - m
               || pz < m || pz > vol_size[2] - m;
        if (_openFlowFaces)
            return yz;
        bool x = px < m || px > vol_size[0] - m;
        return x || yz;
    }

    private bool IsVent(float px, float py, float pz)
    {
        if (InAabb(px, py, pz, _ventCenter, _ventHalf)) return true;
        if (extraVents != null)
            for (int k = 0; k < extraVents.Length; k++)
                if (InAabb(px, py, pz, extraVents[k].center, extraVents[k].half)) return true;
        return false;
    }

    private bool IsSolidAt(float px, float py, float pz)
    {
        // Источник и вытяжки прорезают отверстия (не твёрдые).
        if (InAabb(px, py, pz, _inletCenter, _inletHalf)) return false;
        if (IsVent(px, py, pz)) return false;
        if (IsOuterWall(px, py, pz)) return true;
        if (_obstacleHalf.sqrMagnitude > 1e-6f && InAabb(px, py, pz, _obstacleCenter, _obstacleHalf)) return true;
        if (wallSegments != null)
            for (int k = 0; k < wallSegments.Length; k++)
                if (InAabb(px, py, pz, wallSegments[k].center, wallSegments[k].half)) return true;
        return false;
    }

    /// <summary>Перестроить маски на лету при изменении конфигурации стен/впуска (рантайм-щелчки в Play).</summary>
    private void RebuildMasksIfNeeded()
    {
        if (_openFlowFaces == _lastOpenFlowFaces
            && _inflowFullFace == _lastInflowFullFace
            && Mathf.Approximately(_inflowDepth, _lastInflowDepth))
            return;

        ReleaseMasks();
        BuildMasks();
        _lastOpenFlowFaces = _openFlowFaces;
        _lastInflowFullFace = _inflowFullFace;
        _lastInflowDepth = _inflowDepth;
    }

    private void ReleaseMasks()
    {
        _solidMask?.Release();
        _visualSolidMask?.Release();
        _inflowMask?.Release();
        _sourceMask?.Release();
        _sinkMask?.Release();
    }

    /// <summary>Единые маски (стены/источник/сток): один источник правды для ψ-границы, переноса α и оптики.</summary>
    private void BuildMasks()
    {
        int num = _isf.num;
        var solid = new int[num];
        var visual = new int[num];
        var src = new int[num];
        var sink = new int[num];
        var inflow = new int[num];
        for (int i = 0; i < num; i++)
        {
            float px = _isf.pxCPU[i], py = _isf.pyCPU[i], pz = _isf.pzCPU[i];
            bool inSrc = InAabb(px, py, pz, _inletCenter, _inletHalf);
            bool inSink = IsVent(px, py, pz);
            bool isSolid = !inSrc && !inSink && IsSolidAt(px, py, pz);
            solid[i] = isSolid ? 1 : 0;
            // Визуальная маска: по умолчанию только внутренние препятствия (без внешних стен), иначе ничего не видно.
            visual[i] = isSolid && (_renderOuterWalls || !IsOuterWall(px, py, pz)) ? 1 : 0;
            src[i] = inSrc ? 1 : 0;
            sink[i] = inSink ? 1 : 0;
            // Впуск: либо всё входное сечение (полный сквозняк, авто-масштаб), либо локализованная струя.
            bool inInflow = _inflowFullFace
                ? (px <= _inflowDepth)
                : InAabb(px, py, pz, _inflowCenter, _inflowHalf);
            inflow[i] = (inInflow && !isSolid && !inSink) ? 1 : 0;
        }
        _solidMask = new ComputeBuffer(num, sizeof(int)); _solidMask.SetData(solid);
        _visualSolidMask = new ComputeBuffer(num, sizeof(int)); _visualSolidMask.SetData(visual);
        _sourceMask = new ComputeBuffer(num, sizeof(int)); _sourceMask.SetData(src);
        _sinkMask = new ComputeBuffer(num, sizeof(int)); _sinkMask.SetData(sink);
        _inflowMask = new ComputeBuffer(num, sizeof(int)); _inflowMask.SetData(inflow);
    }

    #endregion

    #region Diagnostics

    private void LogMetrics()
    {
        if (!_debugMetrics || _debugEverySteps <= 0 || iterator % _debugEverySteps != 0)
            return;
        int n = _isf.num;
        var a = new float[n];
        _alpha.Alpha.GetData(a);
        double sum = 0; float mx = 0; int occupied = 0;
        for (int i = 0; i < a.Length; i++)
        {
            sum += a[i];
            if (a[i] > mx) mx = a[i];
            if (a[i] > 0.05f) occupied++;
        }

        // Скорость + проверяемые метрики: кин. энергия, ошибка представимости (RMS дивергенции восстановленной u).
        var vx = new float[n]; var vy = new float[n]; var vz = new float[n];
        _vel.vx.GetData(vx); _vel.vy.GetData(vy); _vel.vz.GetData(vz);
        int rx = _isf.resX, ry = _isf.resY, rz = _isf.resZ;
        float dx = _isf.dx, dy = _isf.dy, dz = _isf.dz;
        double uSum = 0, uSq = 0, divSq = 0; float uMax = 0;
        for (int i = 0; i < rx; i++)
        for (int j = 0; j < ry; j++)
        for (int k = 0; k < rz; k++)
        {
            int idx = i * ry * rz + j * rz + k;
            float ux = vx[idx], uy = vy[idx], uz = vz[idx];
            float m = Mathf.Sqrt(ux * ux + uy * uy + uz * uz);
            uSum += m; uSq += ux * ux + uy * uy + uz * uz; if (m > uMax) uMax = m;
            // Дискретная дивергенция (обратные разности, как в Div-ядре ISF) — невязка несжимаемости.
            int im = i > 0 ? idx - ry * rz : idx;
            int jm = j > 0 ? idx - rz : idx;
            int km = k > 0 ? idx - 1 : idx;
            double div = (ux - vx[im]) / dx + (uy - vy[jm]) / dy + (uz - vz[km]) / dz;
            divSq += div * div;
        }
        double uRMS = Math.Sqrt(uSq / n);
        double divRMS = Math.Sqrt(divSq / n);
        double cell = (dx + dy + dz) / 3.0;
        double repErr = uRMS > 1e-9 ? divRMS * cell / uRMS : 0;   // безразмерная: относит. дивергенция на ячейку
        double Ekin = 0.5 * uSq * (dx * dy * dz);                 // полная кин. энергия поля

        string consv = "";
        if (_conservationTest)
        {
            if (_alphaBaseline < 0) _alphaBaseline = sum;
            double drift = _alphaBaseline > 1e-9 ? (sum - _alphaBaseline) / _alphaBaseline * 100.0 : 0;
            consv = $" | CONS ∫α={sum:F2} base={_alphaBaseline:F2} drift={drift:+0.000;-0.000}%";
        }

        Debug.Log($"[Hybrid3D] step={iterator} hbar={hbar:0.###} alphaSum={sum:F1} alphaMax={mx:F3} " +
                  $"occupied={occupied}/{n}({100.0 * occupied / n:F1}%) |u|mean={uSum / n:F3} |u|max={uMax:F3} " +
                  $"Ekin={Ekin:F3} repErr={repErr:E2} divRMS={divRMS:E2} particles={_particlesCount}/{_tracerPopulationCap}{consv}");
    }

    /// <summary>Снимок поля скорости в бинарный файл VelocityDumps/ для оффлайн-расчёта спектра E(k) (spectrum.py).</summary>
    private void DumpVelocity()
    {
        int n = _isf.num;
        var vx = new float[n]; var vy = new float[n]; var vz = new float[n];
        _vel.vx.GetData(vx); _vel.vy.GetData(vy); _vel.vz.GetData(vz);
        string dir = Path.Combine(Application.dataPath, "..", "VelocityDumps");
        Directory.CreateDirectory(dir);
        string path = Path.Combine(dir, $"vel_hbar{hbar:0.####}_res{_isf.resX}x{_isf.resY}x{_isf.resZ}_step{iterator}.bin");
        using (var w = new BinaryWriter(File.Open(path, FileMode.Create)))
        {
            w.Write(_isf.resX); w.Write(_isf.resY); w.Write(_isf.resZ);
            w.Write(_isf.dx); w.Write(_isf.dy); w.Write(_isf.dz); w.Write(hbar);
            var bytes = new byte[n * 4];
            Buffer.BlockCopy(vx, 0, bytes, 0, n * 4); w.Write(bytes);
            Buffer.BlockCopy(vy, 0, bytes, 0, n * 4); w.Write(bytes);
            Buffer.BlockCopy(vz, 0, bytes, 0, n * 4); w.Write(bytes);
        }
        Debug.Log($"[Hybrid3D] velocity dump → {path}");
    }

    #endregion

    #region Raymarch sources

    public bool TryGetScalarField(out ComputeBuffer alpha,
        out Vector3 volumeMinWorld, out Vector3 volumeSizeWorld,
        out int resX, out int resY, out int resZ)
    {
        alpha = null; volumeMinWorld = default; volumeSizeWorld = default; resX = resY = resZ = 0;
        if (!_initialized || _alpha == null) return false;
        alpha = _alpha.Alpha;
        volumeMinWorld = transform.position;
        volumeSizeWorld = new Vector3(vol_size[0], vol_size[1], vol_size[2]);
        resX = _isf.resX; resY = _isf.resY; resZ = _isf.resZ;
        return true;
    }

    public bool TryGetSolidMask(out ComputeBuffer solidMask, out int resX, out int resY, out int resZ)
    {
        // Для рендера отдаём ВИЗУАЛЬНУЮ маску (без внешних стен) — физика потока по-прежнему на полной _solidMask.
        solidMask = null; resX = resY = resZ = 0;
        if (!_initialized || _visualSolidMask == null) return false;
        solidMask = _visualSolidMask;
        resX = _isf.resX; resY = _isf.resY; resZ = _isf.resZ;
        return true;
    }

    #endregion

    #region Gizmos

    private void OnDrawGizmos()
    {
        var o = transform.position;
        var vol = new Vector3(vol_size[0], vol_size[1], vol_size[2]);
        Gizmos.color = Color.yellow;
        Gizmos.DrawWireCube(o + vol / 2f, vol);
        Gizmos.color = new Color(0.9f, 0.9f, 0.2f, 0.5f);
        Gizmos.DrawWireCube(o + _inflowCenter, _inflowHalf * 2f);
        Gizmos.color = new Color(1f, 0.55f, 0.1f, 0.9f);
        Gizmos.DrawWireCube(o + _inletCenter, _inletHalf * 2f);
        Gizmos.color = new Color(0.2f, 0.85f, 0.35f, 0.9f);
        Gizmos.DrawWireCube(o + _ventCenter, _ventHalf * 2f);
        if (_obstacleHalf.sqrMagnitude > 1e-6f)
        {
            Gizmos.color = new Color(0.4f, 0.5f, 0.7f, 0.9f);
            Gizmos.DrawWireCube(o + _obstacleCenter, _obstacleHalf * 2f);
        }
        Gizmos.color = new Color(0.55f, 0.6f, 0.7f, 0.9f);
        if (wallSegments != null)
            foreach (var w in wallSegments)
                Gizmos.DrawWireCube(o + w.center, w.half * 2f);
        Gizmos.color = new Color(0.2f, 0.85f, 0.35f, 0.9f);
        if (extraVents != null)
            foreach (var v in extraVents)
                Gizmos.DrawWireCube(o + v.center, v.half * 2f);
    }

    #endregion
}
