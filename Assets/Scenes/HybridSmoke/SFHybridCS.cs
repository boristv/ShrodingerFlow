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

    [Header("Трассеры — опциональный A/B-режим (для Billboard/Shaded)")]
    [Tooltip("Считать частицы-трассеры параллельно α. Несутся ТОЙ ЖЕ скоростью ISF+LES. Не влияют на α-поле. Видны в режимах Billboard/Shaded.")]
    [SerializeField] private bool _spawnTracers = true;
    [SerializeField, Range(0, 1000)] private int _tracerPerStep = 80;
    [Tooltip("Авто-срок жизни = время пересечения комнаты (длина / скорость впуска) × запас. Само масштабируется под размер сцены — вручную менять не нужно.")]
    [SerializeField] private bool _autoTracerLifetime = true;
    [SerializeField, Range(1f, 6f)] private float _tracerLifetimeSafety = 3f;
    [Tooltip("Ручной срок жизни (шагов) — используется, только если авто выключен.")]
    [SerializeField, Range(0, 4000)] private int _tracerLifetime = 600;
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

    private CSISF _isf;
    private CSVelocity _vel;
    private CSScalarField _alpha;
    private ComputeBuffer _solidMask, _sourceMask, _sinkMask;
    private ComputeBuffer _inflowMask;      // широкое −X сечение: задаёт скорость впуска (Дирихле)
    private ComputeBuffer _visualSolidMask; // только то, что рисуем (мебель; внешние стены опц.)
    private bool _initialized;

    // Параметры струи входа (бегущая фаза, как в сценарии Jet).
    private float _kInX, _kInY, _kInZ, _omega;

    // Трассеры (A/B), полностью независимы от α.
    private CSParticles _particles;
    private ParticleGpuBuffers _particleBuffers;
    private Vector3[] _renderPos, _renderVel, _prevPos, _displayVelSmooth;
    private float[] _pxArr, _pyArr, _pzArr, _age;
    private int _particlesCount, _maxParticles;
    private int _lifeSteps; // эффективный срок жизни трассера (авто из времени пересечения или ручной)
    private const float DisplayVelBlend = 0.32f;

    private void Start()
    {
        _isf = new CSISF();
        _isf.Init(_kernelsShader, _fftShader, _lesShader, vol_size, vol_res, hbar, dt);
        _vel = new CSVelocity(_isf.resX, _isf.resY, _isf.resZ);
        _alpha = new CSScalarField(_scalarShader, _isf.resX, _isf.resY, _isf.resZ, _isf.dx, _isf.dy, _isf.dz);

        _kInX = _inletVelocity.x / hbar;
        _kInY = _inletVelocity.y / hbar;
        _kInZ = _inletVelocity.z / hbar;
        _omega = _inletVelocity.sqrMagnitude / (2f * hbar);

        BuildMasks();
        InitPsiPlaneWave(_inletVelocity * 0.2f);
        RunInitBoundary(8);

        InitTracers();
        BuildWallMeshes();

        _initialized = true;
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

        // Срок жизни из физики: сколько шагов трассер летит через комнату при скорости впуска, × запас.
        float speed = Mathf.Max(0.05f, _inletVelocity.magnitude);
        int crossSteps = Mathf.CeilToInt((vol_size[0] / speed) / Mathf.Max(1e-4f, dt));
        _lifeSteps = _autoTracerLifetime
            ? Mathf.CeilToInt(crossSteps * _tracerLifetimeSafety)
            : _tracerLifetime;

        // Бюджет = приток за всё время жизни (+запас); масштабируется вместе со сроком жизни.
        _maxParticles = Mathf.Clamp(_tracerPerStep * (_lifeSteps + 8), 256, 400000);
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
        _age = new float[_maxParticles];
    }

    private void Update()
    {
        if (!_initialized || _paused) return;
        for (int s = 0; s < _stepsPerFrame; s++)
        {
            iterator++;
            Step();
        }
        UpdateTracerDisplay();
    }

    private void OnDestroy()
    {
        _solidMask?.Release();
        _visualSolidMask?.Release();
        _inflowMask?.Release();
        _sourceMask?.Release();
        _sinkMask?.Release();
        _particles?.Dispose();
        _alpha?.Dispose();
        _vel?.Dispose();
        _isf?.Dispose();
    }

    #region Step

    private void Step()
    {
        _isf.kinematicViscosity = _kinematicViscosity;
        _isf.UpdateSpace(_useLES, null);

        // Непрерывная струя: бегущая фаза -ω·t на входе (постоянная подкачка импульса, как сопло в Jet).
        float invH = 1f / hbar;
        float jetPhase = -_omega * dt * iterator;
        for (int b = 0; b < _boundaryIters; b++)
        {
            _isf.ApplyJetBoundary(_solidMask, 0f, 0f, 0f, 0f);
            _isf.ApplyJetBoundary(_inflowMask, _kInX, _kInY, _kInZ, jetPhase);
            _isf.ApplyJetBoundary(_sinkMask,
                _ventVelocity.x * invH, _ventVelocity.y * invH, _ventVelocity.z * invH, 0f);
            _isf.PressureProject();
        }

        // Восстановление стабилизированной скорости и перенос α той же скоростью (гл. 5.1.4–5.1.5).
        _isf.UpdateVelocities(_vel);
        _alpha.Step(_vel, _solidMask, _sourceMask, _sinkMask,
            dt, _alphaDiffusion, _alphaSourceValue, _alphaSinkFactor, _alphaDiffuseIters);

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
        int before = _particles.Size;
        _particles.AddParticles(xx, yy, zz, n);
        _particlesCount = _particles.Size;
        for (int i = before; i < _particlesCount; i++)
            _age[i] = 0f;
    }

    private void CompactTracers()
    {
        if (_particlesCount == 0) return;
        _particles.ReadPositions(_pxArr, _pyArr, _pzArr);
        Vector3 vmin = _ventCenter - _ventHalf, vmax = _ventCenter + _ventHalf;
        for (int i = 0; i < _particlesCount; i++)
        {
            _age[i] += 1f;
            bool tooOld = _lifeSteps > 0 && _age[i] > _lifeSteps;
            bool inVent = _pxArr[i] >= vmin.x && _pxArr[i] <= vmax.x
                       && _pyArr[i] >= vmin.y && _pyArr[i] <= vmax.y
                       && _pzArr[i] >= vmin.z && _pzArr[i] <= vmax.z;
            if (tooOld || inVent)
                _pxArr[i] = _pyArr[i] = _pzArr[i] = -1f;
        }
        _particles.WritePositions(_pxArr, _pyArr, _pzArr, _particlesCount);
        _particles.CompactParticles(_pxArr, _pyArr, _pzArr,
            vol_size[0], vol_size[1], vol_size[2], _prevPos, _displayVelSmooth, _age);
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

        // Скорость: затухает поле после стартового плюма или держится? (ключевая диагностика)
        var vx = new float[n]; var vy = new float[n]; var vz = new float[n];
        _vel.vx.GetData(vx); _vel.vy.GetData(vy); _vel.vz.GetData(vz);
        double uSum = 0; float uMax = 0;
        for (int i = 0; i < n; i++)
        {
            float m = Mathf.Sqrt(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
            uSum += m; if (m > uMax) uMax = m;
        }

        Debug.Log($"[Hybrid3D] step={iterator} alphaSum={sum:F1} alphaMax={mx:F3} occupied(>0.05)={occupied}/{n} ({100.0 * occupied / n:F1}%) |u|mean={uSum / n:F3} |u|max={uMax:F3} lifeSteps={_lifeSteps} particles={_particlesCount}");
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
