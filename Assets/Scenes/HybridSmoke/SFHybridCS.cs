using UnityEngine;
using ComputeShaderSF;
using ShrodingerFlow.Particles;

/// <summary>
/// Гибридная волново-плотностная модель ISF (диссертация Тюлькина): ψ — вихревое ядро, α — поле состояния
/// (плотность дыма / интерфейс / оптика), переносимое скоростью ũ из ψ+LES. v1: ИЗОТЕРМИЧЕСКАЯ 3D-комната.
/// Шаг (гл. 5.1): ISF(ψ)→u→LES→проекция→advect(α)→границы. Силы вводятся через состояние (фазовый импульс на входе).
/// Единые маски (гл. 4.4.9): стены/источник/сток видимы и ψ-границе, и переносу α, и оптике.
/// </summary>
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
    [Tooltip("Рисовать внешние стены в раймарче. Обычно выкл: они закрывают обзор. Границей потока остаются в любом случае.")]
    [SerializeField] private bool _renderOuterWalls;

    [Header("Поток / ветер (объёмная продувка через состояние)")]
    [Tooltip("Однородная сила-«ветер» (фаза на обе компоненты ψ): создаёт тягу вход→вытяжка.")]
    [SerializeField] private Vector3 _wind = new Vector3(0.4f, 0f, 0f);
    [SerializeField, Range(1, 8)] private int _boundaryIters = 4;
    [SerializeField, Range(0, 240)] private int _rampSteps = 60;

    [Header("Плотностно-фазовое поле α")]
    [Tooltip("Целевое α в источнике (концентрация дыма на входе).")]
    [SerializeField, Range(0f, 1f)] private float _alphaSourceValue = 1f;
    [Tooltip("Множитель α в вытяжке за шаг (<1 — откачка).")]
    [SerializeField, Range(0f, 1f)] private float _alphaSinkFactor = 0.6f;
    [Tooltip("Регуляризующая диффузия D_α (мягкость интерфейса, заполнение застойных зон).")]
    [SerializeField] private float _alphaDiffusion = 0.0015f;
    [SerializeField, Range(0, 4)] private int _alphaDiffuseIters = 1;

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
    private ComputeBuffer _visualSolidMask; // только то, что рисуем (мебель; внешние стены опц.)
    private bool _initialized;

    private void Start()
    {
        _isf = new CSISF();
        _isf.Init(_kernelsShader, _fftShader, _lesShader, vol_size, vol_res, hbar, dt);
        _vel = new CSVelocity(_isf.resX, _isf.resY, _isf.resZ);
        _alpha = new CSScalarField(_scalarShader, _isf.resX, _isf.resY, _isf.resZ, _isf.dx, _isf.dy, _isf.dz);

        BuildMasks();
        InitPsiPlaneWave(_wind * 0.5f);
        RunInitBoundary(8);

        _initialized = true;
    }

    private void Update()
    {
        if (!_initialized || _paused) return;
        for (int s = 0; s < _stepsPerFrame; s++)
        {
            iterator++;
            Step();
        }
    }

    private void OnDestroy()
    {
        _solidMask?.Release();
        _visualSolidMask?.Release();
        _sourceMask?.Release();
        _sinkMask?.Release();
        _alpha?.Dispose();
        _vel?.Dispose();
        _isf?.Dispose();
    }

    #region Step

    private void Step()
    {
        _isf.kinematicViscosity = _kinematicViscosity;
        _isf.UpdateSpace(_useLES, null);

        float ramp = _rampSteps > 0 ? Mathf.Clamp01((float)iterator / _rampSteps) : 1f;
        ramp = ramp * ramp * (3f - 2f * ramp);

        if (_wind.sqrMagnitude > 1e-8f)
            _isf.ApplyUniformForce(_wind * ramp);

        float invH = ramp / hbar;
        for (int b = 0; b < _boundaryIters; b++)
        {
            _isf.ApplyJetBoundary(_solidMask, 0f, 0f, 0f, 0f);
            _isf.ApplyJetBoundary(_sourceMask,
                _inletVelocity.x * invH, _inletVelocity.y * invH, _inletVelocity.z * invH, 0f);
            _isf.ApplyJetBoundary(_sinkMask,
                _ventVelocity.x * invH, _ventVelocity.y * invH, _ventVelocity.z * invH, 0f);
            _isf.PressureProject();
        }

        // Восстановление стабилизированной скорости и перенос α той же скоростью (гл. 5.1.4–5.1.5).
        _isf.UpdateVelocities(_vel);
        _alpha.Step(_vel, _solidMask, _sourceMask, _sinkMask,
            dt, _alphaDiffusion, _alphaSourceValue, _alphaSinkFactor, _alphaDiffuseIters);

        LogMetrics();
    }

    private void RunInitBoundary(int iterations)
    {
        float invH = 1f / hbar;
        for (int i = 0; i < iterations; i++)
        {
            _isf.ApplyJetBoundary(_solidMask, 0f, 0f, 0f, 0f);
            _isf.ApplyJetBoundary(_sourceMask,
                _inletVelocity.x * invH, _inletVelocity.y * invH, _inletVelocity.z * invH, 0f);
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
        return px < m || px > vol_size[0] - m
            || py < m || py > vol_size[1] - m
            || pz < m || pz > vol_size[2] - m;
    }

    private bool IsSolidAt(float px, float py, float pz)
    {
        if (InAabb(px, py, pz, _inletCenter, _inletHalf)) return false;
        if (InAabb(px, py, pz, _ventCenter, _ventHalf)) return false;
        bool obstacle = _obstacleHalf.sqrMagnitude > 1e-6f && InAabb(px, py, pz, _obstacleCenter, _obstacleHalf);
        return IsOuterWall(px, py, pz) || obstacle;
    }

    /// <summary>Единые маски (стены/источник/сток): один источник правды для ψ-границы, переноса α и оптики.</summary>
    private void BuildMasks()
    {
        int num = _isf.num;
        var solid = new int[num];
        var visual = new int[num];
        var src = new int[num];
        var sink = new int[num];
        for (int i = 0; i < num; i++)
        {
            float px = _isf.pxCPU[i], py = _isf.pyCPU[i], pz = _isf.pzCPU[i];
            bool inSrc = InAabb(px, py, pz, _inletCenter, _inletHalf);
            bool inSink = InAabb(px, py, pz, _ventCenter, _ventHalf);
            bool isSolid = !inSrc && !inSink && IsSolidAt(px, py, pz);
            solid[i] = isSolid ? 1 : 0;
            // Визуальная маска: по умолчанию только внутренние препятствия (без внешних стен), иначе ничего не видно.
            visual[i] = isSolid && (_renderOuterWalls || !IsOuterWall(px, py, pz)) ? 1 : 0;
            src[i] = inSrc ? 1 : 0;
            sink[i] = inSink ? 1 : 0;
        }
        _solidMask = new ComputeBuffer(num, sizeof(int)); _solidMask.SetData(solid);
        _visualSolidMask = new ComputeBuffer(num, sizeof(int)); _visualSolidMask.SetData(visual);
        _sourceMask = new ComputeBuffer(num, sizeof(int)); _sourceMask.SetData(src);
        _sinkMask = new ComputeBuffer(num, sizeof(int)); _sinkMask.SetData(sink);
    }

    #endregion

    #region Diagnostics

    private void LogMetrics()
    {
        if (!_debugMetrics || _debugEverySteps <= 0 || iterator % _debugEverySteps != 0)
            return;
        var a = new float[_isf.num];
        _alpha.Alpha.GetData(a);
        double sum = 0; float mx = 0; int occupied = 0;
        for (int i = 0; i < a.Length; i++)
        {
            sum += a[i];
            if (a[i] > mx) mx = a[i];
            if (a[i] > 0.05f) occupied++;
        }
        Debug.Log($"[Hybrid3D] step={iterator} alphaSum={sum:F1} alphaMax={mx:F3} occupied(>0.05)={occupied}/{a.Length} ({100.0 * occupied / a.Length:F1}%)");
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
        Gizmos.color = new Color(1f, 0.55f, 0.1f, 0.9f);
        Gizmos.DrawWireCube(o + _inletCenter, _inletHalf * 2f);
        Gizmos.color = new Color(0.2f, 0.85f, 0.35f, 0.9f);
        Gizmos.DrawWireCube(o + _ventCenter, _ventHalf * 2f);
        if (_obstacleHalf.sqrMagnitude > 1e-6f)
        {
            Gizmos.color = new Color(0.4f, 0.5f, 0.7f, 0.9f);
            Gizmos.DrawWireCube(o + _obstacleCenter, _obstacleHalf * 2f);
        }
    }

    #endregion
}
