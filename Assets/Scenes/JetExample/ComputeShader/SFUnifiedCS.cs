using UnityEngine;
using ComputeShaderSF;
using ShrodingerFlow.Particles;
#if UNITY_EDITOR
using UnityEditor;
#endif

// «Apply Scenario Defaults» берёт числа из SFUnifiedScenarioPresets (ассет или встроенная копия),
// а не из switch в этом файле. Референс: JetExampleCS.unity, UnifiedCS.unity — см. SFUnifiedScenarioPresets.cs.

public class SFUnifiedCS : SFBase, ISimulationParticleSizeSource
{
    public enum ScenarioType
    {
        Jet,
        SphereObstacle,
        CylinderObstacle,
        TwoSpheres,
        LeapfrogRings
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

    [Header("Частицы")]
    [SerializeField] private int _nParticles = 50;
    [SerializeField] private float _particleSize = 0.1f;

    [Header("Начальное расположение частиц (Box-спавн)")]
    [SerializeField] private Vector2 _boxSpawnX = new Vector2(0.3f, 0.3f);
    [SerializeField] private Vector2 _boxSpawnY = new Vector2(0.5f, 1.5f);
    [SerializeField] private Vector2 _boxSpawnZ = new Vector2(0.5f, 1.5f);

    [Header("Управление")]
    [SerializeField] private bool _useLES = false;
    [SerializeField] private bool _paused;
    [SerializeField, Range(1, 20)] private int _stepsPerFrame = 3;

    private CSISF _isf;
    private CSParticles _particles;
    private CSVelocity _vel;
    private ComputeBuffer _maskBuf1, _maskBuf2;

    private ParticleGpuBuffers _particleBuffers;
    private ParticleDisplay3D _particleDisplay;
    private Vector3[] _renderPos;
    private Vector3[] _renderVel;
    private float[] _pxArr, _pyArr, _pzArr;
    private Vector3[] _prevPos;
    private int _particlesCount;

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
                             || _scenario == ScenarioType.TwoSpheres;
        int maxParticles = oneTimeParticles ? _nParticles : _nParticles * 1000;

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

            if (_spawnEachStep && _scenario == ScenarioType.Jet)
            {
                _compactCounter++;
                if (_compactCounter >= 60)
                {
                    _compactCounter = 0;
                    _particles.CompactParticles(_pxArr, _pyArr, _pzArr,
                        vol_size[0], vol_size[1], vol_size[2]);
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
        _isf.UpdateSpace(_useLES);

        if (_boundaryEachStep)
        {
            float phaseOffset = (_scenario == ScenarioType.Jet)
                ? -_omega * dt * iterator
                : 0f;
            _isf.ApplyJetBoundary(_maskBuf1, _kvecX, _kvecY, _kvecZ, phaseOffset);
            _isf.PressureProject();
        }

        if (_spawnEachStep)
            SpawnNozzleParticles();

        _isf.UpdateVelocities(_vel);
        _particles.CalculateMovement(_vel);

        if (_scenario != ScenarioType.Jet)
            _particles.WrapPositions(vol_size[0], vol_size[1], vol_size[2]);
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

        bool ring = _scenario != ScenarioType.Jet;
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
        bool cull = _scenario == ScenarioType.Jet;
        float maxX = vol_size[0], maxY = vol_size[1], maxZ = vol_size[2];
        float velThreshold = maxX * maxX + maxY * maxY + maxZ * maxZ;
        int visible = 0;

        for (int i = 0; i < _particlesCount; i++)
        {
            float px = _pxArr[i], py = _pyArr[i], pz = _pzArr[i];

            if (cull && (px < 0f || px > maxX || py < 0f || py > maxY || pz < 0f || pz > maxZ))
                continue;

            var pos = new Vector3(px, py, pz) + offset;
            var lastPos = _prevPos[i];
            _prevPos[i] = pos;

            var vel = pos - lastPos;
            if (vel.sqrMagnitude > velThreshold)
                vel = Vector3.zero;

            _renderPos[visible] = pos;
            _renderVel[visible] = vel;
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
        }
    }

    #endregion
}
