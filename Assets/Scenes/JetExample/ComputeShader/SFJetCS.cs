using UnityEngine;
using ComputeShaderSF;
using ShrodingerFlow.Particles;

public class SFJetCS : SFBase, ISimulationParticleSizeSource
{
    [Header("Compute Shaders")]
    [SerializeField] private ComputeShader _kernelsShader;
    [SerializeField] private ComputeShader _fftShader;
    [SerializeField] private ComputeShader _particlesShader;
    [SerializeField] private ComputeShader _lesShader;

    [Header("Начальные условия")]
    [SerializeField] private int[] vol_size = { 4, 2, 2 };
    [SerializeField] private int[] vol_res = { 128, 64, 64 };
    [SerializeField] private float hbar = 0.02f;
    [SerializeField] private float dt = 1f / 48f;

    [SerializeField] private Vector3 jetVelocity = new Vector3(1f, 0f, 0f);
    [SerializeField] private Vector3 nozzleCen = new Vector3(0.3f, 0.9656632f, 1.0659939f);
    [SerializeField] private float nozzleLen = 0.5f;
    [SerializeField] private float nozzleRad = 0.3f;
    [SerializeField] private int n_particles = 50;

    [Header("Дополнительные настройки")]
    [SerializeField] private float _particleSize = 0.05f;
    [SerializeField] private bool _useLES = false;

    [Header("Скорость симуляции")]
    [SerializeField] private bool _paused;
    [SerializeField, Range(1, 20)] private int _stepsPerFrame = 1;

    private CSISF _isf;
    private CSParticles _particles;
    private CSVelocity _vel;
    private ComputeBuffer _isJetBuf;

    private ParticleGpuBuffers _particleBuffers;
    private ParticleDisplay3D _particleDisplay;
    private Vector3[] _renderPos;
    private Vector3[] _renderVel;
    private float[] _pxArr, _pyArr, _pzArr;
    private Vector3[] _prevPos;
    private int _particlesCount;

    private float _kvecX, _kvecY, _kvecZ;
    private float _omega;
    private Vector3 _volSizeV3;
    private bool _initialized;
    private int _compactCounter;
    private float _particleSizeSyncedForDisplay = float.NaN;

    public float SimulationParticleSize => _particleSize;

    private void Start()
    {
        _volSizeV3 = new Vector3(vol_size[0], vol_size[1], vol_size[2]);
        _particleBuffers = GetComponent<ParticleGpuBuffers>();
        _particleDisplay = GetComponent<ParticleDisplay3D>();

        _isf = new CSISF();
        _isf.Init(_kernelsShader, _fftShader, _lesShader, vol_size, vol_res, hbar, dt);

        _particles = new CSParticles();
        _particles.Init(_particlesShader, n_particles * 1000, _isf);

        _vel = new CSVelocity(_isf.resX, _isf.resY, _isf.resZ);

        int maxCloud = n_particles * 1000;
        _particleBuffers?.EnsureCapacity(maxCloud);
        _renderPos = new Vector3[maxCloud];
        _renderVel = new Vector3[maxCloud];
        _pxArr = new float[maxCloud];
        _pyArr = new float[maxCloud];
        _pzArr = new float[maxCloud];
        _prevPos = new Vector3[maxCloud];

        InitPsi();
        BuildIsJetMask();
        InitJetFlow();

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

    private void InitPsi()
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

    private void BuildIsJetMask()
    {
        int num = _isf.num;
        var mask = new int[num];
        for (int i = 0; i < num; i++)
        {
            float px = _isf.pxCPU[i];
            float py = _isf.pyCPU[i];
            float pz = _isf.pzCPU[i];
            bool inJet = Mathf.Abs(px - nozzleCen.x) <= nozzleLen / 2f
                && (py - nozzleCen.y) * (py - nozzleCen.y)
                 + (pz - nozzleCen.z) * (pz - nozzleCen.z) <= nozzleRad * nozzleRad;
            mask[i] = inJet ? 1 : 0;
        }
        _isJetBuf = new ComputeBuffer(num, sizeof(int));
        _isJetBuf.SetData(mask);

        _kvecX = jetVelocity.x / hbar;
        _kvecY = jetVelocity.y / hbar;
        _kvecZ = jetVelocity.z / hbar;
        _omega = (jetVelocity.x * jetVelocity.x
                + jetVelocity.y * jetVelocity.y
                + jetVelocity.z * jetVelocity.z) / (2f * hbar);
    }

    private void InitJetFlow()
    {
        for (int iter = 0; iter < 10; iter++)
        {
            _isf.ApplyJetBoundary(_isJetBuf, _kvecX, _kvecY, _kvecZ, 0f);
            _isf.PressureProject();
        }
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

            _compactCounter++;
            if (_compactCounter >= 60)
            {
                _compactCounter = 0;
                _particles.CompactParticles(_pxArr, _pyArr, _pzArr,
                    vol_size[0], vol_size[1], vol_size[2]);
                _particlesCount = _particles.Size;
            }
        }

        UpdateParticleSystem();
    }

    private void SimulationStep()
    {
        _isf.UpdateSpace(_useLES);

        float phaseOffset = -_omega * dt * iterator;
        _isf.ApplyJetBoundary(_isJetBuf, _kvecX, _kvecY, _kvecZ, phaseOffset);
        _isf.PressureProject();

        SpawnParticles();

        _isf.UpdateVelocities(_vel);
        _particles.CalculateMovement(_vel);
    }

    private void SpawnParticles()
    {
        var xx = new float[n_particles];
        var yy = new float[n_particles];
        var zz = new float[n_particles];
        for (int i = 0; i < n_particles; i++)
        {
            float rt = Random.value * 2f * Mathf.PI;
            xx[i] = nozzleCen.x;
            yy[i] = nozzleCen.y + 0.9f * nozzleRad * Mathf.Cos(rt);
            zz[i] = nozzleCen.z + 0.9f * nozzleRad * Mathf.Sin(rt);
        }
        _particles.AddParticles(xx, yy, zz, n_particles);
        _particlesCount = _particles.Size;
    }

    private void UpdateParticleSystem()
    {
        if (_particlesCount == 0 || _particleBuffers == null) return;

        _particles.ReadPositions(_pxArr, _pyArr, _pzArr);

        var offset = transform.position;
        float maxX = vol_size[0], maxY = vol_size[1], maxZ = vol_size[2];
        int visible = 0;
        for (int i = 0; i < _particlesCount; i++)
        {
            float px = _pxArr[i], py = _pyArr[i], pz = _pzArr[i];
            if (px < 0f || px > maxX || py < 0f || py > maxY || pz < 0f || pz > maxZ)
                continue;

            var pos = new Vector3(px, py, pz) + offset;
            var vel = pos - _prevPos[i];
            _prevPos[i] = pos;

            _renderPos[visible] = pos;
            _renderVel[visible] = vel;
            visible++;
        }

        _particleBuffers.Upload(_renderPos, _renderVel, visible);
    }

    private void OnDestroy()
    {
        _isJetBuf?.Release();
        _vel?.Dispose();
        _particles?.Dispose();
        _isf?.Dispose();
    }

    [ExecuteAlways]
    private void OnDrawGizmos()
    {
        var volV3 = new Vector3(vol_size[0], vol_size[1], vol_size[2]);
        Gizmos.color = Color.yellow;
        Gizmos.DrawWireCube(transform.position + volV3 / 2f, volV3);
    }
}
