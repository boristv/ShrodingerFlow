using System.Numerics;
using ManagedCuda.VectorTypes;
using UnityEngine;
using source.assets.Discrete_space;
using source.assets.Particles;
using source;
using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;
using float3 = Unity.Mathematics.float3;
using Matrix4x4 = UnityEngine.Matrix4x4;
using Quaternion = UnityEngine.Quaternion;
using Vector2 = UnityEngine.Vector2;
using Vector3 = UnityEngine.Vector3;

public class SFJobBurst : MonoBehaviour
{
    float3[] p;
    ParticleSystem.Particle[] cloud;
    
    bool bPointsUpdated = false;
    //PARAMETERS
    [Header("Начальные условия")]
    [SerializeField, Tooltip("Box size")] private int[] vol_size = { 10, 5, 5 };      // box size
    [SerializeField, Tooltip("Volume resolution")] private int[] vol_res = { 128, 64, 64 };    // volume resolution
    [SerializeField, Tooltip("Planck constant")] private float hbar = (float)0.1;           // Planck constant
    [SerializeField, Tooltip("Time step")] private float dt = 1 / (float)12;          // time step
    int tmax = 85;
    [SerializeField, Tooltip("Background velocity")] private float[] background_vel = {(float) -0.2, 0, 0};
        
    
    [SerializeField, Tooltip("First ring radius")] private float r1 = (float)1.5; 
    [SerializeField, Tooltip("Second ring radius")] private float r2 = (float)0.9;              
    [SerializeField] private float[] n1 = {-1,0,0};         
    [SerializeField] private float[] n2 = {-1,0,0};

    [SerializeField, Tooltip("Particles count")] private int n_particles = 100000;

    Velocity vel;
    object lk= new object();

    UnityThreading.ActionThread myThread;

    [Header("Дополнительные настройки")] 
    [SerializeField, Tooltip("Размер частиц")] private float _particleSize = 0.1f;
    
    [Header("Начальное расположение частиц")] 
    [SerializeField] private Vector2 _boxSizeX = new Vector2(5, 5);
    [SerializeField] private Vector2 _boxSizeY = new Vector2(0.5f, 4.5f);
    [SerializeField] private Vector2 _boxSizeZ = new Vector2(0.5f, 4.5f);

    [SerializeField] private Mesh _mesh;
    [SerializeField] private Material _material;

    private ObjectPositionJob _job;
    private NativeArray<float3> _nativeP;
    private NativeArray<float> _nativeCubeYOffsets;
    private NativeArray<Matrix4x4> _nativeMatrices;
    private NativeArray<float3> _nativePositions;
    private RenderParams _rp;

    private bool _isReady;
    
    private void Start()
    {
        p = new float3[n_particles];
        
        _rp = new RenderParams(_material);
        
        var count = n_particles;
        _nativePositions = new NativeArray<float3>(count, Allocator.Persistent);
        _nativeMatrices = new NativeArray<Matrix4x4>(count, Allocator.Persistent);
        _nativeCubeYOffsets = new NativeArray<float>(count, Allocator.Persistent);
        _nativeP = new NativeArray<float3>(count, Allocator.Persistent);
        
        _job = new ObjectPositionJob
        {
            Positions = _nativePositions,
            Size = _particleSize * Vector3.one,
            Rotation = Quaternion.identity
        };
        
        myThread = UnityThreadHelper.CreateThread(Worker);
    }

    private void Worker()
    {
        float[] cen1 = { vol_size[0] / 2f, vol_size[1] / 2f, vol_size[2] / 2f };
        float[] cen2 = { vol_size[0] / 2f, vol_size[1] / 2f, vol_size[2] / 2f };

        //INITIALISATION
        ISF.Init(vol_size, vol_res, hbar, dt);
        Particles.init(n_particles);

        //init psi
        float[] kvec = { background_vel[0] / hbar, background_vel[1] / hbar, background_vel[2] / hbar };
        float phase;
        var tmp1 = new cuFloatComplex[ISF.properties.resx, ISF.properties.resy, ISF.properties.resz];
        var tmp2 = new cuFloatComplex[ISF.properties.resx, ISF.properties.resy, ISF.properties.resz];
        Complex tmp;
        for (int i = 0; i < vol_res[0]; i++)
        {
            for (int j = 0; j < vol_res[1]; j++)
            {
                for (int k = 0; k < vol_res[2]; k++)
                {
                    phase = kvec[0] * ISF.properties.px[i, j, k] +
                            kvec[1] * ISF.properties.py[i, j, k] +
                            kvec[2] * ISF.properties.pz[i, j, k];
                    tmp = Complex.Exp(Complex.ImaginaryOne * phase);
                    tmp1[i, j, k] = new cuFloatComplex((float)tmp.Real, (float)tmp.Imaginary);
                    tmp2[i, j, k] = new cuFloatComplex((float)(tmp.Real * 0.01), (float)(tmp.Imaginary * 0.01));
                }
            }
        }
        float d = ISF.properties.dx * 5;

        ISF.add_circle(tmp1, cen1, n1, r1, d);
        ISF.add_circle(tmp1, cen2, n2, r2, d);

        ISF.psi1.CopyToDevice(tmp1);
        ISF.psi2.CopyToDevice(tmp2);

        ISF.Normalize();
        ISF.PressureProject();

        //init particles
        var x = new float[n_particles];
        var y = new float[n_particles];
        var z = new float[n_particles];
        System.Random rnd = new System.Random();
        for (int i = 0; i < n_particles; i++)
        {
            x[i] = (float) (rnd.NextDouble() * (_boxSizeX.y - _boxSizeX.x) + _boxSizeX.x);
            y[i] = (float) (rnd.NextDouble() * (_boxSizeY.y - _boxSizeY.x) + _boxSizeY.x);
            z[i] = (float) (rnd.NextDouble() * (_boxSizeZ.y - _boxSizeZ.x) + _boxSizeZ.x);
        }

        Particles.add_particles(x, y, z, n_particles);

        vel = new Velocity(ISF.properties.resx, ISF.properties.resy, ISF.properties.resz);

        while (true) {
            //MAIN ITERATION
            //incompressible Schroedinger flow
            ISF.update_space();

            //particle update
            ISF.update_velocities(vel);

            Particles.calculate_movement(vel);

            float[] px = Particles.x;
            float[] py = Particles.y;
            float[] pz = Particles.z;

            for (int ii = 0; ii < n_particles; ++ii)
            {
                p[ii] = new Vector3(px[ii], py[ii], pz[ii]);
            }

            lock (lk)
            {
                _isReady = true;
            }

            if (UnityThreading.ActionThread.CurrentThread.ShouldStop)
                return; // finish
        }
    }
    
    public void Update()
    {
        if (_isReady)
        {
            _job.Matrices = _nativeMatrices;
            for (int ii = 0; ii < n_particles; ii++)
            {
                _nativeP[ii] = p[ii];
            }
            _job.Positions = _nativeP;
            _job.Schedule(_nativeMatrices.Length, 64).Complete();

            Graphics.RenderMeshInstanced(_rp, _mesh, 0, _nativeMatrices);
        }
    }
        
    private void OnDestroy()
    {
        _nativePositions.Dispose();
        _nativeMatrices.Dispose();
        _nativeCubeYOffsets.Dispose();
        _nativeP.Dispose();
    }
}

[BurstCompile]
internal struct ObjectPositionJob : IJobParallelFor
{
    public NativeArray<float3> Positions;
    public NativeArray<Matrix4x4> Matrices;
    public Vector3 Size;
    public Quaternion Rotation;

    public void Execute(int index)
    {
        Matrices[index] = Matrix4x4.TRS(Positions[index], Rotation, Size);
    }
}