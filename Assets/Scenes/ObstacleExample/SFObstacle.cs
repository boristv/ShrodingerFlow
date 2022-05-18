using System.Numerics;
using ManagedCuda.VectorTypes;
using UnityEngine;
using source.assets.Discrete_space;
using source.assets.Particles;
using source;
using Vector2 = UnityEngine.Vector2;
using Vector3 = UnityEngine.Vector3;

public class SFObstacle : SFBase
{
    ParticleSystem.Particle[] cloud;
    bool bPointsUpdated = false;
    private ParticleSystem particleSystem;
    //PARAMETERS
    [Header("Начальные условия")]
    [SerializeField, Tooltip("Box size")] private int[] vol_size = { 4, 2, 2 };      // box size
    [SerializeField, Tooltip("Volume resolution")] private int[] vol_res = { 64, 32, 32 };    // volume resolution
    [SerializeField, Tooltip("Planck constant")] private float hbar = (float)0.1;           // Planck constant
    [SerializeField, Tooltip("Time step")] private float dt = 1 / (float)12;          // time step
    int tmax = 50;

    [SerializeField] private Vector3 obstaclePosition = new Vector3(1.5f, 1, 1);
    [SerializeField] private Vector3 obstacleSize = new Vector3(0.5f, 0.5f, 0.5f);
    [SerializeField, Tooltip("Background velocity")] private float[] background_vel = {(float) -0.2, 0, 0};

    [SerializeField, Tooltip("Particles count")] private int n_particles = 50;

    private int particlesCount;

    Velocity vel;
    object lk= new object();

    UnityThreading.ActionThread myThread;

    [Header("Дополнительные настройки")] 
    [SerializeField, Tooltip("Размер частиц")] private float _particleSize = 0.1f;
    
    [Header("Начальное расположение частиц")] 
    [SerializeField] private Vector2 _boxSizeX = new Vector2(5, 5);
    [SerializeField] private Vector2 _boxSizeY = new Vector2(0.5f, 4.5f);
    [SerializeField] private Vector2 _boxSizeZ = new Vector2(0.5f, 4.5f);
    
    [SerializeField] private Vector3 nozzleCen = new Vector3(2 - 1.7f, 1 - 0.034f, 1 + 0.066f);
    [SerializeField] private float nozzleLen = 0.5f;
    [SerializeField] private float nozzleRad = 0.5f;

    private Vector3 _volSizeV3 = new Vector3(4, 2, 2);
    
    private void Start()
    {
        _volSizeV3 = new Vector3(vol_size[0], vol_size[1], vol_size[2]);
        
        particleSystem = GetComponent<ParticleSystem>();
        myThread = UnityThreadHelper.CreateThread(Worker);
    }

    private void Worker()
    {
        //INITIALISATION
        ISF.Init(vol_size, vol_res, hbar, dt);
        Particles.init(n_particles * 1000);

        //init psi
        bool[,,] isObstacle = new bool[vol_res[0], vol_res[1], vol_res[2]];
        bool[,,] isJet = new bool[vol_res[0], vol_res[1], vol_res[2]];

        float[] kvec = { background_vel[0] / hbar, background_vel[1] / hbar, background_vel[2] / hbar };
        var omega = (Mathf.Pow(background_vel[0], 2) + Mathf.Pow(background_vel[1], 2) + Mathf.Pow(background_vel[2], 2))/(2 * hbar);
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
                    //box
                    /*isObstacle[i, j, k] =
                        Mathf.Abs(ISF.properties.px[i, j, k] - obstaclePosition.x) <= obstacleSize.x / 2f &&
                        Mathf.Abs(ISF.properties.py[i, j, k] - obstaclePosition.y) <= obstacleSize.y / 2f &&
                        Mathf.Abs(ISF.properties.pz[i, j, k] - obstaclePosition.z) <= obstacleSize.z / 2f;//*/
                    //sphere
                    isObstacle[i, j, k] =
                        Mathf.Pow(ISF.properties.px[i, j, k] - obstaclePosition.x, 2) + 
                        Mathf.Pow(ISF.properties.py[i, j, k] - obstaclePosition.y, 2) + 
                        Mathf.Pow(ISF.properties.pz[i, j, k] - obstaclePosition.z, 2)
                        <= Mathf.Pow(obstacleSize.x, 2);//*/

                    phase = kvec[0] * ISF.properties.px[i, j, k] +
                                 kvec[1] * ISF.properties.py[i, j, k] +
                                 kvec[2] * ISF.properties.pz[i, j, k];
                    tmp = Complex.Exp(Complex.ImaginaryOne * phase);
                    tmp1[i, j, k] = new cuFloatComplex((float)tmp.Real, (float)tmp.Imaginary);
                    tmp2[i, j, k] = new cuFloatComplex((float)(tmp.Real * 0.01), (float)(tmp.Imaginary * 0.01));
                }
            }
        }

        ISF.psi1.CopyToDevice(tmp1);
        ISF.psi2.CopyToDevice(tmp2);

        ISF.Normalize();
        //ISF.PressureProject();
        
        
        
        for (int iter = 0;  iter < 10; iter++)
        {
            ISF.psi1.CopyToHost(tmp1);
            ISF.psi2.CopyToHost(tmp2);
            for (int i = 0; i < vol_res[0]; i++)
            {
                for (int j = 0; j < vol_res[1]; j++)
                {
                    for (int k = 0; k < vol_res[2]; k++)
                    {
                        if (isObstacle[i, j, k])
                        {
                            phase = 0;
                            
                            var amp1 = Mathf.Sqrt(Mathf.Pow(tmp1[i, j, k].real, 2) + Mathf.Pow(tmp1[i, j, k].imag, 2));
                            var amp2 = Mathf.Sqrt(Mathf.Pow(tmp2[i, j, k].real, 2) + Mathf.Pow(tmp2[i, j, k].imag, 2));


                            tmp = Complex.Exp(Complex.ImaginaryOne * phase);
                            tmp1[i, j, k] = new cuFloatComplex((float)(amp1 * tmp.Real), (float)(amp1 * tmp.Imaginary));
                            tmp2[i, j, k] = new cuFloatComplex((float)(amp2 * tmp.Real), (float)(amp2 * tmp.Imaginary));
                            //tmp2[i, j, k] = new cuFloatComplex((float)(0.01f), (float)(0));
                        }
                    }
                }
            }
            ISF.psi1.CopyToDevice(tmp1);
            ISF.psi2.CopyToDevice(tmp2);
            ISF.PressureProject();
        }
        

        //init particles
        var x = new float[n_particles];
        var y = new float[n_particles];
        var z = new float[n_particles];

        vel = new Velocity(ISF.properties.resx, ISF.properties.resy, ISF.properties.resz);
        var c = new ParticleSystem.Particle[n_particles * 1000];

        iterator = 0;
        while (true)
        {
            iterator++;
            //MAIN ITERATION
            //incompressible Schroedinger flow
            ISF.update_space();

            ISF.psi1.CopyToHost(tmp1);
            ISF.psi2.CopyToHost(tmp2);
            for (int i = 0; i < vol_res[0]; i++)
            {
                for (int j = 0; j < vol_res[1]; j++)
                {
                    for (int k = 0; k < vol_res[2]; k++)
                    {
                        if (isObstacle[i, j, k])
                        {
                            phase = 0;//- omega * dt * iterator;
                            
                            var amp1 = Mathf.Sqrt(Mathf.Pow(tmp1[i, j, k].real, 2) + Mathf.Pow(tmp1[i, j, k].imag, 2));
                            var amp2 = Mathf.Sqrt(Mathf.Pow(tmp2[i, j, k].real, 2) + Mathf.Pow(tmp2[i, j, k].imag, 2));
                            
                            tmp = Complex.Exp(Complex.ImaginaryOne * phase);
                            
                            tmp1[i, j, k] = new cuFloatComplex((float)(amp1 * tmp.Real), (float)(amp1 * tmp.Imaginary));
                            tmp2[i, j, k] = new cuFloatComplex((float)(amp2 * tmp.Real), (float)(amp2 * tmp.Imaginary));
                            //tmp2[i, j, k] = new cuFloatComplex((float)(0.01f), (float)(0));
                        }
                    }
                }
            }
            
            ISF.psi1.CopyToDevice(tmp1);
            ISF.psi2.CopyToDevice(tmp2);
            
            ISF.PressureProject();
            
            
            //init particles
            x = new float[n_particles];
            y = new float[n_particles];
            z = new float[n_particles];
            System.Random rnd = new System.Random();
            //if (iterator < 20)
            {
                for (int i = 0; i < n_particles; i++)
                {
                    var rt = (float)rnd.NextDouble() * 2 * Mathf.PI;
                    x[i] = nozzleCen.x;
                    y[i] = nozzleCen.y + 0.9f * nozzleRad * Mathf.Cos(rt);
                    z[i] = nozzleCen.z + 0.9f * nozzleRad * Mathf.Sin(rt);
                    //y[i] = (float) (rnd.NextDouble() * (_boxSizeY.y - _boxSizeY.x) + _boxSizeY.x);
                    //z[i] = (float) (rnd.NextDouble() * (_boxSizeZ.y - _boxSizeZ.x) + _boxSizeZ.x);
                }
                Particles.add_particles(x, y, z, n_particles);
                particlesCount += n_particles;
            }
            

            //particle update
            ISF.update_velocities(vel);

            Particles.calculate_movement(vel);
            
            float[] px = Particles.x;
            float[] py = Particles.y;
            float[] pz = Particles.z;

            for (int ii = 0; ii < particlesCount; ++ii) //< n_particles;
            {
                var pos = new Vector3(px[ii], py[ii], pz[ii]);
                var lastPose = c[ii].position;
                c[ii].position = pos;
                c[ii].velocity = (pos - lastPose);
                var cc = c[ii].velocity.normalized;
                var alpha = 1;//cc.magnitude == 0 ? 0 : 1;
                c[ii].color = new Color(cc.x, cc.y, cc.z, alpha);

                c[ii].size = _particleSize;
            }

            lock (pz) {
                cloud = c;
            }

            if (UnityThreading.ActionThread.CurrentThread.ShouldStop)
                return; // finish
        }

    }

    public void Update()
    {
        lock (lk) {
            if (cloud != null) {
                particleSystem.SetParticles(cloud, cloud.Length);
            }
        }
    }
    
    [ExecuteAlways]
    private void OnDrawGizmos()
    {
        // Draw a yellow sphere at the transform's position
        Gizmos.color = Color.yellow;
        Gizmos.DrawWireCube(transform.position + _volSizeV3 / 2f, _volSizeV3);
        Gizmos.color = new Color(1, 1, 0, 0.5f);
        Gizmos.DrawWireSphere(transform.position + obstaclePosition, obstacleSize.x);
    }
}