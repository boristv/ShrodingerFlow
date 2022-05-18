using System.Numerics;
using ManagedCuda.VectorTypes;
using UnityEngine;
using source.assets.Discrete_space;
using source.assets.Particles;
using source;
using Vector2 = UnityEngine.Vector2;
using Vector3 = UnityEngine.Vector3;

public class SFJet : SFBase
{
    ParticleSystem.Particle[] cloud;
    bool bPointsUpdated = false;
    private ParticleSystem particleSystem;
    //PARAMETERS
    [Header("Начальные условия")]
    [SerializeField, Tooltip("Box size")] private int[] vol_size = { 4, 2, 2 };      // box size
    [SerializeField, Tooltip("Volume resolution")] private int[] vol_res = { 64, 32, 32 };    // volume resolution
    [SerializeField, Tooltip("Planck constant")] private float hbar = (float)0.1;           // Planck constant
    [SerializeField, Tooltip("Time step")] private float dt = 1 / (float)48;          // time step
    int tmax = 50;

    [SerializeField] private Vector3 jetVelocity = new Vector3(1f, 0f, 0f);
    [SerializeField] private Vector3 nozzleCen = new Vector3(2 - 1.7f, 1 - 0.034f, 1 + 0.066f);
    [SerializeField] private float nozzleLen = 0.5f;
    [SerializeField] private float nozzleRad = 0.5f;

    [SerializeField, Tooltip("Particles count")] private int n_particles = 50;

    private int particlesCount;

    Velocity vel;
    object lk= new object();

    UnityThreading.ActionThread myThread;

    [Header("Дополнительные настройки")] 
    [SerializeField, Tooltip("Размер частиц")] private float _particleSize = 0.1f;

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
        bool[,,] isJet = new bool[vol_res[0], vol_res[1], vol_res[2]];

        var tmp1 = new cuFloatComplex[ISF.properties.resx, ISF.properties.resy, ISF.properties.resz];
        var tmp2 = new cuFloatComplex[ISF.properties.resx, ISF.properties.resy, ISF.properties.resz];
        Complex tmp;
        for (int i = 0; i < vol_res[0]; i++)
        {
            for (int j = 0; j < vol_res[1]; j++)
            {
                for (int k = 0; k < vol_res[2]; k++)
                {
                    isJet[i,j,k] = Mathf.Abs(ISF.properties.px[i, j, k] - nozzleCen.x) <= nozzleLen/2f && 
                            Mathf.Pow(ISF.properties.py[i, j, k] - nozzleCen.y, 2) + Mathf.Pow(ISF.properties.pz[i, j, k] - nozzleCen.z, 2) <= Mathf.Pow(nozzleRad, 2);

                    tmp = Complex.One;
                    tmp1[i, j, k] = new cuFloatComplex((float)tmp.Real, (float)tmp.Imaginary);
                    tmp2[i, j, k] = new cuFloatComplex((float)(tmp.Real * 0.01), (float)(tmp.Imaginary * 0.01));
                }
            }
        }

        ISF.psi1.CopyToDevice(tmp1);
        ISF.psi2.CopyToDevice(tmp2);

        ISF.Normalize();
        //ISF.PressureProject();
        
        float[] kvec = { jetVelocity.x / hbar, jetVelocity.y / hbar, jetVelocity.z / hbar };
        var omega = (Mathf.Pow(jetVelocity.x, 2) + Mathf.Pow(jetVelocity.y, 2) + Mathf.Pow(jetVelocity.z, 2))/(2 * hbar);
        float phase;
        
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
                        if (isJet[i, j, k])
                        {
                            phase = kvec[0] * ISF.properties.px[i, j, k] +
                                    kvec[1] * ISF.properties.py[i, j, k] +
                                    kvec[2] * ISF.properties.pz[i, j, k];

                            var amp1 = Mathf.Sqrt(Mathf.Pow(tmp1[i, j, k].real, 2) + Mathf.Pow(tmp1[i, j, k].imag, 2));
                            var amp2 = Mathf.Sqrt(Mathf.Pow(tmp2[i, j, k].real, 2) + Mathf.Pow(tmp2[i, j, k].imag, 2));
                        
                            tmp = Complex.Exp(Complex.ImaginaryOne * phase);
                            tmp1[i, j, k] = new cuFloatComplex((float)(amp1 * tmp.Real), (float)(amp1 * tmp.Imaginary));
                            tmp2[i, j, k] = new cuFloatComplex((float)(amp2 * tmp.Real), (float)(amp2 * tmp.Imaginary));
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
        /*for (int i = 0; i < n_particles; i++)
        {
            x[i] = 0;
            y[i] = 0;
            z[i] = 0;
        }

        Particles.add_particles(x, y, z, n_particles);
        particlesCount += n_particles;*/

        vel = new Velocity(ISF.properties.resx, ISF.properties.resy, ISF.properties.resz);
        var c = new ParticleSystem.Particle[n_particles * 1000];
        
        /*for (int i = 0; i < n_particles; i++)
        {
            c[i].size = _particleSize;
        }*/

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
                        if (isJet[i, j, k])
                        {
                            phase = kvec[0] * ISF.properties.px[i, j, k] +
                                    kvec[1] * ISF.properties.py[i, j, k] +
                                    kvec[2] * ISF.properties.pz[i, j, k] - omega * dt * iterator;

                            var amp1 = Mathf.Sqrt(Mathf.Pow(tmp1[i, j, k].real, 2) + Mathf.Pow(tmp1[i, j, k].imag, 2));
                            var amp2 = Mathf.Sqrt(Mathf.Pow(tmp2[i, j, k].real, 2) + Mathf.Pow(tmp2[i, j, k].imag, 2));
                        
                            tmp = Complex.Exp(Complex.ImaginaryOne * phase);
                            tmp1[i, j, k] = new cuFloatComplex((float)(amp1 * tmp.Real), (float)(amp1 * tmp.Imaginary));
                            tmp2[i, j, k] = new cuFloatComplex((float)(amp2 * tmp.Real), (float)(amp2 * tmp.Imaginary));
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
            for (int i = 0; i < n_particles; i++)
            {
                var rt = (float)rnd.NextDouble() * 2 * Mathf.PI;
                x[i] = nozzleCen.x;
                y[i] = nozzleCen.y + 0.9f * nozzleRad * Mathf.Cos(rt);
                z[i] = nozzleCen.z + 0.9f * nozzleRad * Mathf.Sin(rt);
            }
            Particles.add_particles(x, y, z, n_particles);
            particlesCount += n_particles;
            

            //particle update
            ISF.update_velocities(vel);

            Particles.calculate_movement(vel);
            
            float[] px = Particles.x;
            float[] py = Particles.y;
            float[] pz = Particles.z;

            for (int ii = 0; ii < particlesCount; ++ii)
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
    }
}