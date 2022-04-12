using System.Collections;
using System.Collections.Generic;
using source.assets.Particles;
using UnityEngine;
using VoxelSystem;

public class MeshConverter : MonoBehaviour
{
    [SerializeField] private MeshFilter _meshFilter;
    [SerializeField] private bool _manyMeshes;
    [SerializeField] private List<MeshFilter> _meshFilters;
    [SerializeField] private int _particlesCount = 100000;
    [SerializeField] private ParticleSystem _particleSystem;
    [SerializeField] protected ComputeShader voxelizer;
    [SerializeField] protected int resolution = 32;
    [SerializeField] MeshType type = MeshType.Volume;

    [Header("Particles parameters")] 
    [SerializeField] private float _particleSize = 0.1125f;
    [SerializeField] private Color32 _particleColor = new Color32(255, 255, 255, 255);
    
    UnityThreading.ActionThread myThread;
    ParticleSystem.Particle[] cloud;
    object lk= new object();
    
    enum MeshType {
        Volume, Surface
    };
    
    private void Start()
    {
        Convert();
    }

    private void Convert()
    {
        /*Particles.init(_particlesCount);

        //init particles
        var x = new float[_particlesCount];
        var y = new float[_particlesCount];
        var z = new float[_particlesCount];

        var points = PointsFromMesh(_meshFilter.mesh);
            
        for (int i = 0; i < _particlesCount; i++)
        {
            x[i] = points[i].x;
            y[i] = points[i].y;
            z[i] = points[i].z;
        }

        Particles.add_particles(x, y, z, _particlesCount);

        var c = new ParticleSystem.Particle[_particlesCount];*/

        List<Vector3> points;
        if (_manyMeshes)
        {
            points = PointsFromMesh(_meshFilters);
        }
        else
        {
            points = PointsFromMesh(_meshFilter.mesh);
        }
        var particles = new ParticleSystem.Particle[points.Count];
        for (var i = 0; i < particles.Length; i++)
        {
            particles[i].position = points[i];
            particles[i].velocity = Vector3.zero;
            particles[i].startSize = _particleSize;
            particles[i].startColor = _particleColor;
        }

        _particleSystem.SetParticles(particles, particles.Length);
    }

    private List<Vector3> PointsFromMesh(Mesh mesh)
    {
        var list = new List<Vector3>();
        var data = GPUVoxelizer.Voxelize(voxelizer, mesh, resolution, (type == MeshType.Volume));
        foreach (var voxel in data.GetData())
        {
            if (voxel.fill > 0)
            {
                list.Add(voxel.position);
            }
        }
        data.Dispose();
        Debug.Log($"Particles count = {list.Count}");
        return list;
    }
    
    private List<Vector3> PointsFromMesh(List<MeshFilter> meshFilters)
    {
        var list = new List<Vector3>();
        foreach (var meshFilter in meshFilters)
        {
            var data = GPUVoxelizer.Voxelize(voxelizer, meshFilter.mesh, resolution, (type == MeshType.Volume));
            foreach (var voxel in data.GetData())
            {
                if (voxel.fill > 0)
                {
                    list.Add(voxel.position);
                }
            }
            data.Dispose();
        }
        Debug.Log($"Particles count = {list.Count}");
        return list;
    }
}
