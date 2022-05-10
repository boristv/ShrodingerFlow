using System.Collections.Generic;
using UnityEngine;
using VoxelSystem;

public class MeshConverter : MonoBehaviour
{
    [SerializeField] private GameObject _model;
    [SerializeField] private ParticleSystem _particleSystem;
    [SerializeField] protected ComputeShader voxelizer;
    [SerializeField] protected int resolution = 32;
    [SerializeField] MeshType type = MeshType.Volume;

    [Header("Particles parameters")] 
    [SerializeField] private float _particleSize = 0.1125f;
    [SerializeField] private Color32 _particleColor = new Color32(255, 255, 255, 255);
    [SerializeField] private bool _needScaling;
    [SerializeField] private float _scaling = 0.01f;

    private struct MeshFilterAndBounds
    {
        public MeshFilter MeshFilter;
        public Bounds Bounds;
    }
    
    private List<MeshFilterAndBounds> meshFilterAndBoundsList;
    private Bounds FullBounds;
    
    enum MeshType {
        Volume, Surface
    };
    
    private void Start()
    {
        Convert();
    }

    private void Convert()
    {
        var points = PointsFromMesh(_model);
        
        if (points == null || points.Count == 0)
        {
            Debug.Log("Convert failed: no points");
            return;
        }
        
        //Normalizing and scaling
        if (_needScaling)
        {
            var dist = Vector3.Distance(points[0], points[1]);
            
            for (var i = 0; i < points.Count; i++)
            {
                points[i] /= dist;
                points[i] *= _scaling;
            }
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

    private List<Vector3> PointsFromMesh(GameObject meshesParent)
    {
        var points = new List<Vector3>();
        
        MeshFilter[] meshFilters = meshesParent.GetComponentsInChildren<MeshFilter>();

        if (meshFilters.Length == 0)
        {
            Debug.Log("Model does not contain meshes");
            return new List<Vector3>();
        }
        
        meshFilterAndBoundsList = new List<MeshFilterAndBounds>();
        foreach (var meshFilter in meshFilters)
        {
            meshFilter.sharedMesh.RecalculateBounds();
            var meshFilterAndBounds = new MeshFilterAndBounds()
            {
                MeshFilter = meshFilter,
                Bounds = ModifiedGPUVoxelizer.GetRealBounds(meshFilter)
            };
            meshFilterAndBoundsList.Add(meshFilterAndBounds);
        }
        FullBounds = meshFilterAndBoundsList[0].Bounds;
        foreach (var meshFilterAndBounds in meshFilterAndBoundsList)
        {
            FullBounds.Encapsulate(meshFilterAndBounds.Bounds);
        }
        
        foreach (var meshFilterAndBounds in meshFilterAndBoundsList)
        {
            var data = ModifiedGPUVoxelizer.Voxelize(voxelizer, meshFilterAndBounds.MeshFilter, meshFilterAndBounds.Bounds, FullBounds, resolution, (type == MeshType.Volume));
            foreach (var voxel in data.GetData())
            {
                if (voxel.fill > 0)
                {
                    points.Add(voxel.position);
                }
            }
            data.Dispose();
        }
        Debug.Log($"Particles count = {points.Count}");
        return points;
    }

    /*private void OnDrawGizmos()
    {
        if (meshFilterAndBoundsList != null)
        {
            foreach (var element in meshFilterAndBoundsList)
            {
                Gizmos.DrawWireCube(element.Bounds.center, element.Bounds.extents * 2);
            }
        }
        Gizmos.color = Color.yellow;
        Gizmos.DrawWireCube(FullBounds.center, FullBounds.extents * 2);
    }*/
}
