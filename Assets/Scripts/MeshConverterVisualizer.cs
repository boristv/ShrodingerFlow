using System.Collections.Generic;
using UnityEngine;

public enum MeshType {
    Volume, Surface
};

public class MeshConverterVisualizer : MonoBehaviour
{
    [Header("Convert parameters")] 
    [SerializeField] private GameObject _model;
    [SerializeField] protected ComputeShader voxelizer;
    [SerializeField] protected int resolution = 32;
    [SerializeField] MeshType type = MeshType.Volume;

    [Header("Particles parameters")] 
    [SerializeField] private ParticleSystem _particleSystem;
    [SerializeField] private float _particleSize = 0.1125f;
    [SerializeField] private Color32 _particleColor = new Color32(255, 255, 255, 255);

    private void Start()
    {
        var points = MeshConverter.Convert(_model, resolution, voxelizer, type);
        Visualize(points);
    }

    private void Visualize(List<Vector3> points)
    {
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
}
