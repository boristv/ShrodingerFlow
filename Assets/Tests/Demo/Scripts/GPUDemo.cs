using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace VoxelSystem.Demo {

    [RequireComponent (typeof(MeshFilter))]
	public class GPUDemo : MonoBehaviour
	{
		[SerializeField] private Transform _parent;
		[SerializeField] private GameObject _prefab;

        enum MeshType {
            Volume, Surface
        };

        [SerializeField] MeshType type = MeshType.Volume;
		[SerializeField] protected Mesh mesh;
		[SerializeField] protected ComputeShader voxelizer;
		[SerializeField] protected int resolution = 32;
        [SerializeField] protected bool useUV = false;

        [SerializeField] private ParticleSystem _particleSystem;
        [SerializeField] private int n_particles = 100000;

		void Start () {
			var data = GPUVoxelizer.Voxelize(voxelizer, mesh, resolution, (type == MeshType.Volume));
			foreach (var voxel in data.GetData())
			{
				if (voxel.fill > 0)
				{
					var go = Instantiate(_prefab, voxel.position, Quaternion.identity);
					//go.transform.localScale = Vector3.one * 300;
					go.transform.parent = _parent;
				}
			}
			_parent.localScale = 1 * Vector3.one;
			GetComponent<MeshFilter>().sharedMesh = VoxelMesh.Build(data.GetData(), data.UnitLength, useUV);
			data.Dispose();
			
			var c = new ParticleSystem.Particle[n_particles];
        
			for (int i = 0; i < n_particles; i++)
			{
				//c[i].position = 
			}
		}

	}

}
