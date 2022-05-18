using System.Collections.Generic;
using UnityEngine;
using VoxelSystem;

public class MeshConverter
{
    private struct MeshFilterAndBounds
    {
        public MeshFilter MeshFilter;
        public Bounds Bounds;
    }

    public static List<Vector3> Convert(GameObject meshesParent, int resolution, ComputeShader voxelizer, MeshType type)
    {
        var points = new List<Vector3>();
        
        MeshFilter[] meshFilters = meshesParent.GetComponentsInChildren<MeshFilter>();

        if (meshFilters.Length == 0)
        {
            Debug.Log("Model does not contain meshes");
            return new List<Vector3>();
        }
        
        var meshFilterAndBoundsList = new List<MeshFilterAndBounds>();
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
        var fullBounds = meshFilterAndBoundsList[0].Bounds;
        foreach (var meshFilterAndBounds in meshFilterAndBoundsList)
        {
            fullBounds.Encapsulate(meshFilterAndBounds.Bounds);
        }
        
        foreach (var meshFilterAndBounds in meshFilterAndBoundsList)
        {
            var data = ModifiedGPUVoxelizer.Voxelize(voxelizer, meshFilterAndBounds.MeshFilter, meshFilterAndBounds.Bounds, fullBounds, resolution, (type == MeshType.Volume));
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
    
    public static bool TryConvertToGrid(List<Vector3> points, out bool[,,] grid, out Vector3 gridSize)
    {
        if (points == null || points.Count == 0)
        {
            Debug.Log("Convert failed: no points");
            grid = null;
            gridSize = Vector3.zero;
            return false;
        }
        
        int sizeX = 0, sizeY = 0, sizeZ = 0;
        
        var pointsIndex = new List<Vector3Int>();
        var dist = Vector3.Distance(points[0], points[1]);

        var minX = int.MaxValue;
        var minY = int.MaxValue;
        var minZ = int.MaxValue;
        var maxX = int.MinValue;
        var maxY = int.MinValue;
        var maxZ = int.MinValue;
            
        for (var i = 0; i < points.Count; i++)
        {
            points[i] /= dist;
            var x = Mathf.RoundToInt(points[i].x);
            var y = Mathf.RoundToInt(points[i].y);
            var z = Mathf.RoundToInt(points[i].z);
            pointsIndex.Add(new Vector3Int(x, y, z));
            if (x < minX) minX = x;
            if (y < minY) minY = y;
            if (z < minZ) minZ = z;
            if (x > maxX) maxX = x;
            if (y > maxY) maxY = y;
            if (z > maxZ) maxZ = z;
        }

        sizeX = maxX - minX + 1;
        sizeY = maxY - minY + 1;
        sizeZ = maxZ - minZ + 1;

        gridSize = new Vector3(sizeX, sizeY, sizeZ);

        for (int i = 0; i < pointsIndex.Count; i++)
        {
            pointsIndex[i] -= new Vector3Int(minX, minY, minZ);
        }
            
        grid = new bool[sizeX, sizeY, sizeZ];
        Debug.Log($"size: x = {sizeX} | y = {sizeY} | z = {sizeZ}");

        for (int i = 0; i < pointsIndex.Count; i++)
        {
            grid[pointsIndex[i].x, pointsIndex[i].y, pointsIndex[i].z] = true;
        }

        return true;
    }

    public static List<Vector3> Rotate(List<Vector3> points, Vector3 rotation)
    {
        var rotPoints = new List<Vector3>();
        foreach (var point in points)
        {
            rotPoints.Add(ModifiedGPUVoxelizer.Rotate(point, rotation));
        }
        return rotPoints;
    }
}
