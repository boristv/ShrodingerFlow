using System;
using Unity.Mathematics;
using UnityEngine;

public static class MovingController
{
    public const float DEPTH_OFFSET = 1;

    // Variance between mathf & mathematics perlin. Trying to keep them visually similar
    public const float HEIGHT_SCALE = 2.5f;
    public const float BURST_HEIGHT_SCALE = 1.5f;

    public const float NOISE_SCALE = 0.2f;
    private const int SIDE_LENGTH = 100;

    public static readonly quaternion RotGoal = quaternion.Euler(130, 50, 150);

    private static int Depth { get; set; } = 1;
    public static readonly Vector3 ObjectSize = Vector3.one;

    public static int GetCount => SIDE_LENGTH * SIDE_LENGTH * Depth;
    
    public static void SetPositions(Action<int, Vector3> setAction)
    {
        var i = 0;
        for (var y = 0; y < Depth; y++)
        {
            for (var x = 0; x < SIDE_LENGTH; x++)
            {
                for (var z = 0; z < SIDE_LENGTH; z++)
                {
                    setAction(i++, new Vector3(x, y, z));
                }
            }
        }
    }
    
    public static (Vector3 pos, Quaternion rot) CalculatePos(this Vector3 pos, float yOffset, float time)
    {
        var t = Mathf.InverseLerp(yOffset, HEIGHT_SCALE + yOffset, pos.y);
        var rot = Quaternion.Slerp(quaternion.identity, RotGoal, t);
        pos.y = HEIGHT_SCALE * Mathf.PerlinNoise(pos.x * NOISE_SCALE + time, pos.z * NOISE_SCALE + time) + yOffset * DEPTH_OFFSET;
        return (pos, rot);
    }

    public static (float3 pos, Quaternion rot) CalculatePosBurst(this float3 pos, float yOffset, float time)
    {
        var t = math.unlerp(yOffset, BURST_HEIGHT_SCALE + yOffset, pos.y);
        pos.y = BURST_HEIGHT_SCALE * noise.cnoise(new float2(pos.x * NOISE_SCALE + time, pos.z * NOISE_SCALE + time)) +
                yOffset * DEPTH_OFFSET;
        var rot = math.nlerp(quaternion.identity, RotGoal, t);
        return (pos, rot);
    }
}