using System;
using UnityEngine;

/// <summary>
/// Данные пресетов для контекстного меню «Apply Scenario Defaults».
/// Редактируются в SFUnifiedScenarioPresets.asset (источник правды вместо констант в SFUnifiedCS).
/// SphereObstacle — копия полей из ObstacleExample.unity (SFObstacle).
/// </summary>
[Serializable]
public struct SFUnifiedJetPreset
{
    public int[] vol_size;
    public int[] vol_res;
    public float hbar;
    public float dt;
    public Vector3 velocity;
    public Vector3 nozzleCen;
    public float nozzleLen;
    public float nozzleRad;
    public int nParticles;
    public float particleSize;
    public int stepsPerFrame;
    public bool useLES;
}

[Serializable]
public struct SFUnifiedSphereObstaclePreset
{
    public int[] vol_size;
    public int[] vol_res;
    public float hbar;
    public float dt;
    public Vector3 velocity;
    public Vector3 obstaclePos1;
    public float obstacleRadius1;
    public Vector3 nozzleCen;
    public float nozzleLen;
    public float nozzleRad;
    public Vector2 boxSpawnX;
    public Vector2 boxSpawnY;
    public Vector2 boxSpawnZ;
    public int nParticles;
    public float particleSize;
    public int stepsPerFrame;
    public bool useLES;
}

[Serializable]
public struct SFUnifiedCylinderObstaclePreset
{
    public int[] vol_size;
    public int[] vol_res;
    public float hbar;
    public float dt;
    public Vector3 velocity;
    public Vector3 obstaclePos1;
    public float obstacleRadius1;
    public Vector3 nozzleCen;
    public Vector2 boxSpawnX;
    public Vector2 boxSpawnY;
    public Vector2 boxSpawnZ;
    public int nParticles;
    public float particleSize;
    public int stepsPerFrame;
    public bool useLES;
}

[Serializable]
public struct SFUnifiedTwoSpheresPreset
{
    public int[] vol_size;
    public int[] vol_res;
    public float hbar;
    public float dt;
    public Vector3 velocity;
    public Vector3 obstaclePos1;
    public float obstacleRadius1;
    public Vector3 obstaclePos2;
    public float obstacleRadius2;
    public Vector2 boxSpawnX;
    public Vector2 boxSpawnY;
    public Vector2 boxSpawnZ;
    public int nParticles;
    public float particleSize;
    public int stepsPerFrame;
    public bool useLES;
}

[Serializable]
public struct SFUnifiedLeapfrogRingsPreset
{
    public int[] vol_size;
    public int[] vol_res;
    public float hbar;
    public float dt;
    public Vector3 velocity;
    public float ring1Radius;
    public float ring2Radius;
    public Vector3 ring1Normal;
    public Vector3 ring2Normal;
    public Vector2 boxSpawnX;
    public Vector2 boxSpawnY;
    public Vector2 boxSpawnZ;
    public int nParticles;
    public float particleSize;
    public int stepsPerFrame;
    public bool useLES;
}
