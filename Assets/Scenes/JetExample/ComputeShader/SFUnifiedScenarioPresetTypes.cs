using System;
using UnityEngine;

/// <summary>
/// Данные пресетов для контекстного меню «Apply Scenario Defaults».
/// Редактируются в SFUnifiedScenarioPresets.asset (источник правды вместо констант в SFUnifiedCS).
/// SphereObstacle — ObstacleExample (SFObstacle); CylinderObstacle — OtherObstacleExample (OtherSFObstacle); TwoSpheres — дефолты SFUnifiedCS.
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

/// <summary>Пресет «ёмкость»: домен = сосуд; блок жидкости; стенки по периметру; гравитация/вязкость — поля SFUnifiedCS.</summary>
[Serializable]
public struct SFUnifiedRectangularContainerPreset
{
    public int[] vol_size;
    public int[] vol_res;
    public float hbar;
    public float dt;
    public Vector3 fluidMin;
    public Vector3 fluidMax;
    public float wallThickness;
    public bool applyPsi2Gravity;
    public Vector3 psi2Gravity;
    public float kinematicViscosity;
    public int nParticles;
    public float particleSize;
    public int stepsPerFrame;
    public bool useLES;
    public bool useLiquidChiField;
    public float liquidChiThreshold;
}

/// <summary>Пресет «дым в помещении»: 2 комнаты с проёмом, источник дыма (χ + плавучесть), вытяжка (сток χ + подсос).</summary>
[Serializable]
public struct SFUnifiedRoomSmokePreset
{
    public int[] vol_size;
    public int[] vol_res;
    public float hbar;
    public float dt;
    public float wallThickness;
    public float partitionX;
    public float partitionThickness;
    public float doorCenterZ;
    public float doorWidth;
    public float doorHeight;
    public Vector3 sourceCenter;
    public Vector3 sourceHalf;
    public Vector3 emitVelocity;
    public float chiInject;
    public Vector3 ventCenter;
    public Vector3 ventHalf;
    public Vector3 ventSuction;
    public float ventDecay;
    public float buoyancyBeta;
    public Vector3 buoyancyDir;
    public float smokeRiseSpeed;
    public float smokeDiffusion;
    public float tracerDispersion;
    public float turbAmplitude;
    public float turbScale;
    public float kinematicViscosity;
    public int nParticles;
    public float particleSize;
    public int stepsPerFrame;
    public bool useLES;
}
