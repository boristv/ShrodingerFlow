using System.Collections.Generic;
using ShrodingerFlow.Particles;
using UnityEditor;
using UnityEditor.SceneManagement;
using UnityEngine;

/// <summary>
/// Сложный план этажа (как рис. 5.1–5.5 диссертации): здание 8×3×5, очаг слева, перегородки
/// со смещёнными проёмами, колонны, вытяжки справа. Тот же гибридный солвер SFHybridCS.
/// Старая простая комната (HybridSmoke) не трогается. Меню: ShrodingerFlow → Create Hybrid Building Scene.
/// </summary>
public static class HybridBuildingSceneSetup
{
    private const string BaseShaderDir = "Assets/Scenes/JetExample/ComputeShader/";
    private const string ScalarComputePath = "Assets/Scenes/HybridSmoke/Compute/SFHybridScalar.compute";
    private const string ParticlesShaderDir = "Assets/Particles/";
    private const string RaymarchShaderPath = "Assets/Scripts/Rendering/Raymarch/Raymarching.shader";
    private const string PsiDensityComputePath = "Assets/Particles/PsiDensityToVolume3D.compute";
    private const string ParticlesToDensityComputePath = "Assets/Particles/ParticlesToDensityVolume.compute";
    private const string SceneDir = "Assets/Scenes/HybridSmoke";
    private const string ScenePath = "Assets/Scenes/HybridSmoke/HybridBuilding.unity";

    // Здание: X 0..8 (длина), Y 0..3 (высота), Z 0..5 (глубина). Стены на всю высоту (y=1.5, half=1.5).
    private static readonly Vector3 RoomCenter = new Vector3(4f, 1.5f, 2.5f);
    private const float WY = 1.5f; // полувысота стен

    [MenuItem("ShrodingerFlow/Create Hybrid Building Scene")]
    public static void CreateBuildingScene()
    {
        var scene = EditorSceneManager.NewScene(NewSceneSetup.DefaultGameObjects, NewSceneMode.Single);

        SetupCamera();
        Light sun = SetupSun();

        var go = new GameObject("ISF_Building");
        go.transform.position = Vector3.zero;

        ConfigureDisplay(go, sun);
        AttachHybrid(go);

        if (!AssetDatabase.IsValidFolder(SceneDir))
            AssetDatabase.CreateFolder("Assets/Scenes", "HybridSmoke");
        EditorSceneManager.SaveScene(scene, ScenePath);
        Debug.Log($"[HybridBuildingSceneSetup] Scene saved: {ScenePath}");
        Debug.Log("Press Play. Edit the floor plan via SFHybridCS.wallSegments / extraVents (gizmos shown).");
    }

    private static void SetupCamera()
    {
        var cam = Camera.main;
        if (cam == null) return;
        var pos = new Vector3(13.5f, 7.0f, -5.5f);
        cam.transform.position = pos;
        cam.transform.rotation = Quaternion.LookRotation(RoomCenter - pos, Vector3.up);
        cam.clearFlags = CameraClearFlags.SolidColor;
        cam.backgroundColor = new Color(0.70f, 0.75f, 0.82f, 1f);
        cam.fieldOfView = 55f;
        cam.farClipPlane = 200f;
    }

    private static Light SetupSun()
    {
        var lightGo = GameObject.Find("Directional Light") ?? new GameObject("Directional Light");
        var light = lightGo.GetComponent<Light>() ?? lightGo.AddComponent<Light>();
        light.type = LightType.Directional;
        light.intensity = 1.1f;
        lightGo.transform.rotation = Quaternion.Euler(50f, -35f, 0f);
        return light;
    }

    private static void ConfigureDisplay(GameObject go, Light sun)
    {
        go.AddComponent<ParticleGpuBuffers>();
        var display = go.AddComponent<ParticleDisplay3D>();
        display.mode = ParticleDisplay3D.DisplayMode.Raymarch;
        display.raymarchDensitySource = ParticleDisplay3D.RaymarchDensitySource.AlphaField;
        display.shaderBillboard = AssetDatabase.LoadAssetAtPath<Shader>($"{ParticlesShaderDir}ParticleBillboard.shader");
        display.shaderShaded = AssetDatabase.LoadAssetAtPath<Shader>($"{ParticlesShaderDir}Particle3DSurf.shader");
        display.shaderRaymarch = AssetDatabase.LoadAssetAtPath<Shader>(RaymarchShaderPath);
        display.psiDensityToVolume = AssetDatabase.LoadAssetAtPath<ComputeShader>(PsiDensityComputePath);
        display.particlesToDensityVolume = AssetDatabase.LoadAssetAtPath<ComputeShader>(ParticlesToDensityComputePath);

        var so = new SerializedObject(display);
        var sunProp = so.FindProperty("_raymarchSunLight");
        if (sunProp != null) sunProp.objectReferenceValue = sun;
        so.ApplyModifiedPropertiesWithoutUndo();
    }

    private static void AttachHybrid(GameObject go)
    {
        var comp = go.AddComponent<SFHybridCS>();
        var so = new SerializedObject(comp);

        AssignBaseCompute(so, "_kernelsShader", "SFComputeKernels");
        AssignBaseCompute(so, "_fftShader", "SFComputeFFT");
        AssignBaseCompute(so, "_particlesShader", "SFComputeParticles");
        AssignBaseCompute(so, "_lesShader", "SFComputeLES");
        var scalar = AssetDatabase.LoadAssetAtPath<ComputeShader>(ScalarComputePath);
        var scalarProp = so.FindProperty("_scalarShader");
        if (scalarProp != null) scalarProp.objectReferenceValue = scalar;

        SetIntArray(so, "vol_size", 8, 3, 5);
        SetIntArray(so, "vol_res", 128, 32, 64);

        // Очаг слева, впуск-струя, вытяжка справа-вверху.
        SetV3(so, "_inletCenter", new Vector3(0.6f, 0.5f, 2.5f));
        SetV3(so, "_inletHalf", new Vector3(0.22f, 0.4f, 0.45f));
        SetV3(so, "_inletVelocity", new Vector3(0.9f, 0f, 0f));
        // Здание длинное → полный впуск во всё сечение, чтобы дым добивал до дальних комнат.
        SetBool(so, "_inflowFullFace", true);
        SetFloat(so, "_inflowDepth", 0.6f);
        SetV3(so, "_ventCenter", new Vector3(7.9f, 1.8f, 4.2f));
        SetV3(so, "_ventHalf", new Vector3(0.15f, 0.6f, 0.5f));
        SetV3(so, "_ventVelocity", new Vector3(0.6f, 0f, 0f));
        // Одиночное «мебельное» препятствие выключаем — план задают wallSegments.
        SetV3(so, "_obstacleHalf", Vector3.zero);

        so.ApplyModifiedPropertiesWithoutUndo();

        // План этажа — связный «серпантин» (проёмы: низ → верх → середина), без тупиков.
        // Колонны стоят в стороне от проёмов (перемешивают поток, не перекрывают путь).
        var segs = new List<HybridBox>();
        segs.AddRange(WallGapX(2.0f, 0.5f, 1.6f));  // стена 1: проём низкий
        segs.AddRange(WallGapX(4.0f, 3.4f, 4.5f));  // стена 2: проём высокий
        segs.AddRange(WallGapX(6.0f, 2.0f, 3.1f));  // стена 3: проём средний
        segs.Add(Col(3.0f, 3.2f));                  // комната 1, выше пути
        segs.Add(Col(5.0f, 1.5f));                  // комната 2, ниже пути
        segs.Add(Col(7.0f, 3.6f));                  // комната 3
        comp.wallSegments = segs.ToArray();

        comp.extraVents = new[]
        {
            VentBox(7.9f, 1.0f, 0.8f),  // нижняя вытяжка справа
        };
    }

    // Перегородка тонкая по X с проёмом по Z (gz0..gz1): два сегмента до и после проёма.
    private static IEnumerable<HybridBox> WallGapX(float x, float gz0, float gz1)
    {
        const float z0 = 0.1f, z1 = 4.9f; // интерьер по Z (за вычетом внешних стен)
        var list = new List<HybridBox>();
        if (gz0 > z0 + 0.05f)
            list.Add(new HybridBox { center = new Vector3(x, WY, (z0 + gz0) * 0.5f), half = new Vector3(0.08f, WY, (gz0 - z0) * 0.5f) });
        if (gz1 < z1 - 0.05f)
            list.Add(new HybridBox { center = new Vector3(x, WY, (gz1 + z1) * 0.5f), half = new Vector3(0.08f, WY, (z1 - gz1) * 0.5f) });
        return list;
    }

    private static HybridBox Col(float cx, float cz)
        => new HybridBox { center = new Vector3(cx, WY, cz), half = new Vector3(0.18f, WY, 0.18f) };

    private static HybridBox VentBox(float cx, float cy, float cz)
        => new HybridBox { center = new Vector3(cx, cy, cz), half = new Vector3(0.15f, 0.5f, 0.5f) };

    private static void SetV3(SerializedObject so, string name, Vector3 v)
    {
        var p = so.FindProperty(name);
        if (p != null) p.vector3Value = v;
    }

    private static void SetBool(SerializedObject so, string name, bool v)
    {
        var p = so.FindProperty(name);
        if (p != null) p.boolValue = v;
    }

    private static void SetFloat(SerializedObject so, string name, float v)
    {
        var p = so.FindProperty(name);
        if (p != null) p.floatValue = v;
    }

    private static void SetIntArray(SerializedObject so, string name, params int[] vals)
    {
        var p = so.FindProperty(name);
        if (p == null) return;
        p.arraySize = vals.Length;
        for (int i = 0; i < vals.Length; i++)
            p.GetArrayElementAtIndex(i).intValue = vals[i];
    }

    private static void AssignBaseCompute(SerializedObject so, string propertyName, string assetName)
    {
        var guids = AssetDatabase.FindAssets($"{assetName} t:ComputeShader", new[] { BaseShaderDir });
        if (guids.Length == 0) return;
        var path = AssetDatabase.GUIDToAssetPath(guids[0]);
        var shader = AssetDatabase.LoadAssetAtPath<ComputeShader>(path);
        var prop = so.FindProperty(propertyName);
        if (prop != null) prop.objectReferenceValue = shader;
    }
}
