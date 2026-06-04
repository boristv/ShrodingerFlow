using ShrodingerFlow.Particles;
using UnityEditor;
using UnityEditor.SceneManagement;
using UnityEngine;

/// <summary>
/// Сборка сцены гибридной волново-плотностной модели (ψ+α+LES), v1 изотермическая 3D-комната.
/// Меню: ShrodingerFlow → Create Hybrid Smoke Scene.
/// </summary>
public static class HybridSmokeSceneSetup
{
    private const string BaseShaderDir = "Assets/Scenes/JetExample/ComputeShader/";
    private const string ScalarComputePath = "Assets/Scenes/HybridSmoke/Compute/SFHybridScalar.compute";
    private const string ParticlesShaderDir = "Assets/Particles/";
    private const string RaymarchShaderPath = "Assets/Scripts/Rendering/Raymarch/Raymarching.shader";
    private const string PsiDensityComputePath = "Assets/Particles/PsiDensityToVolume3D.compute";
    private const string ParticlesToDensityComputePath = "Assets/Particles/ParticlesToDensityVolume.compute";
    private const string SceneDir = "Assets/Scenes/HybridSmoke";
    private const string ScenePath = "Assets/Scenes/HybridSmoke/HybridSmoke.unity";

    // Комната vol_size = (4,3,4), объект в начале координат.
    private static readonly Vector3 RoomCenter = new Vector3(2f, 1.5f, 2f);

    [MenuItem("ShrodingerFlow/Create Hybrid Smoke Scene")]
    public static void CreateHybridScene()
    {
        var scene = EditorSceneManager.NewScene(NewSceneSetup.DefaultGameObjects, NewSceneMode.Single);

        SetupCamera();
        Light sun = SetupSun();

        var go = new GameObject("ISF_Hybrid");
        go.transform.position = Vector3.zero;

        ConfigureDisplay(go, sun);
        AttachHybrid(go);

        if (!AssetDatabase.IsValidFolder(SceneDir))
            AssetDatabase.CreateFolder("Assets/Scenes", "HybridSmoke");
        EditorSceneManager.SaveScene(scene, ScenePath);
        Debug.Log($"[HybridSmokeSceneSetup] Scene saved: {ScenePath}");
        Debug.Log("Press Play. ParticleDisplay3D: mode=Raymarch, density source=α-field. Tune SFHybridCS in the Inspector.");
    }

    private static void SetupCamera()
    {
        var cam = Camera.main;
        if (cam == null) return;
        var pos = new Vector3(8.5f, 5.5f, -3.5f);
        cam.transform.position = pos;
        cam.transform.rotation = Quaternion.LookRotation(RoomCenter - pos, Vector3.up);
        cam.clearFlags = CameraClearFlags.SolidColor;
        // Фон, который раймарч читает как _CameraOpaqueTexture (светло-голубой).
        cam.backgroundColor = new Color(0.70f, 0.75f, 0.82f, 1f);
        cam.fieldOfView = 55f;
        cam.farClipPlane = 100f;
    }

    private static Light SetupSun()
    {
        var lightGo = GameObject.Find("Directional Light");
        if (lightGo == null)
        {
            lightGo = new GameObject("Directional Light");
            lightGo.AddComponent<Light>().type = LightType.Directional;
        }
        var light = lightGo.GetComponent<Light>();
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

        // _raymarchSunLight — приватное сериализуемое поле.
        var so = new SerializedObject(display);
        var sunProp = so.FindProperty("_raymarchSunLight");
        if (sunProp != null) sunProp.objectReferenceValue = sun;
        so.ApplyModifiedPropertiesWithoutUndo();

        if (display.shaderRaymarch == null)
            Debug.LogWarning($"[HybridSmokeSceneSetup] Shader not found: {RaymarchShaderPath}");
        if (display.psiDensityToVolume == null)
            Debug.LogWarning($"[HybridSmokeSceneSetup] ComputeShader not found: {PsiDensityComputePath}");
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
        if (scalar == null)
            Debug.LogWarning($"[HybridSmokeSceneSetup] ComputeShader not found: {ScalarComputePath}");

        so.ApplyModifiedPropertiesWithoutUndo();
    }

    private static void AssignBaseCompute(SerializedObject so, string propertyName, string assetName)
    {
        var guids = AssetDatabase.FindAssets($"{assetName} t:ComputeShader", new[] { BaseShaderDir });
        if (guids.Length == 0)
        {
            Debug.LogWarning($"[HybridSmokeSceneSetup] ComputeShader '{assetName}' not found in {BaseShaderDir}");
            return;
        }
        var path = AssetDatabase.GUIDToAssetPath(guids[0]);
        var shader = AssetDatabase.LoadAssetAtPath<ComputeShader>(path);
        var prop = so.FindProperty(propertyName);
        if (prop != null) prop.objectReferenceValue = shader;
    }
}
