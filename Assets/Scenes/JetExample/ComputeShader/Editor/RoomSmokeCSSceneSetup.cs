using ShrodingerFlow.Particles;
using UnityEditor;
using UnityEditor.SceneManagement;
using UnityEngine;

/// <summary>
/// Создаёт сцену цифрового двойника распространения дыма в помещении (SFUnifiedCS, сценарий RoomSmoke):
/// две комнаты с дверным проёмом, источник дыма, вытяжка. Привязывает все нужные compute-шейдеры
/// (включая контейнерный ISF и стеночный шейдер частиц), частицы и объёмную визуализацию.
/// </summary>
public static class RoomSmokeCSSceneSetup
{
    private const string ShaderDir = "Assets/Scenes/JetExample/ComputeShader/";
    private const string ParticlesShaderDir = "Assets/Particles/";
    private const string RaymarchShaderPath = "Assets/Scripts/Rendering/Raymarch/Raymarching.shader";
    private const string PsiDensityComputePath = "Assets/Particles/PsiDensityToVolume3D.compute";
    private const string ParticlesToDensityComputePath = "Assets/Particles/ParticlesToDensityVolume.compute";
    private const string ContainerIsfPath = ShaderDir + "SFComputeKernelsContainer.compute";
    private const string ParticlesWallPath = ShaderDir + "SFComputeParticlesWall.compute";
    private const string ScenePath = ShaderDir + "RoomSmokeCS.unity";

    [MenuItem("ShrodingerFlow/Create Room Smoke Scene")]
    public static void CreateRoomSmokeScene()
    {
        var scene = EditorSceneManager.NewScene(NewSceneSetup.DefaultGameObjects, NewSceneMode.Single);

        SetupCamera();

        var go = new GameObject("ISF_RoomSmoke");
        go.transform.position = Vector3.zero;

        ConfigureGpuParticles(go);
        AttachUnifiedCS(go);

        EditorSceneManager.SaveScene(scene, ScenePath);
        Debug.Log($"[RoomSmokeCSSceneSetup] Scene saved: {ScenePath}");
        Debug.Log("Выберите ISF_RoomSmoke и при необходимости вызовите контекстное меню 'Apply Room smoke defaults'.");
    }

    private static void SetupCamera()
    {
        var cam = Camera.main;
        if (cam == null) return;
        cam.transform.position = new Vector3(2f, 4f, -6f);
        cam.transform.rotation = Quaternion.Euler(20f, 0f, 0f);
        cam.clearFlags = CameraClearFlags.SolidColor;
        cam.backgroundColor = new Color(0.05f, 0.06f, 0.09f, 0f);
        cam.fieldOfView = 60f;
    }

    private static void ConfigureGpuParticles(GameObject go)
    {
        go.AddComponent<ParticleGpuBuffers>();

        var display = go.AddComponent<ParticleDisplay3D>();
        display.mode = ParticleDisplay3D.DisplayMode.Billboard;
        display.scale = 0.05f * 50f;

        var bb = AssetDatabase.LoadAssetAtPath<Shader>($"{ParticlesShaderDir}ParticleBillboard.shader");
        var surf = AssetDatabase.LoadAssetAtPath<Shader>($"{ParticlesShaderDir}Particle3DSurf.shader");
        display.shaderBillboard = bb;
        display.shaderShaded = surf;

        display.shaderRaymarch = AssetDatabase.LoadAssetAtPath<Shader>(RaymarchShaderPath);
        display.psiDensityToVolume = AssetDatabase.LoadAssetAtPath<ComputeShader>(PsiDensityComputePath);
        display.particlesToDensityVolume = AssetDatabase.LoadAssetAtPath<ComputeShader>(ParticlesToDensityComputePath);

        if (bb == null)
            Debug.LogWarning($"[RoomSmokeCSSceneSetup] Shader not found: {ParticlesShaderDir}ParticleBillboard.shader");
        if (surf == null)
            Debug.LogWarning($"[RoomSmokeCSSceneSetup] Shader not found: {ParticlesShaderDir}Particle3DSurf.shader");
        if (display.shaderRaymarch == null)
            Debug.LogWarning($"[RoomSmokeCSSceneSetup] Shader not found: {RaymarchShaderPath}");
        if (display.psiDensityToVolume == null)
            Debug.LogWarning($"[RoomSmokeCSSceneSetup] ComputeShader not found: {PsiDensityComputePath}");
        if (display.particlesToDensityVolume == null)
            Debug.LogWarning($"[RoomSmokeCSSceneSetup] ComputeShader not found: {ParticlesToDensityComputePath}");
    }

    private static void AttachUnifiedCS(GameObject go)
    {
        var comp = go.AddComponent<SFUnifiedCS>();
        var so = new SerializedObject(comp);

        AssignComputeShader(so, "_kernelsShader", "SFComputeKernels");
        AssignComputeShader(so, "_fftShader", "SFComputeFFT");
        AssignComputeShader(so, "_particlesShader", "SFComputeParticles");
        AssignComputeShader(so, "_lesShader", "SFComputeLES");

        AssignComputeShaderByPath(so, "_containerIsfShader", ContainerIsfPath);
        AssignComputeShaderByPath(so, "_particlesWallShader", ParticlesWallPath);

        var scenarioProp = so.FindProperty("_scenario");
        if (scenarioProp != null)
            scenarioProp.enumValueIndex = (int)SFUnifiedCS.ScenarioType.RoomSmoke;

        so.ApplyModifiedPropertiesWithoutUndo();

        comp.ApplyRoomSmokeDefaults();
    }

    private static void AssignComputeShader(SerializedObject so, string propertyName, string assetName)
    {
        var guids = AssetDatabase.FindAssets($"{assetName} t:ComputeShader", new[] { ShaderDir });
        if (guids.Length == 0)
        {
            Debug.LogWarning($"[RoomSmokeCSSceneSetup] ComputeShader '{assetName}' not found in {ShaderDir}");
            return;
        }

        var path = AssetDatabase.GUIDToAssetPath(guids[0]);
        AssignComputeShaderByPath(so, propertyName, path);
    }

    private static void AssignComputeShaderByPath(SerializedObject so, string propertyName, string assetPath)
    {
        var shader = AssetDatabase.LoadAssetAtPath<ComputeShader>(assetPath);
        if (shader == null)
        {
            Debug.LogWarning($"[RoomSmokeCSSceneSetup] ComputeShader not found at {assetPath}");
            return;
        }

        var prop = so.FindProperty(propertyName);
        if (prop != null)
            prop.objectReferenceValue = shader;
    }
}
