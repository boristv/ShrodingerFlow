using ShrodingerFlow.Particles;
using UnityEditor;
using UnityEditor.SceneManagement;
using UnityEngine;

public static class SmokeMaze2DSceneSetup
{
    private const string ShaderDir = "Assets/Scenes/JetExample/ComputeShader/";
    private const string ParticlesShaderDir = "Assets/Particles/";
    private const string RaymarchShaderPath = "Assets/Scripts/Rendering/Raymarch/Raymarching.shader";
    private const string PsiDensityComputePath = "Assets/Particles/PsiDensityToVolume3D.compute";
    private const string ParticlesToDensityComputePath = "Assets/Particles/ParticlesToDensityVolume.compute";
    private const string ScenePath = "Assets/Scenes/JetExample/ComputeShader/SmokeMaze2DCS.unity";

    [MenuItem("ShrodingerFlow/Create Smoke Maze 2D Scene")]
    public static void CreateSmokeMaze2DScene()
    {
        var scene = EditorSceneManager.NewScene(NewSceneSetup.DefaultGameObjects, NewSceneMode.Single);
        SetupCamera();

        var go = new GameObject("ISF_SmokeMaze2D");
        go.transform.position = Vector3.zero;

        ConfigureGpuParticles(go);
        var comp = AttachUnifiedCS(go);

        comp.ApplySmokeMaze2DDefaults();

        EditorSceneManager.SaveScene(scene, ScenePath);
        Debug.Log($"[SmokeMaze2DSceneSetup] Scene saved: {ScenePath}");
    }

    private static void SetupCamera()
    {
        var cam = Camera.main;
        if (cam == null) return;
        cam.transform.position = new Vector3(2.5f, 8f, 1.5f);
        cam.transform.rotation = Quaternion.Euler(90f, 0f, 0f);
        cam.clearFlags = CameraClearFlags.SolidColor;
        cam.backgroundColor = new Color(0.15f, 0.16f, 0.18f, 1f);
        cam.fieldOfView = 60f;
        cam.orthographic = true;
        cam.orthographicSize = 3.5f;
    }

    private static void ConfigureGpuParticles(GameObject go)
    {
        go.AddComponent<ParticleGpuBuffers>();

        var display = go.AddComponent<ParticleDisplay3D>();
        display.mode = ParticleDisplay3D.DisplayMode.Billboard;
        display.scale = 0.06f * 50f;

        display.shaderBillboard = AssetDatabase.LoadAssetAtPath<Shader>($"{ParticlesShaderDir}ParticleBillboard.shader");
        display.shaderShaded = AssetDatabase.LoadAssetAtPath<Shader>($"{ParticlesShaderDir}Particle3DSurf.shader");
        display.shaderRaymarch = AssetDatabase.LoadAssetAtPath<Shader>(RaymarchShaderPath);
        display.psiDensityToVolume = AssetDatabase.LoadAssetAtPath<ComputeShader>(PsiDensityComputePath);
        display.particlesToDensityVolume = AssetDatabase.LoadAssetAtPath<ComputeShader>(ParticlesToDensityComputePath);
    }

    private static SFUnifiedCS AttachUnifiedCS(GameObject go)
    {
        var comp = go.AddComponent<SFUnifiedCS>();
        var so = new SerializedObject(comp);

        AssignComputeShader(so, "_kernelsShader", "SFComputeKernels");
        AssignComputeShader(so, "_fftShader", "SFComputeFFT");
        AssignComputeShader(so, "_particlesShader", "SFComputeParticles");
        AssignComputeShader(so, "_lesShader", "SFComputeLES");

        so.ApplyModifiedPropertiesWithoutUndo();
        return comp;
    }

    private static void AssignComputeShader(SerializedObject so, string propertyName, string assetName)
    {
        var guids = AssetDatabase.FindAssets($"{assetName} t:ComputeShader", new[] { ShaderDir });
        if (guids.Length == 0)
        {
            Debug.LogWarning($"[SmokeMaze2DSceneSetup] ComputeShader '{assetName}' not found in {ShaderDir}");
            return;
        }

        var path = AssetDatabase.GUIDToAssetPath(guids[0]);
        var shader = AssetDatabase.LoadAssetAtPath<ComputeShader>(path);
        var prop = so.FindProperty(propertyName);
        if (prop != null)
            prop.objectReferenceValue = shader;
    }
}
