using ShrodingerFlow.Particles;
using UnityEditor;
using UnityEditor.SceneManagement;
using UnityEngine;

public static class UnifiedCSSceneSetup
{
    private const string ShaderDir = "Assets/Scenes/JetExample/ComputeShader/";
    private const string ParticlesShaderDir = "Assets/Particles/";
    private const string ScenePath = "Assets/Scenes/JetExample/ComputeShader/UnifiedCS.unity";

    [MenuItem("ShrodingerFlow/Create Unified CS Scene")]
    public static void CreateUnifiedScene()
    {
        var scene = EditorSceneManager.NewScene(NewSceneSetup.DefaultGameObjects, NewSceneMode.Single);

        SetupCamera();

        var go = new GameObject("ISF_Simulation");
        go.transform.position = new Vector3(0f, 5.45f, 0f);

        ConfigureGpuParticles(go);
        AttachUnifiedCS(go);

        EditorSceneManager.SaveScene(scene, ScenePath);
        Debug.Log($"[UnifiedCSSceneSetup] Scene saved: {ScenePath}");
        Debug.Log("Select the ISF_Simulation object and choose a Scenario in the Inspector.");
    }

    private static void SetupCamera()
    {
        var cam = Camera.main;
        if (cam == null) return;
        cam.transform.position = new Vector3(2.55f, 6.55f, -5f);
        cam.transform.rotation = Quaternion.identity;
        cam.clearFlags = CameraClearFlags.SolidColor;
        cam.backgroundColor = new Color(0.678f, 0.737f, 0.831f, 0f);
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

        if (bb == null)
            Debug.LogWarning($"[UnifiedCSSceneSetup] Shader not found: {ParticlesShaderDir}ParticleBillboard.shader");
        if (surf == null)
            Debug.LogWarning($"[UnifiedCSSceneSetup] Shader not found: {ParticlesShaderDir}Particle3DSurf.shader");
    }

    private static void AttachUnifiedCS(GameObject go)
    {
        var comp = go.AddComponent<SFUnifiedCS>();
        var so = new SerializedObject(comp);

        AssignComputeShader(so, "_kernelsShader", "SFComputeKernels");
        AssignComputeShader(so, "_fftShader", "SFComputeFFT");
        AssignComputeShader(so, "_particlesShader", "SFComputeParticles");
        AssignComputeShader(so, "_lesShader", "SFComputeLES");

        so.ApplyModifiedPropertiesWithoutUndo();
    }

    private static void AssignComputeShader(SerializedObject so, string propertyName, string assetName)
    {
        var guids = AssetDatabase.FindAssets($"{assetName} t:ComputeShader", new[] { ShaderDir });
        if (guids.Length == 0)
        {
            Debug.LogWarning($"[UnifiedCSSceneSetup] ComputeShader '{assetName}' not found in {ShaderDir}");
            return;
        }

        var path = AssetDatabase.GUIDToAssetPath(guids[0]);
        var shader = AssetDatabase.LoadAssetAtPath<ComputeShader>(path);

        var prop = so.FindProperty(propertyName);
        if (prop != null)
            prop.objectReferenceValue = shader;
    }
}
