using ShrodingerFlow.Particles;
using UnityEditor;
using UnityEditor.SceneManagement;
using UnityEngine;
using UnityEngine.SceneManagement;

public static class JetCSSceneSetup
{
    private const string ScenePath = "Assets/Scenes/JetExample/ComputeShader/JetExampleCS.unity";
    private const string ShaderDir = "Assets/Scenes/JetExample/ComputeShader/";
    private const string ParticlesShaderDir = "Assets/Particles/";
    private const string RaymarchShaderPath = "Assets/Scripts/Rendering/Raymarch/Raymarching.shader";
    private const string PsiDensityComputePath = "Assets/Particles/PsiDensityToVolume3D.compute";
    private const string ParticlesToDensityComputePath = "Assets/Particles/ParticlesToDensityVolume.compute";

    [MenuItem("ShrodingerFlow/Create JetExample CS Scene")]
    public static void CreateScene()
    {
        var scene = EditorSceneManager.NewScene(NewSceneSetup.DefaultGameObjects, NewSceneMode.Single);

        SetupCamera();
        var jetGo = CreateJetObject();
        ConfigureGpuParticles(jetGo);
        AttachJetCS(jetGo);

        EditorSceneManager.SaveScene(scene, ScenePath);
        Debug.Log($"[JetCSSceneSetup] Scene saved: {ScenePath}");
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

    private static GameObject CreateJetObject()
    {
        var go = new GameObject("Particles_Jet_CS");
        go.transform.position = new Vector3(0f, 5.45f, 0f);
        return go;
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
            Debug.LogWarning($"[JetCSSceneSetup] Shader not found: {ParticlesShaderDir}ParticleBillboard.shader");
        if (surf == null)
            Debug.LogWarning($"[JetCSSceneSetup] Shader not found: {ParticlesShaderDir}Particle3DSurf.shader");
        if (display.shaderRaymarch == null)
            Debug.LogWarning($"[JetCSSceneSetup] Shader not found: {RaymarchShaderPath}");
        if (display.psiDensityToVolume == null)
            Debug.LogWarning($"[JetCSSceneSetup] ComputeShader not found: {PsiDensityComputePath}");
        if (display.particlesToDensityVolume == null)
            Debug.LogWarning($"[JetCSSceneSetup] ComputeShader not found: {ParticlesToDensityComputePath}");
    }

    private static void AttachJetCS(GameObject go)
    {
        var jet = go.AddComponent<SFJetCS>();

        var so = new SerializedObject(jet);

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
            Debug.LogWarning($"[JetCSSceneSetup] ComputeShader '{assetName}' not found in {ShaderDir}");
            return;
        }

        var path = AssetDatabase.GUIDToAssetPath(guids[0]);
        var shader = AssetDatabase.LoadAssetAtPath<ComputeShader>(path);

        var prop = so.FindProperty(propertyName);
        if (prop != null)
            prop.objectReferenceValue = shader;
    }
}
