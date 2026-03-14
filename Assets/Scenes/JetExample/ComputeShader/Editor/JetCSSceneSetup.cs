using UnityEditor;
using UnityEditor.SceneManagement;
using UnityEngine;
using UnityEngine.SceneManagement;

public static class JetCSSceneSetup
{
    private const string ScenePath = "Assets/Scenes/JetExample/ComputeShader/JetExampleCS.unity";
    private const string ShaderDir = "Assets/Scenes/JetExample/ComputeShader/";
    private const string MatPath   = "Assets/Scenes/JetExample/jet_particle.mat";

    [MenuItem("ShrodingerFlow/Create JetExample CS Scene")]
    public static void CreateScene()
    {
        var scene = EditorSceneManager.NewScene(NewSceneSetup.DefaultGameObjects, NewSceneMode.Single);

        SetupCamera();
        var jetGo = CreateJetObject();
        ConfigureParticleSystem(jetGo);
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

    private static void ConfigureParticleSystem(GameObject go)
    {
        var ps = go.AddComponent<ParticleSystem>();

        var main = ps.main;
        main.maxParticles = 1000;
        main.startLifetime = 9999f;
        main.startSpeed = 0f;
        main.startSize = 0.05f;
        main.simulationSpace = ParticleSystemSimulationSpace.World;
        main.playOnAwake = false;
        main.loop = true;

        var emission = ps.emission;
        emission.enabled = false;

        var shape = ps.shape;
        shape.enabled = false;

        var velocityOverLifetime = ps.velocityOverLifetime;
        velocityOverLifetime.enabled = false;

        var renderer = go.GetComponent<ParticleSystemRenderer>();
        renderer.renderMode = ParticleSystemRenderMode.Billboard;

        var mat = AssetDatabase.LoadAssetAtPath<Material>(MatPath);
        if (mat != null)
            renderer.material = mat;
        else
            Debug.LogWarning($"[JetCSSceneSetup] Material not found: {MatPath}");
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
