using System.Diagnostics;
using System.IO;
using ShrodingerFlow.AnalyticalTests;
using UnityEditor;
using UnityEditor.SceneManagement;
using UnityEngine;
using Debug = UnityEngine.Debug;

public static class AnalyticalTestsSceneSetup
{
    private const string ShaderDir = "Assets/Scenes/JetExample/ComputeShader/";
    private const string ScenePath = "Assets/Tests/Analytical/Scenes/AnalyticalTests.unity";

    [MenuItem("ShrodingerFlow/Analytical Tests/Create or Update Scene")]
    public static void CreateScene()
    {
        var dir = Path.GetDirectoryName(ScenePath);
        if (!string.IsNullOrEmpty(dir) && !Directory.Exists(dir))
            Directory.CreateDirectory(dir);

        var scene = EditorSceneManager.NewScene(NewSceneSetup.DefaultGameObjects, NewSceneMode.Single);

        var go = new GameObject("AnalyticalTests");
        var runner = go.AddComponent<AnalyticalTestRunner>();
        go.AddComponent<VortexRingConservationTest>();
        go.AddComponent<IsfVsLesTest>();
        go.AddComponent<HbarSweepTest>();
        go.AddComponent<DensityTransportTest>();
        go.AddComponent<OpticalConsistencyTest>();
        go.AddComponent<PerformanceBenchmarkTest>();

        var so = new SerializedObject(runner);
        AssignComputeShader(so, "_kernelsShader", "SFComputeKernels");
        AssignComputeShader(so, "_fftShader", "SFComputeFFT");
        AssignComputeShader(so, "_lesShader", "SFComputeLES");
        so.ApplyModifiedPropertiesWithoutUndo();

        EditorSceneManager.SaveScene(scene, ScenePath);
        Debug.Log($"[AnalyticalTests] Сцена создана/обновлена: {ScenePath}");
        Debug.Log("[AnalyticalTests] Выделите GameObject 'AnalyticalTests' и нажмите Play, чтобы прогнать все включённые тесты.");
    }

    [MenuItem("ShrodingerFlow/Analytical Tests/Open Scene")]
    public static void OpenScene()
    {
        if (!File.Exists(ScenePath))
        {
            Debug.LogWarning($"[AnalyticalTests] Сцена не найдена: {ScenePath}. Сначала вызовите 'Create or Update Scene'.");
            return;
        }
        EditorSceneManager.OpenScene(ScenePath, OpenSceneMode.Single);
    }

    [MenuItem("ShrodingerFlow/Analytical Tests/Open Results Folder")]
    public static void OpenResultsFolder()
    {
        string repoRoot = Directory.GetParent(Application.dataPath).FullName;
        string resultsRoot = Path.Combine(repoRoot, "AnalyticalTestResults");
        if (!Directory.Exists(resultsRoot))
            Directory.CreateDirectory(resultsRoot);
        OpenInFileBrowser(resultsRoot);
    }

    private static void AssignComputeShader(SerializedObject so, string propertyName, string assetName)
    {
        var guids = AssetDatabase.FindAssets($"{assetName} t:ComputeShader", new[] { ShaderDir });
        if (guids.Length == 0)
        {
            Debug.LogWarning($"[AnalyticalTests] ComputeShader '{assetName}' не найден в {ShaderDir}");
            return;
        }
        var path = AssetDatabase.GUIDToAssetPath(guids[0]);
        var shader = AssetDatabase.LoadAssetAtPath<ComputeShader>(path);
        var prop = so.FindProperty(propertyName);
        if (prop != null)
            prop.objectReferenceValue = shader;
    }

    private static void OpenInFileBrowser(string path)
    {
#if UNITY_EDITOR_OSX
        Process.Start("open", $"\"{path}\"");
#elif UNITY_EDITOR_WIN
        Process.Start("explorer.exe", path.Replace('/', '\\'));
#else
        Process.Start("xdg-open", path);
#endif
    }
}
