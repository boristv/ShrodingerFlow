#if UNITY_EDITOR
using UnityEditor;
using UnityEngine;

/// <summary>Создаёт SFUnifiedScenarioPresets.asset — правьте в инспекторе вместо кода.</summary>
public static class SFUnifiedScenarioPresetsMenu
{
    [MenuItem("ShrodingerFlow/Create SF Unified Scenario Presets Asset")]
    public static void CreateOrOverwrite()
    {
        const string path = SFUnifiedScenarioPresets.DefaultAssetPath;
        if (System.IO.File.Exists(path))
        {
            if (!EditorUtility.DisplayDialog(
                    "SF Unified Scenario Presets",
                    $"Файл уже есть:\n{path}\nПерезаписать встроенными значениями?",
                    "Перезаписать",
                    "Отмена"))
                return;
        }

        var s = SFUnifiedScenarioPresets.CreateBuiltIn();
        AssetDatabase.CreateAsset(s, path);
        AssetDatabase.SaveAssets();
        EditorUtility.FocusProjectWindow();
        Selection.activeObject = s;
        Debug.Log($"[SFUnifiedScenarioPresets] Created {path}");
    }
}
#endif
