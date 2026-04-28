using UnityEditor;
using UnityEngine;
using UnityEngine.Rendering.Universal;

/// <summary>
/// Добавляет Renderer Feature «Raymarch Fluid» — без него URP не выполняет полноэкранный реймарш (будут только билборд и т.д.).
/// </summary>
public static class RaymarchFluidRendererFeatureInstaller
{
    const string DefaultRendererPath = "Assets/Settings/New Universal Render Pipeline Asset_Renderer.asset";

    [MenuItem("ShrodingerFlow/Rendering/Add Raymarch Fluid Pass To URP Renderer")]
    public static void InstallToDefaultRenderer()
    {
        InstallToRendererAtPath(DefaultRendererPath);
    }

    public static void InstallToRendererAtPath(string rendererAssetPath)
    {
        var rendererData = AssetDatabase.LoadAssetAtPath<UniversalRendererData>(rendererAssetPath);
        if (rendererData == null)
        {
            Debug.LogError($"[RaymarchFluid] UniversalRendererData не найден: {rendererAssetPath}");
            return;
        }

        foreach (var f in rendererData.rendererFeatures)
        {
            if (f is ShrodingerFlow.Particles.RaymarchFluidRendererFeature)
            {
                Debug.Log("[RaymarchFluid] Renderer Feature уже добавлен.");
                return;
            }
        }

        var feature = ScriptableObject.CreateInstance<ShrodingerFlow.Particles.RaymarchFluidRendererFeature>();
        feature.name = "RaymarchFluid";
        AssetDatabase.AddObjectToAsset(feature, rendererData);

        var so = new SerializedObject(rendererData);
        SerializedProperty prop = so.FindProperty("m_RendererFeatures");
        if (prop == null || !prop.isArray)
        {
            Debug.LogError("[RaymarchFluid] Не удалось найти m_RendererFeatures — версия URP другая.");
            Object.DestroyImmediate(feature);
            return;
        }

        prop.arraySize++;
        prop.GetArrayElementAtIndex(prop.arraySize - 1).objectReferenceValue = feature;
        so.ApplyModifiedPropertiesWithoutUndo();
        EditorUtility.SetDirty(rendererData);
        AssetDatabase.SaveAssets();
        Debug.Log($"[RaymarchFluid] Renderer Feature добавлен в {rendererAssetPath}");
    }
}
