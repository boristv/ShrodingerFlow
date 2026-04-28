using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.Universal;

namespace ShrodingerFlow.Particles
{
    /// <summary>
    /// Полноэкранный реймарш в цвет камеры URP. Должен быть добавлен через <see cref="RaymarchFluidRendererFeature"/>.
    /// </summary>
    internal sealed class RaymarchFluidPass : ScriptableRenderPass
    {
        readonly ProfilingSampler _profilingSampler = new ProfilingSampler("ShrodingerFlow Raymarch");

        internal RaymarchFluidPass(RenderPassEvent evt)
        {
            renderPassEvent = evt;
        }

        public override void OnCameraSetup(CommandBuffer cmd, ref RenderingData renderingData)
        {
            ConfigureTarget(renderingData.cameraData.renderer.cameraColorTargetHandle);
        }

        public override void Execute(ScriptableRenderContext context, ref RenderingData renderingData)
        {
            var display = RaymarchFluidBridge.Active;
            if (display == null || display.mode != ParticleDisplay3D.DisplayMode.Raymarch || !display.isActiveAndEnabled)
                return;

            Camera cam = renderingData.cameraData.camera;
            if (cam.cameraType != CameraType.Game && cam.cameraType != CameraType.SceneView)
                return;
            if (cam.cameraType == CameraType.Game && Camera.main != null && cam != Camera.main)
                return;

            if (!display.TryPrepareRaymarchPipeline(cam))
                return;

            Material mat = display.RaymarchMaterialInternal;
            Mesh mesh = ParticleDisplay3D.SharedFullscreenTriangleMesh;
            if (mat == null || mesh == null)
                return;

            CommandBuffer cmd = CommandBufferPool.Get(nameof(RaymarchFluidPass));
            using (new ProfilingScope(cmd, _profilingSampler))
            {
                cmd.DrawMesh(mesh, Matrix4x4.identity, mat, 0);
            }

            context.ExecuteCommandBuffer(cmd);
            CommandBufferPool.Release(cmd);
        }
    }
}
