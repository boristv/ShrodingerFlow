using UnityEngine.Rendering.Universal;

namespace ShrodingerFlow.Particles
{
    /// <summary>
    /// Добавь на Universal Renderer Data (меню ShrodingerFlow → Rendering или вручную).
    /// Без этого проход реймарша не выполняется — на экране будет только билборд и т.п.
    /// </summary>
    public sealed class RaymarchFluidRendererFeature : ScriptableRendererFeature
    {
        RaymarchFluidPass _pass;

        public override void Create()
        {
            // До постобработки: иначе цвет/тонмап могут «перекрывать» или искажать полноэкранный проход.
            _pass = new RaymarchFluidPass(RenderPassEvent.BeforeRenderingPostProcessing);
        }

        public override void AddRenderPasses(ScriptableRenderer renderer, ref RenderingData renderingData)
        {
            if (_pass == null)
                return;
            renderer.EnqueuePass(_pass);
        }
    }
}
