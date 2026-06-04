using UnityEngine;

namespace ShrodingerFlow.Particles
{
    /// <summary>
    /// Источник плотностно-фазового поля α (поле СОСТОЯНИЯ гибридной волново-плотностной модели ISF)
    /// для прямого объёмного рендеринга: один буфер float на ячейку сетки, layout idx = i*resY*resZ + j*resZ + k.
    /// В отличие от частиц-сплатов, плотность дыма переносится той же скоростью, что восстановлена из ψ (гл. 5.1).
    /// </summary>
    public interface IRaymarchScalarFieldSource
    {
        bool TryGetScalarField(out ComputeBuffer alpha,
            out Vector3 volumeMinWorld, out Vector3 volumeSizeWorld,
            out int resX, out int resY, out int resZ);
    }
}
