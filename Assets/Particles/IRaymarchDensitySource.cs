using UnityEngine;

namespace ShrodingerFlow.Particles
{
    /// <summary>
    /// Источник ψ-буферов и размеров объёма для режима <see cref="ParticleDisplay3D.DisplayMode.Raymarch"/>.
    /// </summary>
    public interface IRaymarchDensitySource
    {
        bool TryGetPsiVolume(out ComputeBuffer psi1, out ComputeBuffer psi2,
            out Vector3 volumeMinWorld, out Vector3 volumeSizeWorld,
            out int resX, out int resY, out int resZ);
    }
}
