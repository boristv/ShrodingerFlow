using UnityEngine;

namespace ShrodingerFlow.Particles
{
    /// <summary>
    /// Источник статической маски твёрдых препятствий (стен) для раймарча: тот же буфер, что у симуляции,
    /// поэтому стены в кадре идеально совпадают с границей потока. Layout idx = i*resY*resZ + j*resZ + k.
    /// </summary>
    public interface IRaymarchSolidMaskSource
    {
        /// <summary>true и непустой <paramref name="solidMask"/> (int, 1 = твёрдое), если маска доступна для текущей сцены.</summary>
        bool TryGetSolidMask(out ComputeBuffer solidMask, out int resX, out int resY, out int resZ);
    }
}
