namespace ShrodingerFlow.Particles
{
    /// <summary>Поставщик поля «Размер частиц» симуляции для автоматической подстройки scale на <see cref="ParticleDisplay3D"/>.</summary>
    public interface ISimulationParticleSizeSource
    {
        float SimulationParticleSize { get; }
    }
}
