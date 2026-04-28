namespace ShrodingerFlow.Particles
{
    /// <summary>
    /// Активный <see cref="ParticleDisplay3D"/> в режиме Raymarch — читает <see cref="RaymarchFluidPass"/>.
    /// </summary>
    internal static class RaymarchFluidBridge
    {
        static ParticleDisplay3D _active;

        internal static void Register(ParticleDisplay3D owner)
        {
            _active = owner;
        }

        internal static void Unregister(ParticleDisplay3D owner)
        {
            if (_active == owner)
                _active = null;
        }

        internal static ParticleDisplay3D Active => _active;
    }
}
