using UnityEngine;
using UnityEngine.Rendering;

namespace ShrodingerFlow.Particles
{
    /// <summary>
    /// Отрисовка частиц через Graphics.DrawMeshInstancedIndirect (billboard или сферы).
    /// Ожидает <see cref="ParticleGpuBuffers"/> с тем же GameObject или назначенный вручную.
    /// </summary>
    public class ParticleDisplay3D : MonoBehaviour
    {
        public enum DisplayMode
        {
            None,
            Shaded3D,
            Billboard
        }

        [Header("Settings")]
        public DisplayMode mode = DisplayMode.Billboard;
        public float scale = 5f;
        public Gradient colourMap;
        public int gradientResolution = 64;
        public float velocityDisplayMax = 5f;
        [Header("References")]
        public ParticleGpuBuffers buffers;
        public Shader shaderShaded;
        public Shader shaderBillboard;

        Mesh _mesh;
        Material _mat;
        ComputeBuffer _argsBuffer;
        Texture2D _gradientTexture;
        DisplayMode _modeOld;
        bool _needsUpdate = true;

        static readonly int ColourMapId = Shader.PropertyToID("_ColourMap");

        void Awake()
        {
            if (buffers == null)
                buffers = GetComponent<ParticleGpuBuffers>();
        }

        void LateUpdate()
        {
            if (buffers == null || buffers.PositionBuffer == null || buffers.ActiveCount <= 0)
                return;

            UpdateSettings();

            if (mode != DisplayMode.None && _mesh != null && _mat != null && _argsBuffer != null)
            {
                var bounds = new Bounds(Vector3.zero, Vector3.one * 10000f);
                var cam = Camera.main;
                Graphics.DrawMeshInstancedIndirect(_mesh, 0, _mat, bounds, _argsBuffer, 0, null,
                    ShadowCastingMode.Off, false, gameObject.layer, cam);
            }
        }

        void UpdateSettings()
        {
            if (buffers == null) return;

            if (_modeOld != mode)
            {
                _modeOld = mode;
                if (mode != DisplayMode.None)
                {
                    _mesh = mode == DisplayMode.Billboard
                        ? ParticleMeshUtil.CreateQuadMesh()
                        : ParticleMeshUtil.CreateSphereMesh();

                    IndirectArgsUtil.CreateOrUpdateArgsBuffer(ref _argsBuffer, _mesh, buffers.ActiveCount);

                    _mat = mode switch
                    {
                        DisplayMode.Shaded3D => shaderShaded != null ? new Material(shaderShaded) : null,
                        DisplayMode.Billboard => shaderBillboard != null ? new Material(shaderBillboard) : null,
                        _ => null
                    };

                    if (_mat != null)
                    {
                        _mat.SetBuffer("Positions", buffers.PositionBuffer);
                        _mat.SetBuffer("Velocities", buffers.VelocityBuffer);
                        _needsUpdate = true;
                    }
                }
            }

            if (_mat != null)
            {
                if (_needsUpdate)
                {
                    _needsUpdate = false;
                    ParticleDisplay3D.TextureFromGradient(ref _gradientTexture, gradientResolution, colourMap);
                    _mat.SetTexture(ColourMapId, _gradientTexture);
                }

                _mat.SetFloat("scale", scale * 0.01f);
                _mat.SetFloat("velocityMax", velocityDisplayMax);

                Vector3 s = transform.localScale;
                transform.localScale = Vector3.one;
                Matrix4x4 localToWorld = transform.localToWorldMatrix;
                transform.localScale = s;

                _mat.SetMatrix("localToWorld", localToWorld);
            }

            if (_argsBuffer != null && _mesh != null && mode != DisplayMode.None)
                IndirectArgsUtil.CreateOrUpdateArgsBuffer(ref _argsBuffer, _mesh, buffers.ActiveCount);
        }

        public static void TextureFromGradient(ref Texture2D texture, int width, Gradient gradient,
            FilterMode filterMode = FilterMode.Bilinear)
        {
            if (texture == null)
                texture = new Texture2D(width, 1);
            else if (texture.width != width)
                texture.Reinitialize(width, 1);

            if (gradient == null)
            {
                gradient = new Gradient();
                gradient.SetKeys(
                    new[] { new GradientColorKey(Color.black, 0f), new GradientColorKey(Color.white, 1f) },
                    new[] { new GradientAlphaKey(1f, 0f), new GradientAlphaKey(1f, 1f) });
            }

            texture.wrapMode = TextureWrapMode.Clamp;
            texture.filterMode = filterMode;

            var cols = new Color[width];
            float denom = Mathf.Max(1f, cols.Length - 1f);
            for (int i = 0; i < cols.Length; i++)
            {
                float t = i / denom;
                cols[i] = gradient.Evaluate(t);
            }

            texture.SetPixels(cols);
            texture.Apply();
        }

        void OnValidate()
        {
            _needsUpdate = true;
        }

        void OnDestroy()
        {
            if (_argsBuffer != null)
            {
                _argsBuffer.Release();
                _argsBuffer = null;
            }

            if (_mat != null)
                Destroy(_mat);
        }
    }

    internal static class IndirectArgsUtil
    {
        public static void CreateOrUpdateArgsBuffer(ref ComputeBuffer buffer, Mesh mesh, int instanceCount)
        {
            uint ic = (uint)Mathf.Max(0, instanceCount);
            var args = new uint[5]
            {
                mesh.GetIndexCount(0),
                ic,
                mesh.GetIndexStart(0),
                (uint)mesh.GetBaseVertex(0),
                0
            };

            if (buffer == null)
                buffer = new ComputeBuffer(1, sizeof(uint) * 5, ComputeBufferType.IndirectArguments);
            buffer.SetData(args);
        }
    }

    internal static class ParticleMeshUtil
    {
        public static Mesh CreateQuadMesh()
        {
            var m = new Mesh { name = "ParticleQuad" };
            m.vertices = new[]
            {
                new Vector3(-0.5f, -0.5f, 0f),
                new Vector3(0.5f, -0.5f, 0f),
                new Vector3(0.5f, 0.5f, 0f),
                new Vector3(-0.5f, 0.5f, 0f)
            };
            m.uv = new[]
            {
                new Vector2(0, 0),
                new Vector2(1, 0),
                new Vector2(1, 1),
                new Vector2(0, 1)
            };
            m.triangles = new[] { 0, 1, 2, 0, 2, 3 };
            m.RecalculateNormals();
            m.RecalculateBounds();
            return m;
        }

        public static Mesh CreateSphereMesh()
        {
            var temp = GameObject.CreatePrimitive(PrimitiveType.Sphere);
            var mf = temp.GetComponent<MeshFilter>();
            Mesh shared = mf.sharedMesh;
            var copy = Object.Instantiate(shared);
            copy.name = "ParticleSphere";
            Object.Destroy(temp);
            return copy;
        }
    }
}
