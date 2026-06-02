using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.Universal;

namespace ShrodingerFlow.Particles
{
    /// <summary>
    /// Отрисовка: instanced billboard/сферы или полноэкранный реймарш по объёму плотности из ψ (<see cref="DisplayMode.Raymarch"/>).
    /// Ожидает <see cref="ParticleGpuBuffers"/> для частиц и компонент с <see cref="IRaymarchDensitySource"/> для реймарша.
    /// </summary>
    public class ParticleDisplay3D : MonoBehaviour
    {
        public enum DisplayMode
        {
            None,
            Shaded3D,
            Billboard,
            Raymarch
        }

        /// <summary>
        /// Чем заполняется 3D-текстура для реймарша: |ψ|² по сетке симуляции или сплаты по GPU-позициям частиц.
        /// </summary>
        public enum RaymarchDensitySource
        {
            [InspectorName("Вероятность |ψ|² (сетка симуляции)")]
            PsiProbabilityDensity,
            [InspectorName("Частицы (сплаты по позициям)")]
            ParticleSplats
        }

        [Header("Settings")]
        public DisplayMode mode = DisplayMode.Billboard;
        public float scale = 5f;
        [SerializeField, InspectorName("Автоматический размер"),
         Tooltip("Если включено — поле Scale выставляется из симуляции на этом же GameObject (поле «Размер частиц» × множитель ниже).")]
        private bool _automaticSimulationScale;
        [SerializeField, InspectorName("Множитель"),
         Tooltip("Только при автоматическом размере: Scale = размер частиц симуляции × множитель (по умолчанию как раньше — 50).")]
        private float _particleSizeToScaleMultiplier = 50f;
        public Gradient colourMap;
        public int gradientResolution = 64;
        public float velocityDisplayMax = 5f;
        [Header("References")]
        public ParticleGpuBuffers buffers;
        public Shader shaderShaded;
        public Shader shaderBillboard;

        [Header("Raymarch")]
        [Tooltip("Fluid/Raymarching + ψ→3D. Нужен Renderer Feature «Raymarch Fluid» на URP Renderer (меню ShrodingerFlow → Rendering).")]
        public Shader shaderRaymarch;
        public RaymarchDensitySource raymarchDensitySource = RaymarchDensitySource.ParticleSplats;
        [Tooltip("PsiDensityToVolume3D — |ψ₁|²+|ψ₂|² (режим «Вероятность |ψ|²»).")]
        public ComputeShader psiDensityToVolume;
        [Tooltip("ParticlesToDensityVolume — сплаты частиц в объём (режим «Частицы»).")]
        public ComputeShader particlesToDensityVolume;
        [SerializeField, Tooltip("Радиус капли частицы в МИРОВЫХ единицах (режим «Частицы»). Форма одинакова на всех сценах, не зависит от разрешения/анизотропии сетки. ~ размер частицы.")]
        float _splatRadiusWorld = 0.07f;
        [SerializeField, Tooltip("Вес одного сплата (чем больше — ярче объём). Доля нормализации считается автоматически.")]
        uint _splatterWeightFixed = 65000;
        [SerializeField, Tooltip("Сколько частиц в вокселе соответствует «полной» плотности (=1). Меньше — ярче/плотнее дым. НЕ зависит от общего числа частиц — поэтому яркость не плывёт при росте/убыли популяции.")]
        float _splatNormalizeParticlesPerVoxel = 8f;

        [SerializeField] Light _raymarchSunLight;

        [Header("Raymarch — вид (одна кривая плотности + свет)")]
        [SerializeField, Tooltip("Гамма на сырую плотность: ниже 1 — мягче края, выше 1 — контрастнее.")]
        float _raymarchDensityGamma = 0.55f;
        [SerializeField, Tooltip("Общая «толщина» объёма по освещению (начни с 6–14).")]
        float _raymarchOpticalDensity = 11f;
        [SerializeField, Tooltip("Поглощение по лучу (Beer–Lambert). Выше — темнее силуэт.")]
        float _raymarchAbsorption = 0.42f;
        [SerializeField, Range(0f, 2f), Tooltip("Равномерное рассеяние (подсветка изнутри).")]
        float _raymarchScatterAmbient = 0.24f;
        [SerializeField, Range(0f, 3f), Tooltip("Рассеяние к направлению солнца.")]
        float _raymarchScatterSun = 1.2f;
        [SerializeField, Tooltip("Резкость блика по солнцу.")]
        float _raymarchSunPhasePower = 2.2f;
        [SerializeField] Color _fluidAmbient = new Color(0.1f, 0.22f, 0.38f);
        [SerializeField] Color _fluidSunTint = new Color(0.52f, 0.78f, 1f);

        [Header("Raymarch — стены/препятствия")]
        [SerializeField, Tooltip("Рисовать стены лабиринта как твёрдые непрозрачные тела (берётся та же маска, что задаёт границу потока).")]
        bool _renderWalls = true;
        [SerializeField, Tooltip("Цвет стен.")]
        Color _wallColor = new Color(0.20f, 0.23f, 0.30f);
        [SerializeField, Range(0f, 1f), Tooltip("Базовая подсветка стен (доля цвета вне зависимости от направления на солнце).")]
        float _wallAmbient = 0.4f;

        [Header("Raymarch — ψ только")]
        [SerializeField, Tooltip("Вычитается из ρ только для режима «Вероятность |ψ|²».")]
        float _raymarchDensityOffset = 5e-8f;

        [Header("Raymarch — шаг луча")]
        [SerializeField] float _raymarchStepSize = 0.018f;
        [SerializeField, Tooltip("Отладка: весь экран сиреневый, если луч пересёк объём. Выключи — будет обычный реймарш по плотности.")]
        bool _raymarchDebugTintWhenRayHitsBounds;

        public enum RaymarchDebugShaderOutput
        {
            Normal = 0,
            [InspectorName("DEBUG: весь экран зелёный")]
            SolidGreen = 1,
            [InspectorName("DEBUG: весь экран пурпурный")]
            SolidMagenta = 2
        }

        [SerializeField, Tooltip("Проверка: если при SolidGreen/SolidMagenta экран не меняется — отрисовка идёт не этим материалом/не тем проходом URP.")]
        RaymarchDebugShaderOutput _raymarchDebugShaderOutput;

        Mesh _mesh;
        Material _mat;
        ComputeBuffer _argsBuffer;
        Texture2D _gradientTexture;
        DisplayMode _modeOld;
        bool _needsUpdate = true;

        Material _raymarchMat;
        RenderTexture _densityVolumeRt;
        IRaymarchDensitySource _psiDensitySource;
        IRaymarchSolidMaskSource _solidMaskSource;
        RenderTexture _wallVolumeRt;
        int _wallVolRx = -1;
        int _wallVolRy;
        int _wallVolRz;
        int _densityVolRx = -1;
        int _densityVolRy;
        int _densityVolRz;

        ComputeBuffer _splatterScratch;
        int _splatterScratchVoxels = -1;

        static Mesh _fullscreenTriangleShared;
        internal static Mesh SharedFullscreenTriangleMesh => CreateFullscreenTriangleMesh();
        static readonly int ColourMapId = Shader.PropertyToID("_ColourMap");
        static readonly int RayDebugHitBoundsId = Shader.PropertyToID("_RayDebugHitBounds");
        static readonly int RaymarchDebugForceOutputId = Shader.PropertyToID("_RaymarchDebugForceOutput");

        const string PsiDensityKernel = "PsiToDensity";

        public bool AutomaticSimulationScale => _automaticSimulationScale;

        /// <summary>
        /// При паузе Play Mode в редакторе <see cref="LateUpdate"/> не вызывается, но колбэки URP
        /// (<see cref="RenderPipelineManager.beginCameraRendering"/>) всё ещё идут при перерисовке Game View.
        /// </summary>
        bool _usesScriptableRenderPipeline;

        void Awake()
        {
            if (buffers == null)
                buffers = GetComponent<ParticleGpuBuffers>();
            _psiDensitySource = GetComponent<IRaymarchDensitySource>();
            _solidMaskSource = GetComponent<IRaymarchSolidMaskSource>();
        }

        void OnEnable()
        {
            RefreshPipelineUsage();
            if (_usesScriptableRenderPipeline)
                RenderPipelineManager.beginCameraRendering += OnBeginCameraRendering;
            RefreshRaymarchBridgeRegistration();
        }

        void OnDisable()
        {
            RaymarchFluidBridge.Unregister(this);
            if (_usesScriptableRenderPipeline)
                RenderPipelineManager.beginCameraRendering -= OnBeginCameraRendering;
        }

        void RefreshPipelineUsage()
        {
            _usesScriptableRenderPipeline = GraphicsSettings.renderPipelineAsset != null;
        }

        /// <summary>
        /// В URP <see cref="LateUpdate"/> не вызывает отрисовку; режим Raymarch собирает материал здесь, чтобы при переключении
        /// режима не было кадра без <see cref="_raymarchMat"/>.
        /// </summary>
        void Update()
        {
            if (!_usesScriptableRenderPipeline || mode != DisplayMode.Raymarch)
                return;
            UpdateSettings();
        }

        void LateUpdate()
        {
            if (_usesScriptableRenderPipeline)
                return;

            var cam = Camera.main;
            if (cam == null)
                return;

            if (mode == DisplayMode.Raymarch)
                IssueRaymarchFullscreen(cam);
            else
                IssueDrawMeshInstancedIndirect(cam);
        }

        void OnBeginCameraRendering(ScriptableRenderContext context, Camera camera)
        {
            if (!ShouldDrawParticlesForCamera(camera))
                return;

            if (mode == DisplayMode.Raymarch)
                return;

            IssueDrawMeshInstancedIndirect(camera);
        }

        /// <summary>
        /// Та же логика, что раньше только для билборда: Main Camera или Scene View в редакторе.
        /// </summary>
        static bool ShouldDrawParticlesForCamera(Camera camera)
        {
            if (camera == null)
                return false;

            Camera main = Camera.main;
            if (main == null)
                return true;
            if (camera == main)
                return true;
#if UNITY_EDITOR
            if (camera.cameraType == CameraType.SceneView)
                return true;
#endif
            return false;
        }

        void IssueDrawMeshInstancedIndirect(Camera cam)
        {
            if (cam == null)
                return;

            UpdateSettings();

            if (buffers == null || buffers.PositionBuffer == null || buffers.ActiveCount <= 0)
                return;

            if (mode != DisplayMode.None && _mesh != null && _mat != null && _argsBuffer != null)
            {
                var bounds = new Bounds(Vector3.zero, Vector3.one * 10000f);
                Graphics.DrawMeshInstancedIndirect(_mesh, 0, _mat, bounds, _argsBuffer, 0, null,
                    ShadowCastingMode.Off, false, gameObject.layer, cam);
            }
        }

        /// <summary>Сборка ψ→объём и uniform’ы; вызывай перед <see cref="IssueRaymarchFullscreen"/>.</summary>
        internal bool TryPrepareRaymarchPipeline(Camera cam)
        {
            UpdateSettings();

            if (_raymarchMat == null || shaderRaymarch == null)
                return false;

            if (_psiDensitySource == null)
                _psiDensitySource = GetComponent<IRaymarchDensitySource>();
            if (_psiDensitySource == null ||
                !_psiDensitySource.TryGetPsiVolume(out ComputeBuffer p1, out ComputeBuffer p2,
                    out Vector3 volumeMinWorld, out Vector3 volumeSizeWorld,
                    out int rx, out int ry, out int rz))
                return false;

            bool splats = raymarchDensitySource == RaymarchDensitySource.ParticleSplats;
            if (!splats)
                ReleaseSplatterScratch();

            EnsureDensityVolume(rx, ry, rz);

            if (splats)
            {
                if (particlesToDensityVolume == null || buffers == null || buffers.PositionBuffer == null ||
                    buffers.ActiveCount <= 0)
                    return false;

                DispatchParticleSplatterDensity(rx, ry, rz, volumeMinWorld, volumeSizeWorld);
            }
            else
            {
                if (psiDensityToVolume == null)
                    return false;

                int k = psiDensityToVolume.FindKernel(PsiDensityKernel);
                psiDensityToVolume.SetBuffer(k, "Psi1", p1);
                psiDensityToVolume.SetBuffer(k, "Psi2", p2);
                psiDensityToVolume.SetTexture(k, "DensityOut", _densityVolumeRt);
                psiDensityToVolume.SetInt("_ResX", rx);
                psiDensityToVolume.SetInt("_ResY", ry);
                psiDensityToVolume.SetInt("_ResZ", rz);
                int gx = (rx + 7) / 8;
                int gy = (ry + 7) / 8;
                int gz = (rz + 7) / 8;
                psiDensityToVolume.Dispatch(k, gx, gy, gz);
            }

            bool walls = TryBuildWallVolume(rx, ry, rz);

            ApplyRaymarchUniforms(cam, volumeMinWorld, volumeSizeWorld, splats, walls, rx, ry, rz);
            return true;
        }

        /// <summary>Маска стен (из источника симуляции) → 3D-текстура для раймарча. Стены статичны, но дёшево пересобрать каждый кадр.</summary>
        bool TryBuildWallVolume(int rx, int ry, int rz)
        {
            if (!_renderWalls || psiDensityToVolume == null)
                return false;
            if (_solidMaskSource == null)
                _solidMaskSource = GetComponent<IRaymarchSolidMaskSource>();
            if (_solidMaskSource == null
                || !_solidMaskSource.TryGetSolidMask(out ComputeBuffer mask, out int mrx, out int mry, out int mrz)
                || mask == null || mrx != rx || mry != ry || mrz != rz)
                return false;

            EnsureWallVolume(rx, ry, rz);

            int k = psiDensityToVolume.FindKernel("MaskToDensity");
            psiDensityToVolume.SetBuffer(k, "SolidMask", mask);
            psiDensityToVolume.SetTexture(k, "DensityOut", _wallVolumeRt);
            psiDensityToVolume.SetInt("_ResX", rx);
            psiDensityToVolume.SetInt("_ResY", ry);
            psiDensityToVolume.SetInt("_ResZ", rz);
            psiDensityToVolume.Dispatch(k, (rx + 7) / 8, (ry + 7) / 8, (rz + 7) / 8);
            return true;
        }

        void EnsureWallVolume(int rx, int ry, int rz)
        {
            if (_wallVolumeRt != null && _wallVolRx == rx && _wallVolRy == ry && _wallVolRz == rz)
                return;

            ReleaseWallVolume();
            _wallVolRx = rx;
            _wallVolRy = ry;
            _wallVolRz = rz;

            // RHalf: UAV-совместим и фильтруется линейно на Metal (R8 как random-write не гарантирован).
            var desc = new RenderTextureDescriptor(rx, ry, RenderTextureFormat.RHalf, 0)
            {
                dimension = TextureDimension.Tex3D,
                volumeDepth = rz,
                enableRandomWrite = true,
                msaaSamples = 1
            };
            _wallVolumeRt = new RenderTexture(desc)
            {
                name = "WallMaskVolume3D",
                filterMode = FilterMode.Bilinear,
                wrapMode = TextureWrapMode.Clamp
            };
            _wallVolumeRt.Create();
        }

        void ReleaseWallVolume()
        {
            if (_wallVolumeRt != null)
            {
                _wallVolumeRt.Release();
                Destroy(_wallVolumeRt);
                _wallVolumeRt = null;
            }
            _wallVolRx = -1;
        }

        void DispatchParticleSplatterDensity(int rx, int ry, int rz, Vector3 volumeMinWorld, Vector3 volumeSizeWorld)
        {
            int total = rx * ry * rz;
            EnsureSplatterScratch(total);

            ComputeShader cs = particlesToDensityVolume;
            int kClear = cs.FindKernel("ClearScratch");
            cs.SetBuffer(kClear, "DensityScratch", _splatterScratch);
            cs.SetInt("_ResX", rx);
            cs.SetInt("_ResY", ry);
            cs.SetInt("_ResZ", rz);
            cs.SetInt("_TotalVoxels", total);
            cs.Dispatch(kClear, Mathf.Max(1, (total + 63) / 64), 1, 1);

            int kSplat = cs.FindKernel("SplatterParticles");
            cs.SetBuffer(kSplat, "DensityScratch", _splatterScratch);
            cs.SetBuffer(kSplat, "ParticlePositions", buffers.PositionBuffer);
            cs.SetInt("_ResX", rx);
            cs.SetInt("_ResY", ry);
            cs.SetInt("_ResZ", rz);
            cs.SetInt("_TotalVoxels", total);
            cs.SetInt("_ParticleCount", buffers.ActiveCount);
            cs.SetVector("_VolumeMinWorld", volumeMinWorld);
            cs.SetVector("_VolumeSizeWorld", volumeSizeWorld);
            cs.SetFloat("_SplatRadiusWorld", _splatRadiusWorld);
            cs.SetInt("_SplatWeightFixed", (int)_splatterWeightFixed);

            // Нормировка на ЛОКАЛЬНУЮ концентрацию (частиц/воксель), а не на общее число частиц.
            // Старое `activeCount * weight * 0.22` делало дым тем тусклее, чем больше частиц всего
            // (накопление дыма → темнее вместо плотнее) и заставляло яркость «плыть» по ходу симуляции.
            float denom = Mathf.Max(400f, _splatterWeightFixed * Mathf.Max(0.25f, _splatNormalizeParticlesPerVoxel));
            cs.SetFloat("_ScratchDenom", denom);

            cs.Dispatch(kSplat, Mathf.Max(1, (buffers.ActiveCount + 255) / 256), 1, 1);

            int kToTex = cs.FindKernel("ScratchToDensity");
            cs.SetBuffer(kToTex, "DensityScratch", _splatterScratch);
            cs.SetTexture(kToTex, "DensityOut", _densityVolumeRt);
            cs.SetInt("_ResX", rx);
            cs.SetInt("_ResY", ry);
            cs.SetInt("_ResZ", rz);
            cs.SetInt("_TotalVoxels", total);
            cs.SetFloat("_ScratchDenom", denom);

            int gx = (rx + 7) / 8;
            int gy = (ry + 7) / 8;
            int gz = (rz + 7) / 8;
            cs.Dispatch(kToTex, gx, gy, gz);
        }

        void EnsureSplatterScratch(int totalVoxels)
        {
            if (_splatterScratch != null && _splatterScratchVoxels == totalVoxels)
                return;

            ReleaseSplatterScratch();
            _splatterScratch = new ComputeBuffer(Mathf.Max(1, totalVoxels), sizeof(uint));
            _splatterScratchVoxels = totalVoxels;
        }

        void ReleaseSplatterScratch()
        {
            if (_splatterScratch != null)
            {
                _splatterScratch.Release();
                _splatterScratch = null;
            }

            _splatterScratchVoxels = -1;
        }

        internal Material RaymarchMaterialInternal => _raymarchMat;

        void IssueRaymarchFullscreen(Camera cam)
        {
            if (!TryPrepareRaymarchPipeline(cam))
                return;

            var mesh = SharedFullscreenTriangleMesh;
            Graphics.DrawMesh(mesh, Matrix4x4.identity, _raymarchMat, gameObject.layer, cam, 0,
                null, ShadowCastingMode.Off, false);
        }

        void EnsureDensityVolume(int rx, int ry, int rz)
        {
            if (_densityVolumeRt != null && _densityVolRx == rx && _densityVolRy == ry && _densityVolRz == rz)
                return;

            ReleaseDensityVolume();

            _densityVolRx = rx;
            _densityVolRy = ry;
            _densityVolRz = rz;

            var desc = new RenderTextureDescriptor(rx, ry, RenderTextureFormat.ARGBFloat, 0)
            {
                dimension = TextureDimension.Tex3D,
                volumeDepth = rz,
                enableRandomWrite = true,
                msaaSamples = 1
            };
            // Bilinear (трилинейно для 3D) — мягкие капли. Серое поле раньше давала не фильтрация,
            // а нормировка на общее число частиц + точечный LOAD в шейдере; то и другое исправлено.
            _densityVolumeRt = new RenderTexture(desc) { name = "PsiDensityVolume3D", filterMode = FilterMode.Bilinear, wrapMode = TextureWrapMode.Clamp };
            _densityVolumeRt.Create();
        }

        void ReleaseDensityVolume()
        {
            ReleaseSplatterScratch();
            ReleaseWallVolume();

            if (_densityVolumeRt != null)
            {
                _densityVolumeRt.Release();
                Destroy(_densityVolumeRt);
                _densityVolumeRt = null;
            }

            _densityVolRx = -1;
        }

        void ApplyRaymarchUniforms(Camera cam, Vector3 volumeMinWorld, Vector3 volumeSizeWorld, bool particleSplats,
            bool walls, int resX, int resY, int resZ)
        {
            if (_raymarchMat.shader != shaderRaymarch && shaderRaymarch != null)
                _raymarchMat.shader = shaderRaymarch;

            Vector3 cp = cam.transform.position;
            _raymarchMat.SetVector("_RayWorldSpaceCameraPos", new Vector4(cp.x, cp.y, cp.z, 1f));

            float zFar = cam.farClipPlane;
            Vector3 bl = cam.ViewportToWorldPoint(new Vector3(0f, 0f, zFar));
            Vector3 br = cam.ViewportToWorldPoint(new Vector3(1f, 0f, zFar));
            Vector3 tl = cam.ViewportToWorldPoint(new Vector3(0f, 1f, zFar));
            Vector3 tr = cam.ViewportToWorldPoint(new Vector3(1f, 1f, zFar));
            _raymarchMat.SetVector("_RayViewport_BL", new Vector4(bl.x, bl.y, bl.z, 0f));
            _raymarchMat.SetVector("_RayViewport_BR", new Vector4(br.x, br.y, br.z, 0f));
            _raymarchMat.SetVector("_RayViewport_TL", new Vector4(tl.x, tl.y, tl.z, 0f));
            _raymarchMat.SetVector("_RayViewport_TR", new Vector4(tr.x, tr.y, tr.z, 0f));

            _raymarchMat.SetTexture("_DensityMap", _densityVolumeRt);
            // Обязательно разрешение текущего ψ/volume — иначе в шейдере dim=(1,1,1) и все точки читают один воксель → серый прямоугольник.
            _raymarchMat.SetVector("_DensityRes",
                new Vector4(Mathf.Max(1, resX), Mathf.Max(1, resY), Mathf.Max(1, resZ), 0f));
            _raymarchMat.SetVector("boundsSize", volumeSizeWorld);
            _raymarchMat.SetVector("volumeMin", volumeMinWorld);
            _raymarchMat.SetFloat("volumeValueOffset", particleSplats ? 0f : _raymarchDensityOffset);

            float optical = _raymarchOpticalDensity * (particleSplats ? 1f : 0.35f);
            _raymarchMat.SetFloat("_DensityGamma", _raymarchDensityGamma);
            _raymarchMat.SetFloat("_OpticalDensity", optical);
            _raymarchMat.SetFloat("_Absorption", _raymarchAbsorption);
            _raymarchMat.SetFloat("_ScatterAmbient", _raymarchScatterAmbient);
            _raymarchMat.SetFloat("_ScatterSun", _raymarchScatterSun);
            _raymarchMat.SetFloat("_SunPhasePower", _raymarchSunPhasePower);
            _raymarchMat.SetVector("_FluidAmbient", new Vector4(_fluidAmbient.r, _fluidAmbient.g, _fluidAmbient.b, 1f));
            _raymarchMat.SetVector("_FluidSunTint", new Vector4(_fluidSunTint.r, _fluidSunTint.g, _fluidSunTint.b, 1f));

            _raymarchMat.SetFloat("viewMarchStepSize", _raymarchStepSize);

            // Стены: текстуру биндим всегда (валидную), включаем флагом — иначе семплинг несвязанной 3D-текстуры даёт мусор/варнинг.
            _raymarchMat.SetTexture("_WallMap", walls && _wallVolumeRt != null ? _wallVolumeRt : _densityVolumeRt);
            _raymarchMat.SetFloat("_WallEnabled", walls && _wallVolumeRt != null ? 1f : 0f);
            _raymarchMat.SetVector("_WallColor", new Vector4(_wallColor.r, _wallColor.g, _wallColor.b, 1f));
            _raymarchMat.SetFloat("_WallAmbient", _wallAmbient);

            Vector3 sunDir = ResolveRaymarchSunDirection();
            _raymarchMat.SetVector("dirToSun", sunDir);
            _raymarchMat.SetFloat(RayDebugHitBoundsId, _raymarchDebugTintWhenRayHitsBounds ? 1f : 0f);
            _raymarchMat.SetFloat(RaymarchDebugForceOutputId, (float)_raymarchDebugShaderOutput);
        }

        Vector3 ResolveRaymarchSunDirection()
        {
            if (_raymarchSunLight != null)
                return (-_raymarchSunLight.transform.forward).normalized;
            Light sun = RenderSettings.sun;
            if (sun != null)
                return (-sun.transform.forward).normalized;
            return new Vector3(0.35f, -0.85f, 0.25f).normalized;
        }

        static Mesh CreateFullscreenTriangleMesh()
        {
            if (_fullscreenTriangleShared != null)
                return _fullscreenTriangleShared;

            _fullscreenTriangleShared = new Mesh { name = "RaymarchFullscreenTriangle" };
            _fullscreenTriangleShared.vertices = new[]
            {
                new Vector3(-1f, -1f, 0f),
                new Vector3(-1f, 3f, 0f),
                new Vector3(3f, -1f, 0f)
            };
            _fullscreenTriangleShared.uv = new[]
            {
                new Vector2(0f, 0f),
                new Vector2(0f, 2f),
                new Vector2(2f, 0f)
            };
            _fullscreenTriangleShared.triangles = new[] { 0, 1, 2 };
            _fullscreenTriangleShared.RecalculateBounds();
            _fullscreenTriangleShared.UploadMeshData(true);
            return _fullscreenTriangleShared;
        }

        void UpdateSettings()
        {
            if (buffers == null && mode != DisplayMode.Raymarch)
                return;

            if (_modeOld != mode)
            {
                _modeOld = mode;

                if (_mat != null)
                {
                    Destroy(_mat);
                    _mat = null;
                }

                if (_raymarchMat != null)
                {
                    Destroy(_raymarchMat);
                    _raymarchMat = null;
                }

                if (mode != DisplayMode.None)
                {
                    if (mode == DisplayMode.Raymarch)
                    {
                        _mesh = null;
                        if (shaderRaymarch != null)
                            _raymarchMat = new Material(shaderRaymarch);
                        _needsUpdate = true;
                    }
                    else
                    {
                        if (buffers != null)
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

            if (_argsBuffer != null && _mesh != null && mode != DisplayMode.None && mode != DisplayMode.Raymarch)
                IndirectArgsUtil.CreateOrUpdateArgsBuffer(ref _argsBuffer, _mesh, buffers.ActiveCount);

            RefreshRaymarchBridgeRegistration();
        }

        void RefreshRaymarchBridgeRegistration()
        {
            if (!_usesScriptableRenderPipeline)
                return;
            EnsureMainCameraOpaqueCopyForRaymarch();
            if (isActiveAndEnabled && mode == DisplayMode.Raymarch)
                RaymarchFluidBridge.Register(this);
            else
                RaymarchFluidBridge.Unregister(this);
        }

        /// <summary>
        /// Нужно для чтения _CameraOpaqueTexture в Fluid/Raymarching (фон кадра без аналитического неба).
        /// </summary>
        void EnsureMainCameraOpaqueCopyForRaymarch()
        {
            if (!isActiveAndEnabled || mode != DisplayMode.Raymarch)
                return;
            var cam = Camera.main;
            if (cam == null)
                return;
            var urp = cam.GetUniversalAdditionalCameraData();
            if (urp != null)
                urp.requiresColorTexture = true;
        }

        /// <summary>Выставляет <see cref="scale"/> из размера частиц симуляции (вызывается компонентом симуляции).</summary>
        public void ApplyAutomaticScaleFromSimulation(float simulationParticleSize)
        {
            if (!_automaticSimulationScale) return;
            scale = simulationParticleSize * _particleSizeToScaleMultiplier;
        }

        void TryApplyAutomaticScaleFromSimulation()
        {
            if (!_automaticSimulationScale) return;
            var src = GetComponent<ISimulationParticleSizeSource>();
            if (src != null)
                ApplyAutomaticScaleFromSimulation(src.SimulationParticleSize);
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
            TryApplyAutomaticScaleFromSimulation();
            RefreshPipelineUsage();
            if (_usesScriptableRenderPipeline && mode == DisplayMode.Raymarch)
            {
                EnsureMainCameraOpaqueCopyForRaymarch();
                RefreshRaymarchBridgeRegistration();
            }
        }

        void OnDestroy()
        {
            if (_argsBuffer != null)
            {
                _argsBuffer.Release();
                _argsBuffer = null;
            }

            ReleaseDensityVolume();

            if (_mat != null)
                Destroy(_mat);
            if (_raymarchMat != null)
                Destroy(_raymarchMat);
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
