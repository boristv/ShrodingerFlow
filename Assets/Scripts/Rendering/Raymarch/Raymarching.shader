Shader "Fluid/Raymarching"
{
    Properties
    {
        [HideInInspector] _DensityMap("Density 3D", 3D) = "" {}
    }

    SubShader
    {
        Tags
        {
            "RenderPipeline" = "UniversalPipeline"
            "Queue" = "Transparent+520"
            "RenderType" = "Transparent"
        }

        Pass
        {
            Name "PsiRaymarchVolume"
            Tags { "LightMode" = "UniversalForward" }

            ZWrite Off
            ZTest Always
            Cull Off
            Blend One Zero

            HLSLPROGRAM
            #pragma vertex Vert
            #pragma fragment Frag
            #pragma target 4.5
            #pragma exclude_renderers gles2 gles3

            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/ShaderVariablesFunctions.hlsl"
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/DeclareOpaqueTexture.hlsl"

            TEXTURE3D(_DensityMap);
            SAMPLER(sampler_DensityMap);

            TEXTURE3D(_WallMap);
            SAMPLER(sampler_WallMap);
            float _WallEnabled;
            float4 _WallColor;
            float _WallAmbient;

            float4 _RayViewport_BL;
            float4 _RayViewport_BR;
            float4 _RayViewport_TL;
            float4 _RayViewport_TR;
            float4 _RayWorldSpaceCameraPos;

            float3 boundsSize;
            float3 volumeMin;
            float volumeValueOffset;

            float _DensityGamma;
            float _OpticalDensity;
            float _Absorption;
            float _ScatterAmbient;
            float _ScatterSun;
            float _SunPhasePower;
            float3 _FluidAmbient;
            float3 _FluidSunTint;

            float viewMarchStepSize;
            float3 dirToSun;
            float _RayDebugHitBounds;
            float _RaymarchDebugForceOutput;

            // Разрешение 3D-текстуры плотности (для texel Load без линейной фильтрации).
            float4 _DensityRes;

            struct Attributes
            {
                float4 positionOS : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct Varyings
            {
                float4 positionCS : SV_POSITION;
                float2 uv : TEXCOORD0;
            };

            Varyings Vert(Attributes v)
            {
                Varyings o;
                o.positionCS = float4(v.positionOS.xy, 0.0, 1.0);
                o.uv = v.uv;
                return o;
            }

            float2 RayBoxDst(float3 boundsMin, float3 boundsMax, float3 rayOrigin, float3 rayDir)
            {
                float3 invRayDir = 1.0 / rayDir;
                float3 t0 = (boundsMin - rayOrigin) * invRayDir;
                float3 t1 = (boundsMax - rayOrigin) * invRayDir;
                float3 tmin = min(t0, t1);
                float3 tmax = max(t0, t1);

                float dstA = max(max(tmin.x, tmin.y), tmin.z);
                float dstB = min(tmax.x, min(tmax.y, tmax.z));

                float dstToBox = max(0, dstA);
                float dstInsideBox = max(0, dstB - dstToBox);
                return float2(dstToBox, dstInsideBox);
            }

            float SampleDensityVoxel(float3 uvw)
            {
                // Трилинейная фильтрация (sampler берёт filterMode текстуры = Bilinear) — мягкие
                // капли вместо вокселей-кубиков. Точечный LOAD давал «квадрат вокруг частицы».
                return SAMPLE_TEXTURE3D_LOD(_DensityMap, sampler_DensityMap, saturate(uvw), 0).r;
            }

            float SampleWall(float3 uvw)
            {
                return SAMPLE_TEXTURE3D_LOD(_WallMap, sampler_WallMap, saturate(uvw), 0).r;
            }

            // Нормаль поверхности стены — антиградиент маски (указывает наружу из твёрдого).
            float3 WallNormal(float3 uvw)
            {
                float3 e = 1.0 / float3(max(_DensityRes.x, 1), max(_DensityRes.y, 1), max(_DensityRes.z, 1));
                float gx = SampleWall(uvw + float3(e.x, 0, 0)) - SampleWall(uvw - float3(e.x, 0, 0));
                float gy = SampleWall(uvw + float3(0, e.y, 0)) - SampleWall(uvw - float3(0, e.y, 0));
                float gz = SampleWall(uvw + float3(0, 0, e.z)) - SampleWall(uvw - float3(0, 0, e.z));
                float3 g = float3(gx, gy, gz);
                float len = length(g);
                return len > 1e-5 ? -g / len : float3(0, 1, 0);
            }

            half3 ACESFilm(half3 x)
            {
                half a = 2.51h;
                half b = 0.03h;
                half c = 2.43h;
                half d = 0.59h;
                half e = 0.14h;
                return saturate((x * (a * x + b)) / (x * (c * x + d) + e));
            }

            float3 SampleSky(float3 dir)
            {
                const float3 colGround = float3(0.35, 0.3, 0.35) * 0.53;
                const float3 colSkyHorizon = float3(1, 1, 1);
                const float3 colSkyZenith = float3(0.08, 0.37, 0.73);

                float sun = pow(max(0, dot(dir, dirToSun)), 500) * 1;
                float skyGradientT = pow(smoothstep(0, 0.4, dir.y), 0.35);
                float groundToSkyT = smoothstep(-0.01, 0, dir.y);
                float3 skyGradient = lerp(colSkyHorizon, colSkyZenith, skyGradientT);

                return lerp(colGround, skyGradient, groundToSkyT) + sun * (groundToSkyT >= 1);
            }

            half4 Frag(Varyings i) : SV_Target
            {
                if (_RaymarchDebugForceOutput > 1.5)
                    return half4(1, 0, 1, 1);
                if (_RaymarchDebugForceOutput > 0.5)
                    return half4(0, 1, 0, 1);

                float2 sceneUv = GetNormalizedScreenSpaceUV(i.positionCS);
                float3 bgScene = SampleSceneColor(sceneUv);

                float2 uv = i.uv;
#if UNITY_UV_STARTS_AT_TOP
                uv.y = 1.0 - uv.y;
#endif
                float3 rayBottom = lerp(_RayViewport_BL.xyz, _RayViewport_BR.xyz, uv.x);
                float3 rayTop = lerp(_RayViewport_TL.xyz, _RayViewport_TR.xyz, uv.x);
                float3 worldOnFarPlane = lerp(rayBottom, rayTop, uv.y);
                float3 rayDir = normalize(worldOnFarPlane - _RayWorldSpaceCameraPos.xyz);
                float3 rayPos = _RayWorldSpaceCameraPos.xyz;

                float2 bd = RayBoxDst(volumeMin, volumeMin + boundsSize, rayPos, rayDir);
                if (bd.y <= 1e-6)
                    return half4(bgScene, 1);

                if (_RayDebugHitBounds > 0.5)
                    return half4(1, 0, 1, 1);

                float step = max(viewMarchStepSize, 1e-4);
                float marchLen = max(0, bd.y - step * 0.25);
                float distAlong = step * 0.5;

                float3 scattered = 0;
                float transmittance = 1.0;
                float maxRho = 0;

                bool wallHit = false;
                float3 wallShaded = 0;

                uint iter = 0;
                while (distAlong < marchLen && iter < 512)
                {
                    float3 p = rayPos + rayDir * (bd.x + distAlong);
                    float3 uvw = (p - volumeMin) / boundsSize;

                    if (all(uvw >= 0) && all(uvw <= 1))
                    {
                        // Стена — непрозрачная поверхность: всё за ней перекрыто, луч останавливается.
                        if (_WallEnabled > 0.5 && SampleWall(uvw) > 0.5)
                        {
                            float3 n = WallNormal(uvw);
                            float ndl = max(0.0, dot(n, dirToSun));
                            // Полусферический ambient: грани, смотрящие вверх (к небу), светлее нижних —
                            // теневые стороны не проваливаются в чёрный. Плюс ламбертов вклад солнца.
                            float hemi = 0.5 + 0.5 * n.y;
                            float3 ambient = _WallColor.rgb * _WallAmbient * (0.5 + 0.5 * hemi);
                            float3 diffuse = _WallColor.rgb * (1.0 - _WallAmbient) * ndl;
                            wallShaded = ambient + diffuse;
                            wallHit = true;
                            break;
                        }

                        float rhoRaw = SampleDensityVoxel(uvw) - volumeValueOffset;
                        rhoRaw = max(rhoRaw, 0);
                        float rho = pow(saturate(rhoRaw), max(_DensityGamma, 0.01));
                        maxRho = max(maxRho, rho);

                        // Ниже порога не считаем σ — иначе ∫ даёт серый туман по всему силуэту коробки.
                        if (rho < 0.000015)
                            rho = 0;

                        // Единый коэффициент экстинкции σ_t (без отдельного множителя поглощения —
                        // иначе эмиссия > поглощения, альбедо>1, среда «светится» сверх энергии).
                        float sigmaT = rho * max(_OpticalDensity, 0) * (0.5 + max(_Absorption, 0));

                        float sunScatter = pow(max(0.0, dot(dirToSun, -rayDir)), max(_SunPhasePower, 0.01));
                        float3 Li = _FluidAmbient * _ScatterAmbient + _FluidSunTint * (_ScatterSun * sunScatter);

                        // Энергосохраняющее однократное рассеяние: за шаг взаимодействует доля (1−Tstep)
                        // приходящего света. На краю σ_t→0 ⇒ (1−Tstep)→0 ⇒ нет «лишнего» свечения каймы.
                        float Tstep = exp(-sigmaT * step);
                        scattered += transmittance * (1.0 - Tstep) * Li;
                        transmittance *= Tstep;
                    }

                    distAlong += step;
                    iter++;
                }

                // Непрозрачность дыма вдоль луча.
                float alpha = saturate(1.0 - transmittance);

                // Нет ни дыма, ни стены — фон без изменений (совпадает с ранним выходом → нет каймы).
                if (!wallHit && alpha < 0.002)
                    return half4(bgScene, 1);

                // «Фон» для луча: стена (если попали) перекрывает сцену; дым перед ней её затеняет.
                float3 background = wallHit ? wallShaded : bgScene;

                // Собственный цвет дыма (накопленный свет, нормированный на α), тонмапим ТОЛЬКО его.
                // Альфа-композитинг lerp(фон, дым, α): на краю α→0 ⇒ ровно фон (нет ни серой каймы,
                // ни засветки), внутри ⇒ цвет дыма. Фон не складывается с дымом → не выбивает в белый.
                float3 smokeCol = alpha > 1e-4 ? scattered / alpha : float3(0, 0, 0);
                smokeCol = ACESFilm(smokeCol);
                float3 rgb = lerp(background, smokeCol, alpha);
                return half4(rgb, 1);
            }
            ENDHLSL
        }
    }
    Fallback Off
}
