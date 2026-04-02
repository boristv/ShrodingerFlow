Shader "Fluid/Particle3DSurf"
{
    Properties
    {
        [HideInInspector] _ColourMap("ColourMap", 2D) = "white" {}
    }
    SubShader
    {
        Tags
        {
            "RenderType" = "Opaque"
            "Queue" = "Geometry"
            "RenderPipeline" = "UniversalPipeline"
        }

        Pass
        {
            Name "ForwardLit"
            Tags { "LightMode" = "UniversalForward" }

            ZWrite On
            ZTest LEqual
            Cull Back

            HLSLPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma target 4.5

            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Lighting.hlsl"

            StructuredBuffer<float3> Positions;
            StructuredBuffer<float3> Velocities;

            float scale;

            struct Attributes
            {
                float4 positionOS : POSITION;
                float3 normalOS : NORMAL;
            };

            struct Varyings
            {
                float4 positionCS : SV_POSITION;
                float3 normalWS : TEXCOORD0;
                float3 albedo : TEXCOORD1;
            };

            Varyings vert(Attributes v, uint instanceID : SV_InstanceID)
            {
                Varyings o;

                float3 centre = Positions[instanceID];
                float3 worldPos = centre + v.positionOS.xyz * scale;

                float3 worldNormal = normalize(v.normalOS);

                float3 vel = Velocities[instanceID];
                float len = length(vel);
                float3 dir = len > 1e-6 ? vel / len : float3(0, 0, 0);
                float3 col = abs(dir);

                o.positionCS = TransformWorldToHClip(worldPos);
                o.normalWS = worldNormal;
                o.albedo = col;
                return o;
            }

            half4 frag(Varyings i) : SV_Target
            {
                float3 N = normalize(i.normalWS);
                Light mainLight = GetMainLight();

                float NdotL = saturate(dot(N, mainLight.direction));
                float3 ambient = i.albedo * 0.25;
                float3 diffuse = i.albedo * NdotL * mainLight.color;
                float3 rgb = ambient + diffuse;
                return half4(rgb, 1.0);
            }
            ENDHLSL
        }
    }
    FallBack Off
}
