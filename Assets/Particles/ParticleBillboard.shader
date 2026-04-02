Shader "Fluid/ParticleBillboard"
{
    Properties
    {
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
            Name "ForwardUnlit"
            Tags { "LightMode" = "UniversalForward" }

            ZWrite On
            ZTest LEqual
            Cull Off

            HLSLPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma target 4.5

            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"

            StructuredBuffer<float3> Positions;
            StructuredBuffer<float3> Velocities;

            float scale;

            struct Attributes
            {
                float4 positionOS : POSITION;
                float2 uv : TEXCOORD0;
                float3 normalOS : NORMAL;
            };

            struct Varyings
            {
                float4 positionCS : SV_POSITION;
                float3 colour : COLOR0;
            };

            Varyings vert(Attributes v, uint instanceID : SV_InstanceID)
            {
                Varyings o;

                float3 centreWorld = Positions[instanceID];
                float3 objectVertPos = v.positionOS.xyz * scale * 2.0;
                float4 viewPos = mul(UNITY_MATRIX_V, float4(centreWorld, 1.0)) + float4(objectVertPos, 0.0);
                o.positionCS = mul(UNITY_MATRIX_P, viewPos);

                // Как раньше у ParticleSystem: цвет от направления скорости (см. SFJetCS / SFUnifiedCS).
                float3 vel = Velocities[instanceID];
                float len = length(vel);
                float3 dir = len > 1e-6 ? vel / len : float3(0, 0, 0);
                o.colour = abs(dir);
                return o;
            }

            half4 frag(Varyings i) : SV_Target
            {
                return half4(i.colour, 1.0);
            }
            ENDHLSL
        }
    }
    FallBack Off
}
