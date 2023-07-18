Shader "Unlit/TessellationTest"
{
    Properties
    {
        // Input properties.
        _TesselationLevel("Tesselation Level", Float) = 1
        _MinTessellationDist("Min Tessellation Distance", Float) = 8
        _MaxTessellationDist("Max Tessellation Distance", Float) = 64
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100

        Pass
        {
            CGPROGRAM
            #pragma vertex VertexShaderFunction
            #pragma hull HullShaderFunction
            #pragma domain DomainShaderFunction
            #pragma fragment FragmentShaderFunction
            #pragma target 4.6

            #include "UnityCG.cginc"

            // Vertex shader input.
            struct VS_CONSTANT_INPUT
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            // Input control point.
            struct VS_CONTROL_POINT_OUTPUT
            {
                float4 position : INTERNALTESSPOS;
                float2 uv : TEXCOORD0;
            };

            // Vertex shader function.
            VS_CONTROL_POINT_OUTPUT VertexShaderFunction(VS_CONSTANT_INPUT input)
            {
                VS_CONTROL_POINT_OUTPUT output;

                // Pass input data to hull shader.
                output.position = input.vertex;
                output.position = mul(unity_ObjectToWorld, output.position);
                output.uv = input.uv;

                // Return output.
                return output;
            }

            float _MinTessellationDist;
            float _MaxTessellationDist;

            float DistanceTesselation(float4 position)
            {
                float3 direction = position.xyz - _WorldSpaceCameraPos;
                float distance = length(direction);

                float linearFactor = 1.0f - clamp((distance - _MinTessellationDist) / (_MaxTessellationDist - _MinTessellationDist), 0.0f, 1.0f);
                float cosFactor = (1.0f - cos(linearFactor * UNITY_PI)) * 0.5f;

                return lerp(0.0f, 1.0f, cosFactor);
            }

            // Output patch constant data.
            struct PATCH_TESSELATION_OUTPUT
            {
                float edges[3] : SV_TessFactor;
                float inside : SV_InsideTessFactor;
            };

            float _TesselationLevel;

            // Patch Constant Function
            PATCH_TESSELATION_OUTPUT PatchFunction(InputPatch<VS_CONTROL_POINT_OUTPUT, 3> patch, uint id : SV_PrimitiveID)
            {
                PATCH_TESSELATION_OUTPUT output;

                float tessFactor = _TesselationLevel;

                float lerpFactor = (1.0f - cos(DistanceTesselation((patch[2].position + patch[1].position) * 0.5f) * UNITY_PI)) / 2.0f;

                output.edges[0] = lerp(1.0f, tessFactor, DistanceTesselation((patch[2].position + patch[1].position) * 0.5f));
                output.edges[1] = lerp(1.0f, tessFactor, DistanceTesselation((patch[2].position + patch[0].position) * 0.5f));
                output.edges[2] = lerp(1.0f, tessFactor, DistanceTesselation((patch[0].position + patch[1].position) * 0.5f));
                output.inside = lerp(1.0f, tessFactor, DistanceTesselation((patch[0].position + patch[1].position + patch[2].position) / 3.0f));

                return output;
            }

            [UNITY_domain("tri")]
            [UNITY_partitioning("pow2")]
            [UNITY_outputtopology("triangle_cw")]
            [UNITY_outputcontrolpoints(3)]
            [UNITY_patchconstantfunc("PatchFunction")]
            VS_CONTROL_POINT_OUTPUT HullShaderFunction(InputPatch<VS_CONTROL_POINT_OUTPUT, 3> patch, uint id : SV_OutputControlPointID)
            {
                VS_CONTROL_POINT_OUTPUT output;

                output = patch[id];

                return output;
            }

            // Output domain shader data.
            struct DS_CONSTANT_OUTPUT
            {
                float2 uv : TEXCOORD0;
                float4 position : SV_POSITION;
            };

            [UNITY_domain("tri")]
            DS_CONSTANT_OUTPUT DomainShaderFunction(PATCH_TESSELATION_OUTPUT tessFactor, OutputPatch<VS_CONTROL_POINT_OUTPUT, 3> patch, float3 barycentricCoordinates : SV_DomainLocation)
            {
                DS_CONSTANT_OUTPUT output;

                output.position = patch[0].position * barycentricCoordinates.x + patch[1].position * barycentricCoordinates.y + patch[2].position * barycentricCoordinates.z;
                output.uv = patch[0].uv * barycentricCoordinates.x + patch[1].uv * barycentricCoordinates.y + patch[2].uv * barycentricCoordinates.z;

                output.position = mul(UNITY_MATRIX_VP, output.position);

                return output;
            }

            fixed4 FragmentShaderFunction (DS_CONSTANT_OUTPUT input) : SV_Target
            {
                // sample the texture
                fixed4 color = fixed4(1.0f, 1.0f, 1.0f, 1.0f);

                return color;
            }
            ENDCG
        }
    }
}
