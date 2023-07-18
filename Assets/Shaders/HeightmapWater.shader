Shader "Water/Heightmap Water"
{
    Properties
    {
        // Tesselation variables.
        _TessellationLevel("Tesselation Level", Float) = 1
        _MinTessellationDist("Min Tessellation Distance", Float) = 8
        _MaxTessellationDist("Max Tessellation Distance", Float) = 64

        _WaterColour("Water Colour", Color) = (0.5, 0.5, 0.8, 1.0)

        // Texture variables.
        _Heightmap("Heightmap", 2D) = "black" {}
        _Normalmap("Normal Map", 2D) = "black" {}
    }
        SubShader
    {
        Tags { "RenderType" = "Transparent" }
        Cull Off
        LOD 100

        GrabPass
        {
            "_BackgroundTexture"
        }

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
                float3 normal : NORMAL;
                float3 tangent : TANGENT;
            };

            // Output control point.
            struct VS_CONTROL_POINT_OUTPUT
            {
                float4 position : INTERNALTESSPOS;
                float2 uv : TEXCOORD0;
                float3 normal : NORMAL;
                float3 tangent : TANGENT;
            };

            // Vertex shader function.
            VS_CONTROL_POINT_OUTPUT VertexShaderFunction(VS_CONSTANT_INPUT input)
            {
                VS_CONTROL_POINT_OUTPUT output;

                // Pass input data to hull shader.
                output.position = input.vertex;
                output.position = mul(unity_ObjectToWorld, output.position);
                output.uv = input.uv;

                output.normal = input.normal;
                output.tangent = input.tangent;

                // Return output.
                return output;
            }

            // Tessellation variables.
            float _MinTessellationDist;
            float _MaxTessellationDist;

            // Distance based tesselation function.
            float DistanceTessellation(float4 position)
            {
                // Determine vector from camera to target position, calculate distance.
                float3 direction = position.xyz - _WorldSpaceCameraPos.xyz;
                float distance = length(direction);

                // Calculate linear interpolation factor (0 to 1) using distance. Transform linear factor into cosine.
                float linearFactor = 1.0f - clamp((distance - _MinTessellationDist) / (_MaxTessellationDist - _MinTessellationDist), 0.0f, 1.0f);
                float cosFactor = (1.0f - cos(linearFactor * UNITY_PI)) * 0.5f;

                // Return normalise tesselation factor (from 0 to 1).
                return cosFactor;
            }

            // Return true if all points are below the fustrum plane.
            bool IsBelowPlane(float4 a, float4 b, float4 c, int plane)
            {
                return dot(unity_CameraWorldClipPlanes[plane], a) < 0 && dot(unity_CameraWorldClipPlanes[plane], b) < 0 && dot(unity_CameraWorldClipPlanes[plane], c) < 0;
            }

            // Fustrum culling.
            bool FustrumCulling(float4 a, float4 b, float4 c)
            {
                /*return dot(unity_CameraWorldClipPlanes[0], position) > 0 &&
                    dot(unity_CameraWorldClipPlanes[1], position) > 0 &&
                    dot(unity_CameraWorldClipPlanes[2], position) > 0 &&
                    dot(unity_CameraWorldClipPlanes[3], position) > 0 &&
                    dot(unity_CameraWorldClipPlanes[4], position) > 0 &&
                    dot(unity_CameraWorldClipPlanes[5], position) > 0;*/

                return  !IsBelowPlane(a, b, c, 0) && !IsBelowPlane(a, b, c, 0) && !IsBelowPlane(a, b, c, 0) && !IsBelowPlane(a, b, c, 0) && !IsBelowPlane(a, b, c, 0) && !IsBelowPlane(a, b, c, 0);
            }

            float GetPostProjectionSphereExtent(float3 Origin, float Diameter)
            {
                float4 ClipPos = mul(UNITY_MATRIX_VP, float4(Origin, 1.0));
                return abs((Diameter * unity_CameraProjection._m11) / ClipPos.w);
            }

            float CalculateTessellationFactor(float3 Control0, float3 Control1)
            {
                float diameter = distance(Control0, Control1);
                float3 midpoint = (Control0 + Control1) * 0.5f;
                return max(1.0f, 16.0f * GetPostProjectionSphereExtent(midpoint, diameter));
            }

            // Output patch constant data (tri).
            struct PATCH_TESSELATION_OUTPUT
            {
                float edges[3] : SV_TessFactor;
                float inside : SV_InsideTessFactor;
            };

            // Tessellation variables.
            float _TessellationLevel;

            // Patch Constant Function
            PATCH_TESSELATION_OUTPUT PatchFunction(InputPatch<VS_CONTROL_POINT_OUTPUT, 3> patch, uint id : SV_PrimitiveID)
            {
                PATCH_TESSELATION_OUTPUT output;

                float4 a = (patch[2].position + patch[1].position) * 0.5f;
                float4 b = (patch[2].position + patch[0].position) * 0.5f;
                float4 c = (patch[0].position + patch[1].position) * 0.5f;

                // Check if each point of the traingle is within the fustrum.
                //bool drawTriangle = FustrumCulling(patch[0].position) || FustrumCulling(patch[1].position) ||  FustrumCulling(patch[2].position);

                bool drawTriangle = FustrumCulling(patch[0].position, patch[1].position, patch[2].position);

                output.edges[0] = lerp(1.0f, _TessellationLevel, DistanceTessellation((patch[2].position + patch[1].position) * 0.5f)) * drawTriangle;
                output.edges[1] = lerp(1.0f, _TessellationLevel, DistanceTessellation((patch[2].position + patch[0].position) * 0.5f)) * drawTriangle;
                output.edges[2] = lerp(1.0f, _TessellationLevel, DistanceTessellation((patch[0].position + patch[1].position) * 0.5f)) * drawTriangle;
                //output.edges[0] = CalculateTessellationFactor(patch[2].position, patch[1].position) * drawTriangle;
                //output.edges[1] = CalculateTessellationFactor(patch[2].position, patch[0].position) * drawTriangle;
                //output.edges[2] = CalculateTessellationFactor(patch[0].position, patch[1].position) * drawTriangle;
                output.inside = (output.edges[0] + output.edges[1] + output.edges[2]) / 3.0f;

                return output;
            }

            // Output patch constant data (quad).
            struct QUAD_PATCH_TESSELATION_OUTPUT
            {
                float edges[4] : SV_TessFactor;
                float inside[2] : SV_InsideTessFactor;
            };

            // Patch constant function (quad).
            QUAD_PATCH_TESSELATION_OUTPUT QuadPatchFunction(InputPatch<VS_CONTROL_POINT_OUTPUT, 4> patch, uint id : SV_PrimitiveID)
            {
                QUAD_PATCH_TESSELATION_OUTPUT output;

                // Calculate edge tessellation.
                output.edges[0] = 2;
                output.edges[1] = 2;
                output.edges[2] = 2;
                output.edges[3] = 2;

                // Calculate inside tessellation.
                output.inside[0] = 2;
                output.inside[1] = 2;

                return output;
            }

            [UNITY_domain("tri")]
            [UNITY_partitioning("fractional_odd")]
            [UNITY_outputtopology("triangle_cw")]
            [UNITY_outputcontrolpoints(3)]
            [UNITY_patchconstantfunc("PatchFunction")]
            VS_CONTROL_POINT_OUTPUT HullShaderFunction(InputPatch<VS_CONTROL_POINT_OUTPUT, 3> patch, uint id : SV_OutputControlPointID)
            {
                VS_CONTROL_POINT_OUTPUT output;

                output = patch[id];

                return output;
            }

            [UNITY_domain("quad")]
            [UNITY_partitioning("fractional_odd")]
            [UNITY_outputtopology("triangle_cw")]
            [UNITY_outputcontrolpoints(4)]
            [UNITY_patchconstantfunc("QuadPatchFunction")]
            VS_CONTROL_POINT_OUTPUT QuadHullShaderFunction(InputPatch<VS_CONTROL_POINT_OUTPUT, 4> patch, uint id : SV_OutputControlPointID)
            {
                VS_CONTROL_POINT_OUTPUT output;

                output = patch[id];

                return output;
            }

            // Output domain shader data.
            struct DS_CONSTANT_OUTPUT
            {
                float4 position : SV_POSITION;
                float2 uv : TEXCOORD0;
                float3 normal : NORMAL0;
                float3 binormal : BINORMAL0;
                float3 tangent : TANGENT0;
                float3 screenNormal : NORMAL1;
                float3 viewDirection : NORMAL2;
            };

            // Colour variables.
            sampler2D _Heightmap;
            sampler2D _Normalmap;

            float3 WorldToTangent(float3 normal)
            {
                float3x3 tbn =
                {
                    float3(1, 0, 0), float3(0, 1, 0), float3(0, 0, 1)
                };

                // Calculate normal direction in view space.
                normal = mul(UNITY_MATRIX_V, normal);

                // Convert normal from world space to tangent space (using transpose of tbn).
                normal = mul(normal, transpose(tbn));

                return normal;
            }

            [UNITY_domain("tri")]
            DS_CONSTANT_OUTPUT DomainShaderFunction(PATCH_TESSELATION_OUTPUT tessFactor, OutputPatch<VS_CONTROL_POINT_OUTPUT, 3> patch, float3 barycentricCoordinates : SV_DomainLocation)
            {
                DS_CONSTANT_OUTPUT output;

                // Interpolate input patch based on barycentric coordinates.
                output.position = patch[0].position * barycentricCoordinates.x + patch[1].position * barycentricCoordinates.y + patch[2].position * barycentricCoordinates.z;
                output.uv = patch[0].uv * barycentricCoordinates.x + patch[1].uv * barycentricCoordinates.y + patch[2].uv * barycentricCoordinates.z;
                output.normal = patch[0].normal * barycentricCoordinates.x + patch[1].normal * barycentricCoordinates.y + patch[2].normal * barycentricCoordinates.z;
                output.tangent = patch[0].tangent * barycentricCoordinates.x + patch[1].tangent * barycentricCoordinates.y + patch[2].tangent * barycentricCoordinates.z;

                // Calculate binormal from tangent and cross product.
                output.binormal = cross(output.normal, output.tangent);

                // Sample dispacement from heightmap.
                float3 displacement = tex2Dlod(_Heightmap, float4(output.uv, 0.0f, 0.0f)).rgb;

                // Add heightmap sample to position, transform position to clip space.
                output.position.rgb += displacement;
                output.position = mul(output.position, UNITY_MATRIX_VP);

                // Sample normal from normal map, convert to local vector.
                float3 normalSample = tex2Dlod(_Normalmap, float4(output.uv, 0.0f, 0.0f)).rgb;
                normalSample = (normalSample * 2.0f) - float3(1.0f, 1.0f, 1.0f);

                // Create a TBN matrix. This will scale each colour (r, g, b) from the normal along the surface axis.
                float3x3 tbn =
                {
                    output.tangent, output.binormal, output.normal
                };

                // Convert sampled normal from tangent space to world space.
                output.normal = mul(normalSample, tbn);

                // Calculate screen space normal.
                output.screenNormal = WorldToTangent(output.normal);

                return output;
            }

            // Global textures.
            sampler2D _BackgroundTexture;
            sampler2D _CameraDepthTexture;

            // Refraction variables.
            half _RefractionIntensity;

            fixed4 _WaterColour;

            // Fragment shader function.
            fixed4 FragmentShaderFunction(DS_CONSTANT_OUTPUT input) : SV_Target
            {
                // Calculate current sample position of screen vertex.
                half2 samplePosition = input.position.xy / _ScreenParams.xy;

                // Calculate refraction direction, normalised to screen space.
                half refraction = (input.screenNormal.xy / _ScreenParams.xy) * _RefractionIntensity;

                // Calculate depth for water surface and background.
                half surfaceDepth = input.position.z;
                half backgroundDepth = tex2D(_CameraDepthTexture, samplePosition + refraction).r;

                // Calculate refraction mask based on depth.
                int refractMask = ceil(surfaceDepth - backgroundDepth);

                // Sample background colour using refraction.
                fixed4 backgroundColour = tex2D(_BackgroundTexture, samplePosition + refraction * refractMask);

                half3 lightDirection = normalize(half3(1.0f, -1.0f, 1.0f));

                fixed4 ambientLight = fixed4(1.0f, 1.0f, 1.0f, 1.0f);
                fixed4 finalColour = saturate(backgroundColour * _WaterColour);
                //fixed4 finalColour = fixed4(1.0f, 1.0f, 1.0f, 1.0f);

                // Return final colour.
                return finalColour;
            }
            ENDCG
        }
    }
}
