// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

Shader "Water/Gerstner Tessellated (128)"
{
    Properties
    {
        // Tesselation variables.
        _TessellationLevel("Tesselation Level", Float) = 1
        _MinTessellationDist("Min Tessellation Distance", Float) = 8
        _MaxTessellationDist("Max Tessellation Distance", Float) = 64

        // Water variables.
        _TopColour("Top Colour", Color) = (1, 1, 1, 1)
        _BottomColour("Bottom Colour", Color) = (0, 0, 0, 0)
        _GradientHeight("Gradient Height", Float) = 1
        _RefractionIntensity("Refraction Intensity", Float) = 20

        _Skybox("Skybox", Cube) = "" {}
    }
    SubShader
    {
        Tags { "RenderType"="Transparent" }
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
            };

            // Output control point.
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
                float3 normal : NORMAL0;
                float3 binormal : BINORMAL0;
                float3 tangent : TANGENT0;
                float2 uv : TEXCOORD0;
                float3 colour : COLOUR0;
                float3 screenNormal : NORMAL1;
                float3 viewDirection : NORMAL2;
            };

            // Colour variables.
            float4 _TopColour;
            float4 _BottomColour;
            float _GradientHeight;

            // Wave variables.
            int _NumWaves;
            float _Steepness[128];
            float _Amplitude[128];
            float _DirX[128];
            float _DirY[128];
            float _Speed[128];
            float _Wavelength[128];
            float _Frequency[128];

            [UNITY_domain("tri")]
            DS_CONSTANT_OUTPUT DomainShaderFunction(PATCH_TESSELATION_OUTPUT tessFactor, OutputPatch<VS_CONTROL_POINT_OUTPUT, 3> patch, float3 barycentricCoordinates : SV_DomainLocation)
            {
                DS_CONSTANT_OUTPUT output;

                // Interpolate input patch based on barycentric coordinates.
                output.position = patch[0 ].position * barycentricCoordinates.x + patch[1].position * barycentricCoordinates.y + patch[2].position * barycentricCoordinates.z;
                output.uv = patch[0].uv * barycentricCoordinates.x + patch[1].uv * barycentricCoordinates.y + patch[2].uv * barycentricCoordinates.z;

                // Initialise wave offset.
                float4 waveOffset = float4(0.0f, 0.0f, 0.0f, 1.0f);

                // Loop through each wave.
                for (int i = 0; i < _NumWaves; i++)
                {
                    float w = _Frequency[i];
                    float phase = _Speed[i] * w;
                    float WA = w * _Amplitude[i];
                    float Q = _Steepness[i] / (w * _Amplitude[i] * _NumWaves);

                    // Precalculate sine and cosine values.
                    float cosine = cos(w * dot(float2(_DirX[i], _DirY[i]), output.position.xz) + phase * _Time.y);
                    float sine = sin(w * dot(float2(_DirX[i], _DirY[i]), output.position.xz) + phase * _Time.y);

                    // Add wave transformation to total offset.
                    waveOffset.x += Q * _Amplitude[i] * _DirX[i] * cosine;
                    waveOffset.y += _Amplitude[i] * sine;
                    waveOffset.z += Q * _Amplitude[i] * _DirY[i] * cosine;

                    // Calculate surface normal.
                    output.normal.x -= _DirX[i] * WA * cosine;
                    output.normal.y += Q * WA * sine;
                    output.normal.z -= _DirY[i] * WA * cosine;

                    // Calculate surface binormal.
                    output.binormal.x += Q * _DirX[i] * _DirX[i] * WA * sine;
                    output.binormal.y += _DirX[i] * WA * cosine;
                    output.binormal.z -= Q * _DirX[i] * _DirY[i] * WA * sine;
                }

                // Final step for normal and binormal calculation.
                output.normal.y = 1.0 - output.normal.y;
                output.binormal.x = 1.0 - output.binormal.x;

                // Apply wave offset to output position.
                output.position += waveOffset * DistanceTessellation(output.position);

                // Calculate vector from camera to vertex.
                output.viewDirection = normalize(output.position - _WorldSpaceCameraPos);

                // Normalise normal and binormal vectors.
                normalize(output.normal);
                normalize(output.binormal);

                // Interpolate normal based on distance.
                output.normal = lerp(float3(0.0f, 1.0f, 0.0f), output.normal, DistanceTessellation(output.position));
                output.binormal = lerp(float3(1.0f, 0.0f, 0.0f), output.binormal, DistanceTessellation(output.position));

                // Transform output position from world space to clip space.
                output.position = UnityObjectToClipPos(output.position);

                // Calculate tangent based on normal and binormal.
                output.tangent = cross(output.normal, output.binormal);

                // Calculate vertex colour based on height.
                half lerpFactor = (1.0f - cos(clamp(waveOffset.y / _GradientHeight + 1.0f, 0, 1) * UNITY_PI)) / 2.0f;
                output.colour = lerp(_BottomColour, _TopColour, lerpFactor);

                // Create a matrix using the local axis (up, right, and forward), this will be used to translate the normals from world space to tangent space and vice-versa.
                float3x3 tbn =
                {
                    // The x value will be scaled along the binormal, y on the normal, and z along the tangent.
                    float3(1, 0, 0), float3(0, 1, 0), float3(0, 0, 1)
                };

                // Calculate normal direction in view space.
                output.screenNormal = mul(UNITY_MATRIX_MV, output.normal);

                // Convert normal from world space to tangent space (using transpose of tbn).
                output.screenNormal = normalize(mul(transpose(tbn), output.screenNormal));

                return output;
            }

            // Global textures.
            sampler2D _BackgroundTexture;
            sampler2D _CameraDepthTexture;

            // Rerfraction variables.
            half _RefractionIntensity;

            samplerCUBE _Skybox;

            // Fragment shader function.
            fixed4 FragmentShaderFunction (DS_CONSTANT_OUTPUT input) : SV_Target
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

                // Calculated reflected light direction.
                float3 reflection = reflect(input.viewDirection, input.normal);
                fixed4 reflectionColour = texCUBE(_Skybox, reflection);

                fixed4 ambientLight = fixed4(1.0f, 1.0f, 1.0f, 1.0f);
                fixed4 finalColour = saturate(fixed4(input.colour, 1.0f) * backgroundColour * (dot(-lightDirection, input.normal) * 0.5 + 0.5) + reflectionColour * 0.1f);

                // Return final colour.
                return finalColour;
            }
            ENDCG
        }
    }
}
