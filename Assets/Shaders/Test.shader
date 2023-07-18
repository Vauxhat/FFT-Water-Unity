Shader "Water/Test"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
        _TopColour ("Top Colour", Color) = (1, 1, 1, 1)
        _BottomColour("Bottom Colour", Color) = (0, 0, 0, 0)
        _GradientHeight ("Gradient Height", Float) = 1
        _RefractionIntensity ("Refraction Intensity", Float) = 20
    }
    SubShader
    {
        Tags { "RenderType"="Transparent" "PreviewType"="Plane"}
        LOD 100

        GrabPass
        {
            "_BackgroundTexture"
        }

        Pass
        {
            CGPROGRAM
            #pragma vertex vertexShader
            #pragma fragment frag
            #pragma multi_compile_fog

            #include "UnityCG.cginc"

            // Define constant variables.
            static const float PI = 3.14159265f;

             // Input struct for vertex shader.
            struct vertexInput
            {
                float4 vertex : POSITION;
                float3 normal : NORMAL;
                float2 uv : TEXCOORD0;
                UNITY_FOG_COORDS(1)
            };

            // Output struct for vertex shader.
            struct vertexOutput
            {
                float4 vertex : SV_POSITION;
                float3 normal : NORMAL0;
                float3 binormal : BINORMAL;
                float3 tangent : TANGENT;
                float2 uv : TEXCOORD1;
                float4 grabPos : TEXCOORD0;
                float3 screenNormal : COLOUR;
                float3 viewDirection : NORMAL1;
                float3 colour : COLOUR1;
            };

            sampler2D _MainTex;
            float4 _MainTex_ST;

            float4 _TopColour;
            float4 _BottomColour;
            float _GradientHeight;
            
            int _NumWaves;

            float _Steepness[128];
            float _Amplitude[128];
            float _DirX[128];
            float _DirY[128];
            float _Speed[128];
            float _Wavelength[128];
            float _Frequency[128];

            vertexOutput vertexShader(vertexInput input)
            {
                // Initialise output.
                vertexOutput output;

                // Calculate world space position of input vertex.
                float4 worldPos = mul(unity_ObjectToWorld, input.vertex);

                for (int i = 0; i < _NumWaves; i++)
                {
                    // Calculate frequency (w) and phase based on input values.
                    //float w = (2.0f * UNITY_PI) / _Wavelength[i];
                    float w = _Frequency[i];
                    float phase = _Speed[i] * w;
                    float WA = w * _Amplitude[i];
                    float Q = _Steepness[i] / (w * _Amplitude[i] * _NumWaves);

                    // Precalculate sine and cosine values.
                    float cosine = cos(w * dot(float2(_DirX[i], _DirY[i]), worldPos.xz) + phase * _Time.y);
                    float sine = sin(w * dot(float2(_DirX[i], _DirY[i]), worldPos.xz) + phase * _Time.y);

                    // Apply wave transformation to world position.
                    worldPos.x += Q * _Amplitude[i] * _DirX[i] * cosine;
                    worldPos.y += _Amplitude[i] * sine;
                    worldPos.z += Q * _Amplitude[i] * _DirY[i] * cosine;

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

                // Calculate vector from camera to world position.
                output.viewDirection = normalize(worldPos - _WorldSpaceCameraPos);

                // Calculate vertex colour based on heigh.
                half lerpFactor = (1.0f - cos(clamp((worldPos.y - input.vertex.y) / _GradientHeight + 1.0f, 0, 1) * UNITY_PI)) / 2.0f;
                output.colour = lerp(_BottomColour, _TopColour, lerpFactor);

                // Transform output vertex into clip space.
                output.vertex = mul(UNITY_MATRIX_VP, worldPos);

                // Normalise normal and binormal vectors.
                normalize(output.normal);
                normalize(output.binormal);

                // Calculate tangent based on normal and binormal.
                output.tangent = cross(output.normal, output.binormal);

                output.uv = TRANSFORM_TEX(input.uv, _MainTex);
                UNITY_TRANSFER_FOG(output, output.vertex);

                output.grabPos = ComputeGrabScreenPos(output.vertex);

                // Create a matrix using the local axis (up, right, and forward), this will be used to translate the normals from world space to tangent space and vice-versa.
                float3x3 tbn = 
                { 
                    // The x value will be scaled along the binormal, y on the normal, and z along the tangent.
                    //output.tangent, output.binormal, output.normal
                    float3(1, 0, 0), float3(0, 1, 0), float3(0, 0, 1)
                };

                // Calculate normal direction in view space.
                output.screenNormal = mul(UNITY_MATRIX_MV, output.normal);

                // Convert normal from world space to tangent space (using transpose of tbn).
                output.screenNormal = normalize(mul(transpose(tbn), output.screenNormal));

                // Return finished output.
                return output;
            }

            // Structure for hull shader.
            struct tessFactor
            {
                float edge[3] : SV_TessFactor;
                float inside : SV_InsideTessFactor;
            };

            struct domainOutput
            {

            };

            tessFactor patchFunction(InputPatch<vertexInput, 3> patch)
            {
                tessFactor output;
                output.edge[0] = 1;
                output.edge[1] = 1;
                output.edge[2] = 1;
                output.inside = 1;
                return output;
            }

            [UNITY_domain("tri")]
            [UNITY_outputcontrolpoints(3)]
            [UNITY_outputtopology("triangle_cw")]
            [UNITY_partitioning("integer")]
            [UNITY_patchconstantfunc("patchFunction")]
            vertexOutput hullShader(InputPatch<vertexOutput, 3> patch, uint patchID : SV_OutputControlPointID)
            {
                return patch[patchID];
            }

            [UNITY_domain("tri")]
            domainOutput domainShader(tessFactor factors, OutputPatch<vertexInput, 3> patch, float3 uv : SV_DomainLocation)
            {
                domainOutput output;



                return output;
            }

            float4 tesselationSahder()
            {
                return 4.f;
            }

            sampler2D _BackgroundTexture;
            sampler2D _CameraDepthTexture;

            half _RefractionIntensity;

            fixed4 frag(vertexOutput i) : SV_Target
            {
                

                // Calculate current sample position of screen vertex.
                half2 samplePosition = i.vertex.xy / _ScreenParams.xy;

                
                //half backgroundDepth = tex2D(_CameraDepthTexture, samplePosition).r;

                //half depthFactor = min((surfaceDepth - backgroundDepth) / 0.01f, 1.0f);
                //half refraction = (i.screenNormal.xz / _ScreenParams.xy) * refractionIntensity;

                //fixed4 backgroundColour = tex2D(_BackgroundTexture, samplePosition + refraction * refractionIntensity);
                //fixed4 fogColour = fixed4(0.2f, 0.2f, 0.8f, 1.0f);

                // Calculate refraction direction, normalised to screen space.
                half refraction = (i.screenNormal.xy / _ScreenParams.xy) * _RefractionIntensity;

                // Calculate depth for water surface and background.
                half surfaceDepth = i.vertex.z;
                half backgroundDepth = tex2D(_CameraDepthTexture, samplePosition + refraction).r;

                // Calculate refraction mask based on depth.
                int refractMask = ceil(surfaceDepth - backgroundDepth);

                // Sample background colour using refraction.
                fixed4 backgroundColour = tex2D(_BackgroundTexture, samplePosition + refraction * refractMask);

                half3 lightDirection = normalize(half3(1.0f, -1.0f, 1.0f));

                half4 lightColour = half4(0.05, 0.2, 0.3, 1.0) * dot(-lightDirection, i.normal);
                //float4 lightColour = float4(1, 1, 1, 1) * dot(-lightDirection, i.normal);

                half shininess = 20.0f;
                half specularIntensity = 0.01f;
                half specularFactor = pow(max(dot(reflect(-lightDirection, i.normal), i.viewDirection), 0.0f), shininess);
                half4 specularColour = half4(1.0, 1.0, 1.0, 1.0) * specularFactor * specularIntensity;

                //half4 bgcolor = tex2D(_BackgroundTexture, i.grabPos + float4(i.screenNormal, 0.f) * 0.1f);
                half4 bgcolor = tex2D(_BackgroundTexture, samplePosition);
                //half4 bgcolor = tex2D(_CameraDepthTexture, samplePosition);
                //bgcolor = bgcolor * half4(0.0f, 0.0f, 1.0f, 1.0f);
                //float3 colour = (i.screenNormal + float3(1.0, 1.0, 1.0)) * 0.5;
                float3 colour = float3(0.8f, 0.8f, 1.0f);
                //fixed4 fogColour = half4(0.2f, 0.2f, 0.8f, 1.0f);
                // sample the texture
                //fixed4 col = float4(colour, 1.0) * bgcolor;
               // fixed4 col = float4(i.vertex.z, i.vertex.z, i.vertex.z, 1.0) * bgcolor;
                //fixed4 finalColour = saturate(lightColour + specularColour); //lerp(backgroundColour* fogColour, fogColour, depthFactor);
                fixed4 ambientLight = fixed4(1.0f, 1.0f, 1.0f, 1.0f);
                fixed4 finalColour = saturate(fixed4(i.colour, 1.0f) * backgroundColour * (dot(-lightDirection, i.normal) * 0.5 + 0.5));
                //fixed4 finalColour = fixed4((i.screenNormal + float3(1.0, 1.0, 1.0)) * 0.5f, 1.0f);
                // apply fog
                UNITY_APPLY_FOG(i.fogCoord, col);
                //return col;
                return finalColour;
            }
            ENDCG
        }
    }
}
