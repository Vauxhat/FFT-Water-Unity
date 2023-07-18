Shader "Developer/Checkerboard Local"
{
    Properties
    {
        _Color("Color", Color) = (1,1,1,1)
        _Glossiness("Smoothness", Range(0,1)) = 1.0
        _Metallic("Metallic", Range(0,1)) = 0.0
    }
        SubShader
    {
        Tags { "RenderType" = "Opaque" }
        LOD 200

        CGPROGRAM
        // Physically based Standard lighting model, and enable shadows on all light types
        #pragma surface surf Standard fullforwardshadows vertex:vert

        // Use shader model 3.0 target, to get nicer looking lighting
        #pragma target 3.0

        struct Input
        {
            float2 uv_MainTex;
            float3 worldPos;
            float3 worldNormal;
            float3 localPos;
        };

        half _Glossiness;
        half _Metallic;
        fixed4 _Color;

        // Add instancing support for this shader. You need to check 'Enable Instancing' on materials that use the shader.
        // See https://docs.unity3d.com/Manual/GPUInstancing.html for more information about instancing.
        // #pragma instancing_options assumeuniformscaling
        UNITY_INSTANCING_BUFFER_START(Props)
            // put more per-instance properties here
        UNITY_INSTANCING_BUFFER_END(Props)

        void vert(inout appdata_full v, out Input o)
        {
            UNITY_INITIALIZE_OUTPUT(Input, o);

            // Store model matrix.
            float4x4 world = unity_ObjectToWorld;

            // Extract scale components from model matrix.
            half x = length(float3(world._11, world._21, world._31));
            half y = length(float3(world._12, world._22, world._32));
            half z = length(float3(world._13, world._23, world._33));

            // Create final scale matrix.
            float4x4 scale = {
                x, 0, 0, 0,
                0, y, 0, 0,
                0, 0, z, 0,
                0, 0, 0, 1
            };

            // Apply scale to vertex, ignoring rotation and translation.
            o.localPos = mul(v.vertex, scale);
        }

        void surf(Input IN, inout SurfaceOutputStandard o)
        {
            float3 blendWeight = fixed3(0.2f, 0.2f, 0.3f);
            blendWeight = abs(IN.worldNormal);

            // Calculate checkerboard mask in the range of 0 (even tiles) to 1 (odd tiles).
            int checkerMask = (abs(floor(IN.worldPos.x)) + abs(floor(IN.worldPos.y)) + abs(floor(IN.worldPos.z))) % 2.0f;

            if (blendWeight.x > blendWeight.y && blendWeight.x > blendWeight.z)
            {
                checkerMask = (abs(floor(IN.localPos.y)) + abs(floor(IN.localPos.z))) % 2.0f;
            }
            else if (blendWeight.y > blendWeight.z && blendWeight.y > blendWeight.x)
            {
                checkerMask = (abs(floor(IN.localPos.x)) + abs(floor(IN.localPos.z))) % 2.0f;
            }
            else
            {
                checkerMask = (abs(floor(IN.localPos.x)) + abs(floor(IN.localPos.y))) % 2.0f;
            }

            //int maskA = (abs(floor(IN.worldPos.y)) + abs(floor(IN.worldPos.z))) % 2.0f;
            //int maskB = (abs(floor(IN.worldPos.x)) + abs(floor(IN.worldPos.z))) % 2.0f;
            //int maskC = (abs(floor(IN.worldPos.x)) + abs(floor(IN.worldPos.y))) % 2.0f;

            //checkerMask = (abs(IN.uv_MainTex.x + IN.uv_MainTex.y) * 2.0f) % 1.0f;

            // Update albedo based on checkerboard mask, apply colour.
            o.Albedo = lerp(fixed4(0.6f, 0.6f, 0.6f, 1.0f), fixed4(0.4f, 0.4f, 0.4f, 1.0f), checkerMask) * _Color;

            // Metallic and smoothness come from slider variables
            o.Metallic = _Metallic;
            o.Smoothness = lerp(0.6f, 0.15f, checkerMask);
            //o.Alpha = c.a;
        }
        ENDCG
    }
        FallBack "Diffuse"
}
