Shader "Developer/Checkerboard World"
{
    Properties
    {
        _Color("Color", Color) = (1,1,1,1)
        _Glossiness ("Smoothness", Range(0,1)) = 1.0
        _Metallic ("Metallic", Range(0,1)) = 0.0
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 200

        CGPROGRAM
        // Physically based Standard lighting model, and enable shadows on all light types
        #pragma surface surf Standard fullforwardshadows

        // Use shader model 3.0 target, to get nicer looking lighting
        #pragma target 3.0

        struct Input
        {
            float2 uv_MainTex;
            float3 worldPos;
            float3 worldNormal;
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

        void surf (Input IN, inout SurfaceOutputStandard o)
        {
            float3 blendWeight = abs(IN.worldNormal);

            // Calculate checkerboard mask in the range of 0 (even tiles) to 1 (odd tiles).
            int checkerMask;

            if (blendWeight.x > blendWeight.y && blendWeight.x > blendWeight.z)
            {
                checkerMask = (abs(floor(IN.worldPos.y)) + abs(floor(IN.worldPos.z))) % 2.0f;
            }
            else if (blendWeight.y > blendWeight.z && blendWeight.y > blendWeight.x)
            {
                checkerMask = (abs(floor(IN.worldPos.x)) + abs(floor(IN.worldPos.z))) % 2.0f;
            }
            else
            {
                checkerMask = (abs(floor(IN.worldPos.x)) + abs(floor(IN.worldPos.y))) % 2.0f;
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
