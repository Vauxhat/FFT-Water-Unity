Shader "Water/Cascade Ocean"
{
    Properties
    {
        // Input properties.
        _TesselationLevel("Tesselation Level", Float) = 1
        _MinTessellationDist("Min Tessellation Distance", Float) = 8
        _MaxTessellationDist("Max Tessellation Distance", Float) = 64

        _LOD0("Texture", 2D) = "black" {}
        _LOD1("Texture", 2D) = "black" {}
        _LOD2("Texture", 2D) = "black" {}
        _LOD3("Texture", 2D) = "black" {}
    }
        SubShader
    {
        Tags { "RenderType" = "Opaque" "PreviewType" = "Plane" }
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

        output.edges[0] = tessFactor;
        output.edges[1] = tessFactor;
        output.edges[2] = tessFactor;
        output.inside = tessFactor;

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
        float4 colour : COLOR;
    };

    sampler2D _LOD0;
    sampler2D _LOD1;
    sampler2D _LOD2;
    sampler2D _LOD3;

    [UNITY_domain("tri")]
    DS_CONSTANT_OUTPUT DomainShaderFunction(PATCH_TESSELATION_OUTPUT tessFactor, OutputPatch<VS_CONTROL_POINT_OUTPUT, 3> patch, float3 barycentricCoordinates : SV_DomainLocation)
    {
        DS_CONSTANT_OUTPUT output;

        output.position = patch[0].position * barycentricCoordinates.x + patch[1].position * barycentricCoordinates.y + patch[2].position * barycentricCoordinates.z;
        output.uv = patch[0].uv * barycentricCoordinates.x + patch[1].uv * barycentricCoordinates.y + patch[2].uv * barycentricCoordinates.z;

        float2 direction = output.position.xz;
        float distance = length(direction);

        float4 a = (1.0f - min((max((distance - 16.0f), 0.0f) / 16.0f), 1.0f)) * tex2Dlod(_LOD0, float4(direction / 16.0f, 0.0f, 0.0f));
        float4 b = (1.0f - min((max((distance - 128.0f), 0.0f) / 128.0f), 1.0f)) * tex2Dlod(_LOD1, float4(direction / 128.0f, 0.0f, 0.0f));
        float4 c = tex2Dlod(_LOD2, float4(direction / 1024.0f, 0.0f, 0.0f));

        float4 heightSample = a + b + c;

        output.position.x -= heightSample.x;
        output.position.y += heightSample.y;
        output.position.z -= heightSample.z;

        output.position = mul(UNITY_MATRIX_VP, output.position);

        output.colour = heightSample;

        return output;
    }

    fixed4 FragmentShaderFunction(DS_CONSTANT_OUTPUT input) : SV_Target
    {
        // sample the texture
        fixed4 color = tex2D(_LOD1, input.uv);

        return input.colour;
    }
    ENDCG
}
    }
}
