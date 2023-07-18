using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
public class WaterMesh : MonoBehaviour
{
    public float _width = 128.0f;
    public float _height = 128.0f;
    public int _xDiv = 2;
    public int _yDiv = 2;

    private Vector3[] vertices;
    private Vector2[] uv;
    private Vector3[] normals;
    private Vector4[] tangents;
    private int[] triangles;

    // Start is called before the first frame update
    void Start()
    {
        /*float borderSize = 950;
        float distance = 1000;

        vertices =  new Vector3[8];

        vertices[0] = new Vector3(distance, 0, distance);
        vertices[1] = new Vector3(-distance, 0, distance);
        vertices[2] = new Vector3(distance, 0, -distance);
        vertices[3] = new Vector3(-distance, 0, -distance);
        vertices[4] = new Vector3(distance - borderSize, 0, distance - borderSize);
        vertices[5] = new Vector3(-distance + borderSize, 0, distance - borderSize);
        vertices[6] = new Vector3(distance - borderSize, 0, -distance + borderSize);
        vertices[7] = new Vector3(-distance + borderSize, 0, -distance + borderSize);

        triangles = new int[]
        {
            0, 4, 1,
            1, 4, 5,
            1, 5, 3,
            3, 5, 7,
            3, 7, 2,
            2, 7, 6,
            2, 6, 0,
            0, 6, 4
        };*/

        GenerateMesh(_width, _height, _xDiv, _yDiv);

        Mesh mesh = new Mesh();
        mesh.vertices = vertices;
        mesh.triangles = triangles;
        mesh.uv = uv;
        mesh.normals = normals;
        mesh.tangents = tangents;
        mesh.name = "Water Mesh";

        this.GetComponent<MeshFilter>().mesh = mesh;
    }

    // Update is called once per frame
    void Update()
    {

    }

    void GenerateMesh(float width, float height, int x_div, int y_div)
    {
        vertices = new Vector3[(x_div * 2 + 1) * (y_div * 2 + 1)];    // Triangle.
        uv = new Vector2[(x_div * 2 + 1) * (y_div * 2 + 1)];
        normals = new Vector3[(x_div * 2 + 1) * (y_div * 2 + 1)];
        tangents = new Vector4[(x_div * 2 + 1) * (y_div * 2 + 1)];

        //vertices = new Vector3[(x_div + 1) * (y_div + 1)];              // Quad.

        float texelWidth = width / (float)x_div;
        float texelHeight = height / (float)y_div;

        /*for (int y = 0; y < y_div + 1; y++)
        {
            float yPos = y * texelHeight - height * 0.5f;

            for (int x = 0; x < x_div + 1; x++)
            {
                float xPos = x * texelWidth - width * 0.5f;

                vertices[x + y * (x_div + 1)] = new Vector3(xPos, 0.0f, yPos);
            }
        }*/

        // Triangle vertex generator.
        for (int y = 0; y < y_div * 2 + 1; y++)
        {
            float yPos = y * (texelHeight * 0.5f) - height * 0.5f;
            float v = (float)y / (float)y_div * 0.5f;

            for (int x = 0; x < x_div * 2 + 1; x++)
            {
                float xPos = x * (texelWidth * 0.5f) - width * 0.5f;
                float u = (float)x / (float)x_div * 0.5f;

                vertices[x + y * (x_div * 2 + 1)] = new Vector3(xPos, 0.0f, yPos);
                uv[x + y * (x_div * 2 + 1)] = new Vector2(u, v);
                normals[x + y * (x_div * 2 + 1)] = new Vector3(0.0f, 1.0f, 0.0f);
                tangents[x + y * (x_div * 2 + 1)] = new Vector4(1.0f, 0.0f, 0.0f, 1.0f);
            }
        }

        // Quad vertex generator.
        //for (int y = 0; y < y_div + 1; y++)
        //{
        //    float yPos = y * texelHeight - height * 0.5f;

        //    for (int x = 0; x < x_div + 1; x++)
        //    {
        //        float xPos = x * texelWidth - width * 0.5f;

        //        vertices[x + y * (x_div + 1)] = new Vector3(xPos, 0.0f, yPos);
        //    }
        //}

        triangles = new int[x_div * y_div * 8 * 3];   // Triangle.
        //triangles = new int[x_div * y_div * 4];         // Quad.

        /*for (int i = 0; i < triangles.Length; i += 6)
        {
            int gridY = i / 6 / x_div;
            int gridX = (i / 6) - gridY;

            // Top triangle.
            triangles[i] = x + (y * x_div);
            triangles[i + 1] = x + ((y + 1) * x_div);
            triangles[i + 2] = x + 1 + (y * x_div);

            // Bottom triangle.
            triangles[i + 3] = x + 1 + (y * x_div);
            triangles[i + 4] = x + ((y + 1) * x_div);
            triangles[i + 5] = x + 1 + ((y + 1) * x_div);
        }*/

        /*for (int y = 0; y < y_div; y++)
        {
            for (int x = 0; x < x_div; x++)
            {
                // Top triangle.
                triangles[(x + (y * x_div)) * 6 + 0] = x + (y * (x_div + 1));
                triangles[(x + (y * x_div)) * 6 + 1] = x + ((y + 1) * (x_div + 1));
                triangles[(x + (y * x_div)) * 6 + 2] = x + 1 + (y * (x_div + 1));

                // Bottom triangle.
                triangles[(x + (y * x_div)) * 6 + 3] = x + 1 + (y * (x_div + 1));
                triangles[(x + (y * x_div)) * 6 + 4] = x + ((y + 1) * (x_div + 1));
                triangles[(x + (y * x_div)) * 6 + 5] = x + 1 + ((y + 1) * (x_div + 1));
            }
        }*/

        // Triangle topography generator.
        for (int y = 0; y < y_div; y++)
        {
            for (int x = 0; x < x_div; x++)
            {
                // Top left triangle.
                triangles[(x + (y * x_div)) * 24 + 0] = 2 * x + (2 * y * (x_div * 2 + 1));
                triangles[(x + (y * x_div)) * 24 + 1] = 2 * x + 1 + ((2 * y + 1) * (x_div * 2 + 1));
                triangles[(x + (y * x_div)) * 24 + 2] = 2 * x + 1 + (2 * y * (x_div * 2 + 1));

                // Top right triangle.
                triangles[(x + (y * x_div)) * 24 + 3] = 2 * x + 1 + (2 * y * (x_div * 2 + 1));
                triangles[(x + (y * x_div)) * 24 + 4] = 2 * x + 1 + ((2 * y + 1) * (x_div * 2 + 1));
                triangles[(x + (y * x_div)) * 24 + 5] = 2 * x + 2 + (2 * y * (x_div * 2 + 1));

                // Left top triangle.
                triangles[(x + (y * x_div)) * 24 + 6] = 2 * x + (2 * y * (x_div * 2 + 1));
                triangles[(x + (y * x_div)) * 24 + 7] = 2 * x + ((2 * y + 1) * (x_div * 2 + 1));
                triangles[(x + (y * x_div)) * 24 + 8] = 2 * x + 1 + ((2 * y + 1) * (x_div * 2 + 1));

                // Left bottom triangle.
                triangles[(x + (y * x_div)) * 24 + 9] = 2 * x + ((2 * y + 1) * (x_div * 2 + 1));
                triangles[(x + (y * x_div)) * 24 + 10] = 2 * x + ((2 * y + 2) * (x_div * 2 + 1));
                triangles[(x + (y * x_div)) * 24 + 11] = 2 * x + 1 + ((2 * y + 1) * (x_div * 2 + 1));


                // Top left triangle.
                triangles[(x + (y * x_div)) * 24 + 14] = 2 * x + ((2 * y + 2) * (x_div * 2 + 1));
                triangles[(x + (y * x_div)) * 24 + 13] = 2 * x + 1 + ((2 * y + 1) * (x_div * 2 + 1));
                triangles[(x + (y * x_div)) * 24 + 12] = 2 * x + 1 + ((2 * y + 2) * (x_div * 2 + 1));

                // Top right triangle.
                triangles[(x + (y * x_div)) * 24 + 17] = 2 * x + 1 + ((2 * y + 2) * (x_div * 2 + 1));
                triangles[(x + (y * x_div)) * 24 + 16] = 2 * x + 1 + ((2 * y + 1) * (x_div * 2 + 1));
                triangles[(x + (y * x_div)) * 24 + 15] = 2 * x + 2 + ((2 * y + 2) * (x_div * 2 + 1));

                // Left top triangle.
                triangles[(x + (y * x_div)) * 24 + 20] = 2 * x + 2 + (2 * y * (x_div * 2 + 1));
                triangles[(x + (y * x_div)) * 24 + 19] = 2 * x + 2 + ((2 * y + 1) * (x_div * 2 + 1));
                triangles[(x + (y * x_div)) * 24 + 18] = 2 * x + 1 + ((2 * y + 1) * (x_div * 2 + 1));

                // Left bottom triangle.
                triangles[(x + (y * x_div)) * 24 + 23] = 2 * x + 2 + ((2 * y + 1) * (x_div * 2 + 1));
                triangles[(x + (y * x_div)) * 24 + 22] = 2 * x + 2 + ((2 * y + 2) * (x_div * 2 + 1));
                triangles[(x + (y * x_div)) * 24 + 21] = 2 * x + 1 + ((2 * y + 1) * (x_div * 2 + 1));
            }
        }

        // Quad topography generator.
        //for (int y = 0; y < y_div; y++)
        //{
        //    for (int x = 0; x < x_div; x++)
        //    {
        //        triangles[(x + (y * x_div)) * 4 + 0] = x + (y * (x_div + 1));
        //        triangles[(x + (y * x_div)) * 4 + 1] = x + ((y + 1) * (x_div + 1));
        //        triangles[(x + (y * x_div)) * 4 + 2] = x + 1 + ((y + 1) * (x_div + 1));
        //        triangles[(x + (y * x_div)) * 4 + 3] = x + 1 + (y * (x_div + 1));
        //    }
        //}
    }
}
