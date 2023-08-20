using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class FloatingPhysics : MonoBehaviour
{
    public OceanGenerator _waterSimulator;
    private Vector3 _worldPosition;

    private Vector3 _forward;
    private Vector3 _right;

    private Vector3 _bounds;

    // Start is called before the first frame update
    void Start()
    {
        // Initialise world position.
        _worldPosition = this.transform.position;

        // Initialise directional vectors.
        _forward = this.transform.forward;
        _right = this.transform.right;

        // Get object bounds (no rotation) from mesh filter, apply scale.
        _bounds = Vector3.Scale(this.GetComponent<MeshFilter>().mesh.bounds.size, this.transform.localScale);
    }

    // Update is called once per frame
    void Update()
    {
        // Check if water simulator exists.
        if (_waterSimulator)
        {
            // Sample displacements along the forward and right vectors of the object bounds.
            Vector3 a = _worldPosition + _forward * _bounds.z * 0.5f + _waterSimulator.GetDisplacement(_worldPosition + _forward * _bounds.z * 0.5f);
            Vector3 b = _worldPosition - _forward * _bounds.z * 0.5f + _waterSimulator.GetDisplacement(_worldPosition - _forward * _bounds.z * 0.5f);
            Vector3 c = _worldPosition + _right * _bounds.z * 0.5f + _waterSimulator.GetDisplacement(_worldPosition + _right * _bounds.z * 0.5f);
            Vector3 d = _worldPosition - _right * _bounds.z * 0.5f + _waterSimulator.GetDisplacement(_worldPosition - _right * _bounds.z * 0.5f);

            // Calculate forward, up, and right vectors.
            Vector3 forward = Vector3.Normalize(a - b);
            Vector3 right = Vector3.Normalize(c - d);
            Vector3 up = Vector3.Cross(forward, right);

            // Update world position using average displacement.
            this.transform.position = 0.25f * (a + b + c + d);

            // Update rotation based on axis vectors.
            this.transform.rotation = Quaternion.LookRotation(forward, up);
        }
    }
}
