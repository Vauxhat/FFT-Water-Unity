using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class FloatingPhysics : MonoBehaviour
{
    public FourierWaterCPU _waterSimulator;
    private Vector3 _worldPosition;

    // Start is called before the first frame update
    void Start()
    {
        // Initialise world position.
        _worldPosition = this.transform.position;
    }

    // Update is called once per frame
    void Update()
    {
        // Check if water simulator exists.
        if (_waterSimulator)
        {
            // Update object position, using world position and displacement vector.
            this.transform.position = _worldPosition + _waterSimulator.GetDisplacement(_worldPosition);

            // Sample surface normal relative to world position.
            Vector3 normal = _waterSimulator.GetSurfaceNormal(_worldPosition);

            // Calculate rotation vector and angle.
            Vector3 rotationVector = Vector3.Cross(normal, Vector3.up);
            float angle = Vector3.Angle(normal, Vector3.up);

            // Create a quaternion which aligns the object to the surface.
            Quaternion objectToSurface = Quaternion.AngleAxis(angle, rotationVector);

            // Update object rotation, apply dampening.
            this.transform.rotation = Quaternion.Slerp(Quaternion.Euler(Vector3.zero), objectToSurface, 0.2f);
        }
    }
}
