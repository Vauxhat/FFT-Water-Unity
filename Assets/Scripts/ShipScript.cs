using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ShipScript : MonoBehaviour
{
    public WaterGenerator waterGenerator;

    Vector3 velocity = Vector3.zero;
    Vector3 worldPosition = Vector3.zero;
    Vector3 rotation = Vector3.zero;

    Vector3 up = Vector3.zero;

    Vector3 waterPosition = Vector3.zero;
    Quaternion waterRotation;

    public float bobbing = 8;

    // Start is called before the first frame update
    void Start()
    {
        waterPosition = gameObject.transform.position;
    }

    // Update is called once per frame
    void Update()
    {
        if (waterGenerator)
        {
            float range = 4.0f;

            Vector3 front = worldPosition + this.gameObject.transform.forward * range;
            Vector3 back = worldPosition - this.gameObject.transform.forward * range;
            Vector3 left = worldPosition - this.gameObject.transform.right * range;
            Vector3 right = worldPosition + this.gameObject.transform.right * range;

            // Update water position, moving towards target position.
            //Vector3 targetPosition = 0.25f * (waterGenerator.GetWaveOffset(front) + waterGenerator.GetWaveOffset(back) + waterGenerator.GetWaveOffset(left) + waterGenerator.GetWaveOffset(right));
            Vector3 targetPosition = waterGenerator.GetWaveOffset(worldPosition, bobbing);
            //waterPosition = Vector3.MoveTowards(waterPosition, targetPosition, 1.0f * Time.deltaTime);
            waterPosition = targetPosition;



            //Vector3 direction = targetPosition - waterPosition;
            //direction.Normalize();

            //velocity = Vector3.ClampMagnitude(velocity + direction.normalized * 2.0f * Time.deltaTime, 1.5f);
            //waterPosition += velocity * Time.deltaTime;

            // Apply new transform to game object.
            //this.gameObject.transform.position = worldPosition + waterPosition;
            this.gameObject.transform.position = worldPosition + waterPosition;

            Vector3 surfaceNormal = waterGenerator.GetWaveNormal(worldPosition, bobbing);
            Vector3 rotationVector = Vector3.Cross(Vector3.up, surfaceNormal);
            float angle = Vector3.Angle(Vector3.up, surfaceNormal);

            // Create a quaternion which will align the upward vector to the water surface.
            Quaternion objectToSurface = Quaternion.AngleAxis(angle, rotationVector);

            // Calculate new rotation target using local rotation and surface rotation.
            Quaternion targetRotation = Quaternion.Euler(rotation) * objectToSurface;

            // Update water rotation, moving towards new rotation target.
            waterRotation = Quaternion.RotateTowards(waterRotation, targetRotation, 10.0f * Time.deltaTime);

            // Apply new transform to game object.
            this.gameObject.transform.rotation = targetRotation;
        }
    }
}
