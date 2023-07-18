using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CameraController : MonoBehaviour
{
    public enum CameraMode
    {
        FirstPerson,
        ThirdPerson
    }

    [Min(0.0f)]  public float distance = 4.0f;

    public CameraMode cameraMode = CameraMode.FirstPerson;

    Vector3 rotation;
    Vector3 position;

    // Start is called before the first frame update
    void Start()
    {
        // Initialise global variables.
        rotation = this.transform.eulerAngles;
        position = this.transform.localPosition;
    }

    // Update is called once per frame
    void Update()
    {
        // Define look sensitivity.
        float sensitivity = 2.0f;

        // Update camera based on camera mode.
        switch (cameraMode)
        {
            case CameraMode.ThirdPerson:
                {
                    // Update distance to target.
                    distance -= Input.mouseScrollDelta.y * 2.0f;

                    // Check if mouse button is held down.
                    if (Input.GetMouseButton(0))
                    {
                        // Rotate camera based on mosue movement.
                        rotation.y += Input.GetAxis("Mouse X") * sensitivity;
                        rotation.x -= Input.GetAxis("Mouse Y") * sensitivity;

                        // Clamp rotation on the x axis.
                        rotation.x = Mathf.Clamp(rotation.x, -90.0f, 90.0f);

                        // Update rotation.
                        this.transform.localRotation = Quaternion.Euler(rotation);
                    }

                    // Set camera position to orbit around focal point.
                    this.transform.localPosition = position - this.transform.forward * distance;

                    break;
                }
            default:
                {
                    // Rotate camera based on mosue movement.
                    rotation.x -= Input.GetAxis("Mouse Y") * sensitivity;

                    // Clamp rotation on the x axis.
                    rotation.x = Mathf.Clamp(rotation.x, -90.0f, 90.0f);

                    // Update rotation.
                    this.transform.localRotation = Quaternion.Euler(rotation);

                    // Set camera position.
                    this.transform.localPosition = position;

                    break;
                }
        }
    }
}
