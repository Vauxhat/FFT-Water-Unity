using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PlayerScript : MonoBehaviour
{
    Vector3 position;
    Vector3 rotation;
    Vector3 velocity;

    // State flags.
    public bool isSprinting { get; private set; } = false;

    // Movement variables.
    float walkSpeed = 3.0f;
    float runSpeed = 6.0f;

    private Transform playerTransform;
    public Transform boatTransform;

    // Start is called before the first frame update
    void Start()
    {
        position = gameObject.transform.position;
        rotation = gameObject.transform.rotation.eulerAngles;
        velocity = Vector3.zero;
    }

    // Update is called once per frame
    void Update()
    {
        // Define look sensitivity.
        float sensitivity = 2.0f;

        // Rotate camera based on mosue movement.
        rotation.y += Input.GetAxis("Mouse X") * sensitivity;

        this.transform.rotation = Quaternion.Euler(rotation);

        float moveSpeed;

        // Update sprint flag.
        isSprinting = Input.GetKey(KeyCode.LeftShift);

        // Set move speed based on whether the player is sprinting.
        switch (isSprinting)
        {
            case true: moveSpeed = runSpeed; break;
            default: moveSpeed = walkSpeed; break;
        }

        Vector3 moveDirection = Vector3.zero;

        if (Input.GetKey(KeyCode.W))
        {
            //this.transform.position += this.transform.forward * moveSpeed * Time.deltaTime;
            moveDirection += this.transform.forward;
        }
        if (Input.GetKey(KeyCode.S))
        {
            //this.transform.position -= this.transform.forward * moveSpeed * Time.deltaTime;
            moveDirection -= this.transform.forward;
        }
        if (Input.GetKey(KeyCode.A))
        {
            //this.transform.position -= this.transform.right * moveSpeed * Time.deltaTime;
            moveDirection -= this.transform.right;
        }
        if (Input.GetKey(KeyCode.D))
        {
            //this.transform.position += this.transform.right * moveSpeed * Time.deltaTime;
            moveDirection += this.transform.right;
        }

        // Normalise move direction.
        moveDirection.Normalize();

        // Update position based on move direction.
        velocity = moveDirection * moveSpeed;
        position += velocity * Time.deltaTime;

        //ayerTransform.position = position;
        //ayerTransform.rotation = Quaternion.Euler(rotation);

        this.transform.position = boatTransform.localToWorldMatrix.MultiplyPoint(position);
        this.transform.rotation = boatTransform.rotation * Quaternion.Euler(rotation);
    }
}
