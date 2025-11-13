#ifndef CAMERA_H
#define CAMERA_H

#include <glm/glm.hpp>

class Camera
{
public:
    glm::vec3 position;
    glm::vec3 front;
    glm::vec3 up;

    Camera() : position(0.0f), front(0.0f, 0.0f, -1.0f), up(0.0f, 1.0f, 0.0f) {}
    glm::mat4 getViewMatrix() { return glm::mat4(1.0f); }
    glm::mat4 getProjectionMatrix() { return glm::mat4(1.0f); }
};

#endif // CAMERA_H
