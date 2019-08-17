#ifndef D6_GLCAMERA_HH
#define D6_GLCAMERA_HH

#include "gl/opengl.hh"
#include "utils/alg.hh"

#include <Eigen/Geometry>

namespace d6
{

struct Camera
{
    Eigen::Matrix4f viewMatrix = Eigen::Matrix4f::Identity();
    Eigen::Matrix4f projectionMatrix = Eigen::Matrix4f::Identity();

    Eigen::Vector3f target = Eigen::Vector3f::Zero();
    Eigen::Vector3f position = Eigen::Vector3f::Zero();
    Eigen::Vector3f up = Eigen::Vector3f(0,0,1);

    float fieldOfView = M_PI_2;
    float nearPlane = 0.0f;
    float farPlane = 0.0f;
    float aspectRatio = 1.0f;

    bool valid() const 
    {
        return !target.isApprox(position) && nearPlane > 0.0f && farPlane > nearPlane;
    }

    Eigen::Vector3f direction() {
        return (target - position).normalized();
    }

    void lookAt(const Eigen::Vector3f &cameraPosition, const Eigen::Vector3f &cameraTarget, const Eigen::Vector3f &upVector)
    {
        position = cameraPosition;
        target = cameraTarget;
        up = upVector;
        rebuildView();
    }

    void rebuildView()
    {
        Eigen::Matrix3f R;
        R.col(2) = (position - target).normalized();
        R.col(0) = up.cross(R.col(2)).normalized();
        R.col(1) = R.col(2).cross(R.col(0));
        viewMatrix.topLeftCorner<3, 3>() = R.transpose();
        viewMatrix.topRightCorner<3, 1>() = -R.transpose() * position;
        viewMatrix(3, 3) = 1.0f;
    } 

    void setPerspective(float fovY, float aspect, float near, float far)
    {
        fieldOfView  = fovY;
        nearPlane = near;
        farPlane = far;
        aspectRatio = aspect;

        rebuildPerspective();
    }

    void rebuildPerspective()
    {
        float far = farPlane;
        float near = nearPlane;

        float theta = fieldOfView * 0.5;
        float range = far - near;
        float invtan = 1. / tan(theta);

        projectionMatrix(0, 0) = invtan / aspectRatio;
        projectionMatrix(1, 1) = invtan;
        projectionMatrix(2, 2) = -(near + far) / range;
        projectionMatrix(3, 2) = -1;
        projectionMatrix(2, 3) = -2 * near * far / range;
        projectionMatrix(3, 3) = 0;
    }

    void apply() {

        #ifndef GL_CORE
        glMatrixMode(GL_PROJECTION);
        glLoadMatrixf(projectionMatrix.data());
        glMatrixMode(GL_MODELVIEW);
        glLoadMatrixf(viewMatrix.data());
        #endif
    }

    void zoom(float amount)
    {
        if (!valid())
            return;
        const auto &dir = direction();

        position = target + (position - target).dot(dir) * amount * dir;

        nearPlane *= amount;
        farPlane *= amount;

        rebuildView();
        rebuildPerspective();
    }

    void rotate(float yawAngle, float pitchAngle)
    {
        if (!valid())
            return;
        const auto &dir = direction();

        const Eigen::Vector3f pitchAxis = dir.cross(up).normalized();

        Eigen::AngleAxisf yaw(yawAngle, up);
        Eigen::AngleAxisf pitch(pitchAngle, pitchAxis);

        position = target + pitch * yaw * (position - target);       

        rebuildView();
    }
    
    void translate(float xAmount, float yAmount)
    {
        const float depth  = (target - position).norm();
        const float planeHeight = depth * std::tan(fieldOfView * 0.5);

        Eigen::Vector3f worldMotion =
            viewMatrix.topLeftCorner<3, 3>().transpose() *
            Eigen::Vector3f(xAmount * aspectRatio * planeHeight, yAmount * planeHeight, 0.0f);

        target -= worldMotion ;
        position -= worldMotion ;

        rebuildView();
    }

    void frame(const Eigen::Vector3f& box) {
        const float previousDepth = (position - target).norm();

        target = box * 0.5f;

        const float height = 0.75f * std::fabs( viewMatrix.row(1).head<3>().dot(box) );
        const float width =  0.75f * box.lpNorm<Eigen::Infinity>()  / aspectRatio;

        const float depth = std::max(height, width) / std::tan(fieldOfView * 0.5);
        position = target - (depth) * direction();

        nearPlane *= depth/previousDepth;
        farPlane *= depth/previousDepth;
        
        rebuildView();
        rebuildPerspective();
    }
};

}

#endif