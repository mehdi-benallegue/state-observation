#include <state-observation/sensors-simulation/algorithm/rotation-velocity.hpp>

namespace stateObservation
{
    namespace algorithm
    {
        Vector3 RotationVelocity::rotationVelocityMeasure(const Vector3 & angVelocityVector, const Matrix3 & orientation) const
        {
          return Vector3(orientation.transpose()*angVelocityVector);
        }

    }
}
