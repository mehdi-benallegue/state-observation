#include <state-observation/sensors-simulation/algorithm/linear-acceleration.hpp>

namespace stateObservation
{
    namespace algorithm
    {
        Vector3 LinearAcceleration::accelerationMeasure(const Vector3 & acceleration, const Matrix3 & orientation) const
        {
          return  Vector3(orientation.transpose()*(acceleration + cst::gravity));
        }

    }
}
