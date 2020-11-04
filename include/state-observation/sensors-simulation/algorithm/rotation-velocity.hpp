/**
 * \file      rotation-velocity.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief     The implementation of the algorithm of a rotation velocity sensor
 *
 * \details
 *
 *
 */

#ifndef SENSORALGORITHMSROTATIONVELOCITYHPP
#define SENSORALGORITHMSROTATIONVELOCITYHPP

#include <state-observation/api.h>
#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{
namespace algorithm
{
/**
 * \class  RotationVelocity
 * \brief Implements the gyrometer measurement algorithm
 *
 */

class STATE_OBSERVATION_DLLAPI RotationVelocity
{
public:
  /// virtual destructor
  virtual ~RotationVelocity() {}

  /// The angular velocity measurement in the local frame represented by the orientation Matrix
  Vector3 rotationVelocityMeasure(const Vector3 & angVelocityVector, const Matrix3 & orientation) const;

protected:
};
} // namespace algorithm

} // namespace stateObservation

#endif // SENSORALGORITHMSROTATIONVELOCITYHPP
