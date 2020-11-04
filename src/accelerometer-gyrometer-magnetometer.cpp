#include <state-observation/sensors-simulation/accelerometer-gyrometer-magnetometer.hpp>

namespace stateObservation
{
AccelerometerGyrometerMagnetometer::AccelerometerGyrometerMagnetometer()
: r_(Matrix3::Zero()), acc_(Vector3::Zero()), omega_(Vector3::Zero()), magne_(Vector3::Zero()),
  output_(Vector::Zero(measurementSize_, 1)), currentStateSize_(stateSize_)
{
#ifdef STATEOBSERVATION_VERBOUS_CONSTRUCTORS
  std::cout << std::endl << "AccelerometerGyrometerMagnetometer Constructor" << std::endl;
#endif // STATEOBSERVATION_VERBOUS_CONSTRUCTOR

  matrixMode_ = false;

  // This is the magnetic field in Toulouse, in uT.
  magne_[0] = 0;
  magne_[1] = 32;
  magne_[2] = -37;
}

Index AccelerometerGyrometerMagnetometer::getStateSize_() const
{
  return currentStateSize_;
}

Index AccelerometerGyrometerMagnetometer::getMeasurementSize_() const
{
  return measurementSize_;
}

Vector AccelerometerGyrometerMagnetometer::computeNoiselessMeasurement_()
{
  if(!matrixMode_)
  {
    Quaternion q(state_.head<4>());
    r_ = q.toRotationMatrix();

    acc_ = state_.segment(4, 3);
    omega_ = state_.tail(3);
  }
  else
  {
    r_ = Eigen::Map<Matrix3>(&state_[0]);
    acc_ = state_.segment<3>(9);
    omega_ = state_.tail<3>();
  }

  output_.head<3>() = accelerationMeasure(acc_, r_);
  output_.segment<3>(3) = rotationVelocityMeasure(omega_, r_);
  output_.tail<3>() = magneticFieldMeasure(r_);

  return output_;
}

void AccelerometerGyrometerMagnetometer::setMatrixMode(bool matrixMode)
{
  matrixMode_ = matrixMode;
  if(!matrixMode_)
    currentStateSize_ = stateSize_;
  else
    currentStateSize_ = stateSizeMatrix_;
}

} // namespace stateObservation
