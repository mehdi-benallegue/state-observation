#include <state-observation/dynamical-system/imu-magnetometer-dynamical-system.hpp>

#include <state-observation/tools/miscellaneous-algorithms.hpp>

namespace stateObservation
{

IMUMagnetometerDynamicalSystem::IMUMagnetometerDynamicalSystem()
: processNoise_(0x0), dt_(1), orientationVector_(Vector3::Zero()), quaternion_(Quaternion::Identity())
{
#ifdef STATEOBSERVATION_VERBOUS_CONSTRUCTORS
  std::cout << std::endl << "IMUFixedContactDynamicalSystem Constructor" << std::endl;
#endif // STATEOBSERVATION_VERBOUS_CONSTRUCTOR
       // ctor
}

IMUMagnetometerDynamicalSystem::~IMUMagnetometerDynamicalSystem()
{
  // dtor
}

Vector IMUMagnetometerDynamicalSystem::stateDynamics(const Vector & x, const Vector &, TimeIndex)
{
  assertStateVector_(x);

  Vector3 position = x.segment(indexes::pos, 3);
  Vector3 velocity = x.segment(indexes::linVel, 3);
  Vector3 acceleration = x.segment(indexes::linAcc, 3);

  Vector3 orientationV = x.segment(indexes::ori, 3);
  Vector3 angularVelocity = x.segment(indexes::angVel, 3);
  Vector3 angularAcceleration = x.segment(indexes::angAcc, 3);

  Quaternion orientation = computeQuaternion_(orientationV);

  kine::integrateKinematics(position, velocity, acceleration, orientation, angularVelocity, angularAcceleration, dt_);

  // x_{k+1}
  Vector xk1 = Vector::Zero(18, 1);

  xk1.segment(indexes::pos, 3) = position;
  xk1.segment(indexes::linVel, 3) = velocity;

  AngleAxis orientationAA(orientation);

  orientationV = orientationAA.angle() * orientationAA.axis();

  xk1.segment(indexes::ori, 3) = orientationV;
  xk1.segment(indexes::angVel, 3) = angularVelocity;

  xk1.segment(indexes::linAcc, 3).setZero();
  xk1.segment(indexes::angAcc, 3).setZero();

  if(processNoise_ != 0x0)
    return processNoise_->getNoisy(xk1);
  else
    return xk1;
}

Quaternion IMUMagnetometerDynamicalSystem::computeQuaternion_(const Vector3 & x)
{
  if(orientationVector_ != x)
  {
    quaternion_ = kine::rotationVectorToAngleAxis(x);
    orientationVector_ = x;
  }

  return quaternion_;
}

Vector IMUMagnetometerDynamicalSystem::measureDynamics(const Vector & x, const Vector &, TimeIndex k)
{
  assertStateVector_(x);

  Vector3 acceleration = x.segment(indexes::linAcc, 3);

  Vector3 orientationV = x.segment(indexes::ori, 3);
  Vector3 angularVelocity = x.segment(indexes::angVel, 3);

  Quaternion q = computeQuaternion_(orientationV);

  Vector v = Vector::Zero(10, 1);

  v.head<4>() = q.coeffs();

  v.segment(4, 3) = acceleration;
  v.tail(3) = angularVelocity;

  sensor_.setState(v, k);

  return sensor_.getMeasurements();
}

void IMUMagnetometerDynamicalSystem::setProcessNoise(NoiseBase * n)
{
  processNoise_ = n;
}

void IMUMagnetometerDynamicalSystem::resetProcessNoise()
{
  processNoise_ = 0x0;
}

void IMUMagnetometerDynamicalSystem::setMeasurementNoise(NoiseBase * n)
{
  sensor_.setNoise(n);
}
void IMUMagnetometerDynamicalSystem::resetMeasurementNoise()
{
  sensor_.resetNoise();
}

void IMUMagnetometerDynamicalSystem::setSamplingPeriod(double dt)
{
  dt_ = dt;
}

Index IMUMagnetometerDynamicalSystem::getStateSize() const
{
  return stateSize_;
}

Index IMUMagnetometerDynamicalSystem::getInputSize() const
{
  return inputSize_;
}

Index IMUMagnetometerDynamicalSystem::getMeasurementSize() const
{
  return measurementSize_;
}

NoiseBase * IMUMagnetometerDynamicalSystem::getProcessNoise() const
{
  return processNoise_;
}

NoiseBase * IMUMagnetometerDynamicalSystem::getMeasurementNoise() const
{
  return sensor_.getNoise();
}
} // namespace stateObservation
