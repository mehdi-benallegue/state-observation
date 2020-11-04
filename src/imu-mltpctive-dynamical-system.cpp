#include <state-observation/dynamical-system/imu-mltpctive-dynamical-system.hpp>

#include <state-observation/tools/miscellaneous-algorithms.hpp>

namespace stateObservation
{

IMUMltpctiveDynamicalSystem::IMUMltpctiveDynamicalSystem()
: opt_(stateTangentSize_, measurementSize_), processNoise_(0x0), dt_(1)
{
#ifdef STATEOBSERVATION_VERBOUS_CONSTRUCTORS
  std::cout << std::endl << "IMUFixedContactDynamicalSystem Constructor" << std::endl;
#endif // STATEOBSERVATION_VERBOUS_CONSTRUCTOR
  // ctor
}

IMUMltpctiveDynamicalSystem::~IMUMltpctiveDynamicalSystem()
{
  // dtor
}

Vector IMUMltpctiveDynamicalSystem::stateDynamics(const Vector & x, const Vector & u, TimeIndex)
{
  assertStateVector_(x);
  assertInputVector_(u);

  Vector3 position = x.segment(indexes::pos, 3);
  Vector3 velocity = x.segment(indexes::linVel, 3);
  Vector3 acceleration = x.segment(indexes::linAcc, 3);

  Quaternion orientation(x.segment<4>(indexes::ori));
  Vector3 angularVelocity = x.segment(indexes::angVel, 3);
  Vector3 angularAcceleration = x.segment(indexes::angAcc, 3);

  kine::integrateKinematics(position, velocity, acceleration, orientation, angularVelocity, angularAcceleration, dt_);

  double norm2 = orientation.squaredNorm();

  if(norm2 > 1 + kine::quatNormTol || norm2 < 1 - kine::quatNormTol)
  {
    orientation.normalize();
  }

  // x_{k+1}
  Vector xk1 = Vector::Zero(indexes::size, 1);

  xk1.segment(indexes::pos, 3) = position;
  xk1.segment(indexes::linVel, 3) = velocity;

  xk1.segment<4>(indexes::ori) = orientation.coeffs();

  xk1.segment(indexes::angVel, 3) = angularVelocity;

  // inputs
  Vector3 accelerationInput = u.head(3);
  Vector3 angularAccelerationInput = u.tail(3);

  xk1.segment<3>(indexes::linAcc) = accelerationInput;
  xk1.segment<3>(indexes::angAcc) = angularAccelerationInput;

  if(processNoise_ != 0x0)
    return processNoise_->getNoisy(xk1);
  else
    return xk1;
}

Vector IMUMltpctiveDynamicalSystem::measureDynamics(const Vector & x, const Vector &, TimeIndex k)
{
  assertStateVector_(x);

  Vector3 acceleration = x.segment(indexes::linAcc, 3);

  Quaternion q(x.segment<4>(indexes::ori));

  Vector3 angularVelocity = x.segment(indexes::angVel, 3);

  Vector v = Vector::Zero(10, 1);

  v.head<4>() = q.coeffs();

  v.segment(4, 3) = acceleration;
  v.tail(3) = angularVelocity;

  sensor_.setState(v, k);

  return sensor_.getMeasurements();
}

Matrix IMUMltpctiveDynamicalSystem::getAMatrix(const Vector & xh)
{

  opt_.deltaR = xh.segment<3>(indexes::angVel) * dt_ + xh.segment<3>(indexes::angAcc) * dt_ * dt_ / 2;

  kine::derivateRotationMultiplicative(opt_.deltaR, opt_.jRR, opt_.jRv);

  opt_.AJacobian.block<3, 3>(indexesTangent::ori, indexesTangent::ori) = opt_.jRR;

  opt_.AJacobian.block<3, 3>(indexesTangent::pos, indexesTangent::linVel).diagonal().setConstant(dt_);
  opt_.AJacobian.block<3, 3>(indexesTangent::ori, indexesTangent::angVel) = opt_.jRv * dt_;
  opt_.AJacobian.block<6, 6>(indexesTangent::linVel, indexesTangent::linAcc).diagonal().setConstant(dt_);

  opt_.AJacobian.block<3, 3>(indexesTangent::pos, indexesTangent::linAcc).diagonal().setConstant(dt_ * dt_ * 0.5);
  opt_.AJacobian.block<3, 3>(indexesTangent::ori, indexesTangent::angAcc) = opt_.jRv * dt_ * dt_ / 2;

  return opt_.AJacobian;
}

Matrix IMUMltpctiveDynamicalSystem::getCMatrix(const Vector & xp)
{
  opt_.Rt = Quaternion(xp.segment<4>(indexes::ori)).toRotationMatrix().transpose();

  opt_.CJacobian.block<3, 3>(0, indexesTangent::ori).noalias() =
      opt_.Rt * kine::skewSymmetric(xp.segment<3>(indexes::linAcc) + cst::gravity);

  opt_.CJacobian.block<3, 3>(0, indexesTangent::linAcc).noalias() = opt_.Rt;

  opt_.CJacobian.block<3, 3>(3, indexesTangent::ori).noalias() =
      opt_.Rt * kine::skewSymmetric(xp.segment<3>(indexes::angVel));

  opt_.CJacobian.block<3, 3>(3, indexesTangent::angVel).noalias() = opt_.Rt;

  return opt_.CJacobian;
}

void IMUMltpctiveDynamicalSystem::setProcessNoise(NoiseBase * n)
{
  processNoise_ = n;
}

void IMUMltpctiveDynamicalSystem::resetProcessNoise()
{
  processNoise_ = 0x0;
}

void IMUMltpctiveDynamicalSystem::setMeasurementNoise(NoiseBase * n)
{
  sensor_.setNoise(n);
}
void IMUMltpctiveDynamicalSystem::resetMeasurementNoise()
{
  sensor_.resetNoise();
}

void IMUMltpctiveDynamicalSystem::setSamplingPeriod(double dt)
{
  dt_ = dt;
}

Index IMUMltpctiveDynamicalSystem::getStateSize() const
{
  return stateSize_;
}

Index IMUMltpctiveDynamicalSystem::getInputSize() const
{
  return inputSize_;
}

Index IMUMltpctiveDynamicalSystem::getMeasurementSize() const
{
  return measurementSize_;
}

NoiseBase * IMUMltpctiveDynamicalSystem::getProcessNoise() const
{
  return processNoise_;
}

NoiseBase * IMUMltpctiveDynamicalSystem::getMeasurementNoise() const
{
  return sensor_.getNoise();
}

void IMUMltpctiveDynamicalSystem::stateSum(const Vector & stateVector, const Vector & tangentVector, Vector & sum)
{
  sum.resize(indexes::size);
  sum.segment<3>(indexes::pos) = stateVector.segment<3>(indexes::pos) + tangentVector.segment<3>(indexesTangent::pos);
  sum.segment<3>(indexes::linVel) =
      stateVector.segment<3>(indexes::linVel) + tangentVector.segment<3>(indexesTangent::linVel);
  sum.segment<3>(indexes::linAcc) =
      stateVector.segment<3>(indexes::linAcc) + tangentVector.segment<3>(indexesTangent::linAcc);

  sum.segment<4>(indexes::ori) = (kine::rotationVectorToQuaternion(tangentVector.segment<3>(indexesTangent::ori))
                                  * Quaternion(stateVector.segment<4>(indexes::ori)))
                                     .coeffs();
  sum.segment<3>(indexes::angVel) =
      stateVector.segment<3>(indexes::angVel) + tangentVector.segment<3>(indexesTangent::angVel);
  sum.segment<3>(indexes::angAcc) =
      stateVector.segment<3>(indexes::angAcc) + tangentVector.segment<3>(indexesTangent::angAcc);
}

void IMUMltpctiveDynamicalSystem::stateDifference(const Vector & stateVector1,
                                                  const Vector & stateVector2,
                                                  Vector & difference)
{
  difference.resize(indexesTangent::size);
  difference.segment<3>(indexesTangent::pos) =
      stateVector1.segment<3>(indexes::pos) - stateVector2.segment<3>(indexes::pos);
  difference.segment<3>(indexesTangent::linVel) =
      stateVector1.segment<3>(indexes::linVel) - stateVector2.segment<3>(indexes::linVel);
  difference.segment<3>(indexesTangent::linAcc) =
      stateVector1.segment<3>(indexes::linAcc) - stateVector2.segment<3>(indexes::linAcc);

  difference.segment<3>(indexesTangent::ori) =
      kine::quaternionToRotationVector(Quaternion(stateVector1.segment<4>(indexes::ori))
                                       * Quaternion(stateVector2.segment<4>(indexes::ori)).conjugate());
  difference.segment<3>(indexesTangent::angVel) =
      stateVector1.segment<3>(indexes::angVel) - stateVector2.segment<3>(indexes::angVel);
  difference.segment<3>(indexesTangent::angAcc) =
      stateVector1.segment<3>(indexes::angAcc) - stateVector2.segment<3>(indexes::angAcc);
}

} // namespace stateObservation
