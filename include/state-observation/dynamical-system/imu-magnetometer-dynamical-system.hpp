/**
 * \file      imu-magnetometer-dynamical-system.hpp
 * \author    Mehdi Benallegue, Joseph Mirabel
 * \date      2015
 * \brief     The file describes the dynamical system defined by an inertial
 *            measurement unit (IMU) fixed on a rigid body. It includes a
 *            magnetometer.
 *
 * \details
 *
 *
 */

#ifndef IMU_DYNAMICAL_SYSTEM_HPP
#define IMU_DYNAMICAL_SYSTEM_HPP

#include <state-observation/api.h>
#include <state-observation/dynamical-system/dynamical-system-functor-base.hpp>
#include <state-observation/noise/noise-base.hpp>
#include <state-observation/sensors-simulation/accelerometer-gyrometer-magnetometer.hpp>
#include <state-observation/tools/rigid-body-kinematics.hpp>

namespace stateObservation
{

/**
 * \class  IMUMagnetometerDynamicalSystem
 * \brief  The class is an implementation of the dynamical system defined by
 *         an inertial measurement unit (IMU) fixed on a rigid body. The state
 *         is the position velocity and acceleration and the orientation and rotation
 *         velocity and acceleration. The sensors are the accelerometer, the gyrometer
 *         and the magnetometer
 *
 *
 */
class STATE_OBSERVATION_DLLAPI IMUMagnetometerDynamicalSystem : public DynamicalSystemFunctorBase
{
public:
  /// The constructor
  IMUMagnetometerDynamicalSystem();

  /// The virtual destructor
  virtual ~IMUMagnetometerDynamicalSystem();

  /// Description of the state dynamics
  virtual Vector stateDynamics(const Vector & x, const Vector & u, TimeIndex k);

  /// Description of the sensor's dynamics
  virtual Vector measureDynamics(const Vector & x, const Vector & u, TimeIndex k);

  /// Sets a noise which disturbs the state dynamics
  virtual void setProcessNoise(NoiseBase *);
  /// Removes the process noise
  virtual void resetProcessNoise();
  /// Gets the process noise
  virtual NoiseBase * getProcessNoise() const;

  /// Sets a noise which disturbs the measurements
  virtual void setMeasurementNoise(NoiseBase *);
  /// Removes the measurement noise
  virtual void resetMeasurementNoise();
  /// Gets a pointer on the measurement noise
  virtual NoiseBase * getMeasurementNoise() const;

  /// Set the period of the time discretization
  virtual void setSamplingPeriod(double dt);

  /// Gets the state size
  virtual Index getStateSize() const;
  /// Gets the input size
  virtual Index getInputSize() const;
  /// Gets the measurement size
  virtual Index getMeasurementSize() const;

protected:
  typedef kine::indexes<kine::rotationVector> indexes;

  AccelerometerGyrometerMagnetometer sensor_;

  NoiseBase * processNoise_;

  double dt_;

  Vector3Unaligned orientationVector_;
  QuaternionUnaligned quaternion_;

  Quaternion computeQuaternion_(const Vector3 & x);

  static const Index stateSize_ = 18;
  static const Index inputSize_ = 0;
  static const Index measurementSize_ = 9;

private:
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
} // namespace stateObservation
#endif // IMU-DYNAMICAL-SYSTEM_HPP
