/**
 * \file      model-base-ekf-flex-estimator-imu.hpp
 * \author    Mehdi Benallegue
 * \date      2013
 * \brief     Declares the class of the estimation of the flexibility using an
 *            extended Kalman filter and a fixed contact hypothesis
 *
 * \details
 *
 *
 */

#ifndef FLEXBILITYESTMATOR_MODELBASEEKFFLEXIBILITYESTIMATOR_IMU_H
#define FLEXBILITYESTMATOR_MODELBASEEKFFLEXIBILITYESTIMATOR_IMU_H

#include <state-observation/api.h>
#include <state-observation/flexibility-estimation/ekf-flexibility-estimator-base.hpp>
//#include <state-observation/flexibility-estimation/stable-imu-fixed-contact-dynamical-system.hpp>
#include <state-observation/flexibility-estimation/imu-elastic-local-frame-dynamical-system.hpp>
//#include <state-observation/flexibility-estimation/imu-fixed-contact-dynamical-system.hpp>

namespace stateObservation
{
namespace flexibilityEstimation
{

/**
 * \class  ModelBaseEKFFlexEstimatorIMU
 * \brief  This class implements the flexibility estimation of a robot with
 *         the hypothesis that the contact positions do not move. This constraint
 *         is expressed using fictious measurements but the interface is transparent
 *         to this assumption, the state is expressed using classical representation
 *         of position, velocity, acceleration, orientation (using (theta x mu) representation)
 *         angular velocity (omega) and acceleration (omega dot)
 *
 */

class STATE_OBSERVATION_DLLAPI ModelBaseEKFFlexEstimatorIMU : public EKFFlexibilityEstimatorBase,
                                                              private boost::noncopyable
{
public:
  struct contactModel
  {
    /// indexes of the different components of a vector of the input state
    static const unsigned elasticContact = IMUElasticLocalFrameDynamicalSystem::contactModel::elasticContact;
    static const unsigned pendulum = IMUElasticLocalFrameDynamicalSystem::contactModel::pendulum;
  };

  /// The constructor, it requires the value of the time discretization period
  explicit ModelBaseEKFFlexEstimatorIMU(double dt = 0.005);

  /// Virtual destructor
  virtual ~ModelBaseEKFFlexEstimatorIMU();

  /// Sets the number of contacts can be changed online
  void setContactsNumber(unsigned i);

  unsigned getContactsNumber()
  {
    return functor_.getContactsNumber();
  }

  IMUElasticLocalFrameDynamicalSystem getFunctor()
  {
    return functor_;
  }

  virtual stateObservation::Vector computeAccelerations()
  {
    return functor_.computeAccelerations(getFlexibilityVector(), getInput());
  }

  void setContactModel(unsigned nb);

  /// Sets the value of the next sensor measurement y_{k+1}
  virtual void setMeasurement(const Vector & y);

  /// Sets the process covariance matrice
  virtual void setProcessNoiseCovariance(const Matrix & Q);

  /// Sets the measurements covariance matrice
  virtual void setMeasurementNoiseCovariance(const Matrix & R);

  /// gets the covariance matrices for the process noises
  virtual Matrix getProcessNoiseCovariance() const;

  /// gets the covariance matrices for the sensor noises
  virtual Matrix getMeasurementNoiseCovariance() const;

  virtual Vector getMomentaDotFromForces();
  virtual Vector getMomentaDotFromKinematics();
  virtual Vector getForcesAndMoments();

  // get state covariance
  stateObservation::Vector getStateCovariance() const
  {
    stateObservation::Matrix P(ekf_.getStateCovariance());
    stateObservation::Vector Pvec(ekf_.getStateSize());
    for(Index i = 0; i < ekf_.getStateSize(); ++i) Pvec(i) = P(i, i);
    return Pvec;
  }

  virtual void setComBiasGuess(const stateObservation::Vector & x);

  /// Sets a value of the flexibility x_k provided from another source
  /// can be used for initialization of the estimator
  virtual void setFlexibilityGuess(const Matrix & x);

  /// Gets an estimation of the flexibility in the form of a homogeneous matrix
  virtual Matrix4 getFlexibility();

  /// Gets an estimation of the flexibility in the form of a state vector \hat{x_{k+1}}
  virtual const Vector & getFlexibilityVector();

  virtual stateObservation::Matrix & computeLocalObservationMatrix();
  virtual stateObservation::Matrix getAMatrix()
  {
    return ekf_.getA();
  }

  virtual stateObservation::Matrix getCMatrix()
  {
    return ekf_.getC();
  }

  virtual Index getMeasurementSize() const;

  virtual Index getStateSize() const;

  virtual Index getInputSize() const;

  /// sets to whether or not the force mesurements are taken into account
  virtual void setWithForcesMeasurements(bool);

  bool getWithForcesMeasurements();

  virtual void setWithAbsolutePos(bool);

  void setWithUnmodeledForces(bool b);

  bool getWithUnmodeledForces()
  {
    return withUnmodeledForces_;
  }

  bool getWithAbsolutePos()
  {
    return withAbsolutePos_;
  }

  virtual void setWithComBias(bool b);

  virtual bool getWithComBias()
  {
    return withComBias_;
  }

  virtual void setUnmodeledForceVariance(double d);
  virtual void setUnmodeledForceProcessVariance(double d);

  /// sets the force sensor variance, either with a diagonal matrix with
  /// constant diagonal, or using a 3x3 matrix
  virtual void setForceVariance(double d);
  virtual void setForceVariance(const Matrix3 & C);

  virtual void setAbsolutePosVariance(double d);

  /// sets the sampling period
  virtual void setSamplingPeriod(double);

  /// Enable or disable the estimation
  void setOn(bool & b);

  virtual void setKfe(const Matrix3 & m);
  virtual void setKfv(const Matrix3 & m);
  virtual void setKte(const Matrix3 & m);
  virtual void setKtv(const Matrix3 & m);

  virtual void setKfeRopes(const Matrix3 & m);
  virtual void setKfvRopes(const Matrix3 & m);
  virtual void setKteRopes(const Matrix3 & m);
  virtual void setKtvRopes(const Matrix3 & m);

  virtual void setPe(const stateObservation::Vector3 & Pe)
  {
    functor_.setPe(Pe);
  }

  virtual Matrix getKfe() const;
  virtual Matrix getKfv() const;
  virtual Matrix getKte() const;
  virtual Matrix getKtv() const;

  /// Resets the covariance matrices to their original values
  virtual void resetCovarianceMatrices();
  virtual void resetStateCovarianceMatrix();

  virtual void setRobotMass(double m);
  virtual double getRobotMass() const
  {
    return functor_.getRobotMass();
  }

  void setTorquesLimit(const Vector3 & v)
  {
    limitTorques_ = v;
  }

  void setForcesLimit(const Vector3 & v)
  {
    limitForces_ = v;
  }

  virtual stateObservation::Vector3 getForcesLimit() const
  {
    return limitForces_;
  }

  virtual stateObservation::Vector3 getTorquesLimit() const
  {
    return limitTorques_;
  }

  void setLimitOn(const bool & b)
  {
    limitOn_ = b;
  }

  virtual bool getLimitOn() const
  {
    return limitOn_;
  }

  static Matrix getDefaultQ();

  static Matrix6 getDefaultRIMU();

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
  virtual void updateMeasurementCovarianceMatrix_();

  IMUElasticLocalFrameDynamicalSystem functor_;

  Vector x_;

  Matrix R_, Q_, P_;

  const Index stateSize_;

  static const Index measurementSizeBase_ = 12;

  static const Index inputSizeBase_ = IMUElasticLocalFrameDynamicalSystem::input::sizeBase;
  Index inputSize_;

  double dt_; // sampling period
  bool on_;

  double unmodeledForceVariance_;
  Matrix forceVariance_; // force sensor variance
  double absPosVariance_;

  bool useFTSensors_;

  bool withAbsolutePos_;
  bool withComBias_;
  bool withUnmodeledForces_;

  Vector3 limitTorques_;
  Vector3 limitForces_;
  bool limitOn_;

  struct optimization
  {
    stateObservation::Matrix O;
    stateObservation::Matrix CA;
  } op_;

private:
};

} // namespace flexibilityEstimation
} // namespace stateObservation
#endif // FLEXBILITYESTMATOR_MODELBASEEKFFLEXIBILITYESTIMATOR_IMU_H
