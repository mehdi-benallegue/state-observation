
///\file      lipm-dcm-bias-estimator.hpp
///\author    Mehdi Benallegue
///\date      2020
///\brief     Estimation of a bias betweeen the divergent component of motion
///           and the corresponding zero moment point for a linearized inverted
///           pendulum model
///
///\detail
///
///
#ifndef LIPMDCMBIASESTIMATOR_HPP
#define LIPMDCMBIASESTIMATOR_HPP

#include <state-observation/api.h>
#include <state-observation/dynamics-estimators/unidim-lipm-dcm-bias-estimator.hpp>
#include <state-observation/tools/rigid-body-kinematics.hpp>

namespace stateObservation
{

/// \class LipmDcmBiasEstimator
/// \brief Estimation of a bias betweeen the divergent component of motion
///        and the corresponding zero moment point for a linearized inverted
///        pendulum model.
///
/// \details
/// A humanoid robot can be modeled as an inverted pendulum. The dynamics can be
/// linearized to obtain a dynamics with a convergent and a divergent component of motion (DCN).
/// The dynamics of the DCM depends on the Zero Moment Point.
/// The DCM can be measured using the CoM and its velocity, but the CoM position can be biased.
/// This estimator uses Kalman Filtering to estimate this bias in one axis.

class STATE_OBSERVATION_DLLAPI LipmDcmBiasEstimator
{
private:
  constexpr static double defaultDt_ = 0.005;

  /// default expected drift of the bias every second
  constexpr static double defaultBiasDriftSecond_ = 0.002;

  /// default error in the estimation of the sensors
  constexpr static double defaultZmpErrorStd_ = 0.005;
  constexpr static double defaultDcmErrorStd_ = 0.01;

  /// default uncertainty in the initial values of DCM and Bias
  constexpr static double defaultDCMUncertainty = 0.01;
  constexpr static double defaultBiasUncertainty = 0.01;

public:
  /// @brief Construct a new Lipm Dcm Bias Estimator object
  /// @details Use this if no DCM measurements are available or when a good guess of its unbiased position is available
  ///
  /// @param omega_0                the natural frequency of the DCM (rad/s)
  /// @param dt                     the sampling time in seconds
  /// @param biasDriftPerSecondStd  the standard deviation of the drift (m/s)
  /// @param zmpMeasureErrorStd     the standard deviaiton of the zmp estimation error (m)
  /// @param dcmMeasureErrorStd     the standard deviation of the dcm estimation error, NOT including the bias (m)
  /// @param initDCM                the initial value of the DCM
  /// @param initBias               the initial value of the bias
  /// @param tinitDcmUncertainty    the uncertainty in the DCM initial value in meters
  /// @param initBiasUncertainty    the uncertainty in the bias initial value in meters
  LipmDcmBiasEstimator(double omega_0 = 1,
                       double dt = defaultDt_,
                       double biasDriftPerSecondStd = defaultBiasDriftSecond_,
                       const Vector2 & initZMP = Vector2::Zero(),
                       const Vector2 & initDcm = Vector2::Zero(),
                       const Vector2 & initBias = Vector2::Zero(),
                       double dcmMeasureErrorStd = defaultDcmErrorStd_,
                       double zmpMeasureErrorStd = defaultZmpErrorStd_,
                       const Vector2 & initDcmUncertainty = Vector2::Constant(defaultDCMUncertainty),
                       const Vector2 & initBiasUncertainty = Vector2::Constant(defaultBiasUncertainty));

  /// @brief Resets the estimator with first measurements
  /// @details Use this when initializing with an available DCM (biased / or not) measurement
  ///
  /// @param measuredDcm            the the measured position of the DCM in the world frame
  /// @param measuredZMP            the the measured position of the ZMP in the world frame
  /// @param yaw                    the initial yaw angle in the form of a rotation matrix
  /// @param measurementIsWithBias  sets if yes or no the first measurement is biased
  /// @param biasDriftPerSecondStd  the standard deviation of the drift (m/s)
  /// @param zmpMeasureErrorStd     the standard deviaiton of the zmp estimation error (m)
  /// @param dcmMeasureErrorStd     the standard deviation of the dcm estimation error, NOT including the bias (m)
  /// @param initBias               the initial value of the drift
  /// @param initBiasuncertainty    the uncertainty in the bias initial value in meters
  void resetWithMeasurements(const Vector2 & measuredDcm,
                             const Vector2 & measuredZMP,
                             const Matrix2 & yaw = Matrix2::Identity(),
                             bool measurementIsWithBias = true,
                             double biasDriftPerSecondStd = defaultBiasDriftSecond_,
                             double dcmMeasureErrorStd = defaultDcmErrorStd_,
                             double zmpMeasureErrorStd = defaultZmpErrorStd_,
                             const Vector2 & initBias = Vector2::Constant(0),
                             const Vector2 & initBiasuncertainty = Vector2::Constant(defaultBiasUncertainty));

  /// @brief Resets the estimator with first measurements
  /// @details Use this when initializing with an available DCM (biased / or not) measurement
  ///
  /// @param measuredDcm            the the measured position of the DCM in the world frame
  /// @param measuredZMP            the the measured position of the ZMP in the world frame
  /// @param yaw                    the initial yaw angle
  /// @param measurementIsWithBias  sets if yes or no the first measurement is biased
  /// @param biasDriftPerSecondStd  the standard deviation of the drift (m/s)
  /// @param zmpMeasureErrorStd     the standard deviaiton of the zmp estimation error (m)
  /// @param dcmMeasureErrorStd     the standard deviation of the dcm estimation error, NOT including the bias (m)
  /// @param initBias               the initial value of the drift
  /// @param initBiasuncertainty    the uncertainty in the bias initial value in meters
  inline void resetWithMeasurements(const Vector2 & measuredDcm,
                                    const Vector2 & measuredZMP,
                                    double yaw,
                                    bool measurementIsWithBias = true,
                                    double biasDriftPerSecondStd = defaultBiasDriftSecond_,
                                    double dcmMeasureErrorStd = defaultDcmErrorStd_,
                                    double zmpMeasureErrorStd = defaultZmpErrorStd_,
                                    const Vector2 & initBias = Vector2::Constant(0),
                                    const Vector2 & initBiasuncertainty = Vector2::Constant(defaultBiasUncertainty))
  {
    resetWithMeasurements(measuredDcm, measuredZMP, Rotation2D(yaw).toRotationMatrix(), measurementIsWithBias,
                          biasDriftPerSecondStd, dcmMeasureErrorStd, zmpMeasureErrorStd, initBias, initBiasuncertainty);
  }

  /// @brief Resets the estimator with first measurements
  /// @details Use this when initializing with an available DCM (biased / or not) measurement
  ///
  /// @param measuredDcm            the the measured position of the DCM in the world frame
  /// @param measuredZMP            the the measured position of the ZMP in the world frame
  /// @param rotation                the 3d orientation from which the initial yaw angle will be extracted using the
  /// angle agnostic approach. This orientation is from local to global. i.e. bias_global == orientation * bias*local
  /// @param measurementIsWithBias  sets if yes or no the first measurement is biased
  /// @param biasDriftPerSecondStd  the standard deviation of the drift (m/s)
  /// @param zmpMeasureErrorStd     the standard deviaiton of the zmp estimation error (m)
  /// @param dcmMeasureErrorStd     the standard deviation of the dcm estimation error, NOT including the bias (m)
  /// @param initBias               the initial value of the drift
  /// @param initBiasuncertainty    the uncertainty in the bias initial value in meters
  inline void resetWithMeasurements(const Vector2 & measuredDcm,
                                    const Vector2 & measuredZMP,
                                    const Matrix3 & rotation,
                                    bool measurementIsWithBias = true,
                                    double biasDriftPerSecondStd = defaultBiasDriftSecond_,
                                    double zmpMeasureErrorStd = defaultZmpErrorStd_,
                                    double dcmMeasureErrorStd = defaultDcmErrorStd_,
                                    const Vector2 & initBias = Vector2::Constant(0),
                                    const Vector2 & initBiasuncertainty = Vector2::Constant(defaultBiasUncertainty))
  {
    resetWithMeasurements(measuredDcm, measuredZMP, kine::rotationMatrixToYawAxisAgnostic(rotation),
                          measurementIsWithBias, biasDriftPerSecondStd, zmpMeasureErrorStd, dcmMeasureErrorStd,
                          initBias, initBiasuncertainty);
  }

  ///@brief Destroy the Lipm Dcm Bias Estimator object
  ~LipmDcmBiasEstimator();

  ///@brief Set the Lipm Natural Frequency
  ///
  ///@param omega_0  is the sampling time in seconds
  void setLipmNaturalFrequency(double omega_0);

  /// @brief Get the Lipm Natural Frequency
  ///
  /// @return double
  inline double getLipmNaturalFrequency() const
  {
    return omega0_;
  }

  ///@brief Set the Sampling Time
  ///
  ///@param dt sampling time
  void setBias(const Vector2 & bias);

  ///@copydoc setBias(double bias)
  ///
  ///@param the uncertainty you have in this guess in meters
  void setBias(const Vector2 & bias, const Vector2 & uncertainty);

  /// @brief Set the Bias Drift Per Second
  ///
  /// @param driftPerSecond the standard deviation of the drift (m/s)
  void setBiasDriftPerSecond(double driftPerSecond);

  /// @brief set the real DCM position from a guess
  ///
  /// @param dcm guess
  void setDCM(const Vector2 & dcm);

  /// @copydoc setDCM(double dcm)
  ///
  /// @param dcm
  /// @param uncertainty the uncertainty in this guess
  void setDCM(const Vector2 & dcm, const Vector2 & uncertainty);

  /// @brief Set the Zmp Measurement Error Stamdard devbiation
  ///
  void setZmpMeasureErrorStd(double);

  /// @brief Set the Dcm Measurement Error Standard
  ///
  void setDcmMeasureErrorStd(double);

  /// @brief Set the Inputs of the estimator
  /// @details this version is with no orientation. The orientation will be assumed to be constant foir
  /// this sample
  ///
  /// @param dcm measurement of the DCM in the world frame
  /// @param zmp mesaurement of the ZMP in the world frame
  void setInputs(const Vector2 & dcm, const Vector2 & zmp);

  /// @brief Set the Inputs of the estimator.
  /// @details The yaw will be extracted from the orientation using the axis agnostic
  /// approach.
  ///
  /// @param dcm         measurement of the DCM in the world frame
  /// @param zmp         mesaurement of the ZMP in the world frame
  /// @param orientation the 3d orientation from which the yaw will be extracted. This orientation is from local to
  /// global. i.e. bias_global == orientation * bias*local
  ///
  inline void setInputs(const Vector2 & dcm, const Vector2 & zmp, const Matrix3 & orientation)
  {
    setInputs(dcm, zmp, kine::rotationMatrixToYawAxisAgnostic(orientation));
  }

  /// @brief Set the Inputs of the estimator.
  ///
  /// @param dcm measurement of the DCM in the world frame
  /// @param zmp mesaurement of the ZMP in the world frame
  /// @param yaw is the yaw angle to be used. This orientation is from local to global. i.e. bias_global == R *
  /// bias*local
  inline void setInputs(const Vector2 & dcm, const Vector2 & zmp, double yaw)
  {
    setInputs(dcm, zmp, Rotation2D(yaw).toRotationMatrix());
  }

  /// @brief Set the Inputs of the estimator.
  ///
  /// @param dcm  measurement of the DCM in the world frame
  /// @param zmp  mesaurement of the ZMP in the world frame
  /// @param R    the 2x2 Matrix'representing the yaw angle i.e. bias_global == R * bias*local
  void setInputs(const Vector2 & dcm, const Vector2 & zmp, const Matrix2 & R);

  /// @brief Runs the estimation. Needs to be called every timestep
  ///
  /// @return Vector2
  inline Vector4 update()
  {
    A_.bottomRightCorner<2, 2>().setIdentity(); /// reset the rotation part
    return filter_.getEstimatedState(filter_.getMeasurementTime());
  }

  /// @brief Get the Unbiased DCM filtered by the estimator
  ///
  /// @detailt This is the recommended output to take
  /// @return double
  Vector2 getUnbiasedDCM() const;

  /// @brief Get the estimated Bias
  ///
  /// @return double
  Vector2 getBias() const;

  /// @brief Get the estimated Bias expressed in the local frame of the robot
  ///
  /// @return double
  inline Vector2 getLocalBias() const
  {
    return previousOrientation_.transpose() * getBias();
  }

  /// @brief Get the Kalman Filter object
  /// This can be used to run specific Advanced Kalman filter related funcions
  /// @return LinearKalmanFilter&
  inline LinearKalmanFilter & getFilter()
  {
    return filter_;
  }

  /// @copydoc getFilter()
  /// const version
  inline const LinearKalmanFilter & getFilter() const
  {
    return filter_;
  }

protected:
  typedef Eigen::Matrix<double, 4, 2> Matrix42;
  typedef Eigen::Matrix<double, 2, 4> Matrix24;

  /// @brief set Matrices: A, B, Q
  void updateMatricesABQ_();

  double omega0_;
  double dt_;
  double biasDriftStd_;
  double zmpErrorStd_;

  Vector2 previousZmp_;

  LinearKalmanFilter filter_;
  Matrix4 A_;
  Matrix42 B_;
  /// this needs to be transposed
  Matrix24 C_;
  /// measurement noise
  Matrix2 R_;

  /// process noise
  Matrix4 Q_;

  Matrix2 previousOrientation_;

  /// @brief builds a diagonal out of the square valued of the Vec2
  ///
  inline static Matrix2 Vec2ToSqDiag(const Vector2 & v)
  {
    return Vector2(v.array().square()).asDiagonal();
  }

  /// @brief builds a constant 2x2 diagonal from a double
  ///
  inline static Matrix2 dblToDiag(const double & d)
  {
    return Vector2::Constant(d).asDiagonal();
  }

  /// @brief builds a constant 2x2 diagonal from a square of a double
  ///
  inline static Matrix2 dblToSqDiag(const double & d)
  {
    return dblToDiag(d * d);
  }

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // namespace stateObservation

#endif /// LIPMDCMBIASESTIMATOR_HPP
