
///\file      lipm-dcm-estimator.hpp
///\author    Mehdi Benallegue
///\date      2020
///\brief     Filtering of divergent component of motion (DCM) and estimation of a bias betweeen the DCM
///           and the corresponding zero moment point for a linearized inverted
///           pendulum model
///
///\detail
///
///
#ifndef LIPMDCMBIASESTIMATOR_HPP
#define LIPMDCMBIASESTIMATOR_HPP

#include <state-observation/api.h>
#include <state-observation/observer/linear-kalman-filter.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>
#include <state-observation/tools/rigid-body-kinematics.hpp>

namespace stateObservation
{

/// \class LipmDcmBiasEstimator
/// \brief Filtering of divergent component of motion (DCM) and estimation of a bias betweeen the DCM
///           and the corresponding zero moment point for a linearized inverted
///        pendulum model.
///
/// \details
/// A humanoid robot can be modeled as an inverted pendulum. The dynamics can be
/// linearized to obtain a dynamics with a convergent and a divergent component of motion (DCN).
/// The dynamics of the DCM depends on the Zero Moment Point.
/// The DCM can be measured using the CoM and its velocity, but the CoM position can be biased.
/// This estimator uses Kalman Filtering to estimate this bias and give a delay-free filtering of the DCM.

class STATE_OBSERVATION_DLLAPI LipmDcmEstimator
{
private:
  constexpr static double defaultDt_ = 0.005;
  constexpr static double defaultOmega_ = tools::sqrt(cst::gravityConstant);

public:
  /// default expected drift of the bias every second
  constexpr static double defaultBiasDriftSecond = 0.002;

  /// default error in the estimation of the sensors
  constexpr static double defaultZmpErrorStd = 0.005;
  constexpr static double defaultDcmErrorStd = 0.01;

  /// default uncertainty in the initial values of DCM and Bias
  constexpr static double defaultDCMUncertainty = 0.01;
  constexpr static double defaultBiasUncertainty = 0.01;

  /// default value for Bias limit (0 means limitless)
  constexpr static double defaultBiasLimit = 0;

  /// @brief Construct a new Lipm Dcm Estimator object
  /// @details Use this if no DCM measurements are available or when a good guess of its unbiased position is
  /// available
  ///
  /// @param dt                     the sampling time in seconds
  /// @param omega_0                the natural frequency of the DCM (rad/s)
  /// @param biasDriftPerSecondStd  the standard deviation of the drift (m/s)
  /// @param initZMP                the initial value of the DCM
  /// @param initDCM                the initial value of the DCM
  /// @param initBias               the initial value of the bias
  /// @param dcmMeasureErrorStd     the standard deviation of the dcm estimation error, NOT including the bias (m)
  /// @param zmpMeasureErrorStd     the standard deviaiton of the zmp estimation error (m)
  /// @param biasLimit              the X and Y (expressed in local frame) largest accepted absolute values of the bias
  ///                               (zero means no limit)
  /// @param tinitDcmUncertainty    the uncertainty in the DCM initial value in meters
  /// @param initBiasUncertainty    the uncertainty in the bias initial value in meters
  LipmDcmEstimator(double dt = defaultDt_,
                   double omega_0 = defaultOmega_,
                   double biasDriftPerSecondStd = defaultBiasDriftSecond,
                   double dcmMeasureErrorStd = defaultDcmErrorStd,
                   double zmpMeasureErrorStd = defaultZmpErrorStd,
                   const Vector2 & biasLimit = Vector2::Constant(defaultBiasLimit),
                   const Vector2 & initZMP = Vector2::Zero(),
                   const Vector2 & initDcm = Vector2::Zero(),
                   const Vector2 & initBias = Vector2::Zero(),
                   const Vector2 & initDcmUncertainty = Vector2::Constant(defaultDCMUncertainty),
                   const Vector2 & initBiasUncertainty = Vector2::Constant(defaultBiasUncertainty));

  /// @brief Resets the estimator with first measurements
  /// @details Use this when initializing with an available DCM (biased / or not) measurement
  ///
  /// @param measuredDcm            the the measured position of the DCM in the world frame
  /// @param measuredZMP            the the measured position of the ZMP in the world frame
  /// @param yaw                    the initial yaw angle in the form of a rotation matrix
  /// @param measurementIsWithBias  sets if yes or no the first measurement is biased
  /// @param biasLimit              the X and Y (expressed in local frame) largest accepted absolute values of the bias
  ///                               (zero means no limit)
  /// @param initBias               the initial value of the drift
  /// @param initBiasuncertainty    the uncertainty in the bias initial value in meters
  void resetWithMeasurements(const Vector2 & measuredDcm,
                             const Vector2 & measuredZMP,
                             const Matrix2 & yaw = Matrix2::Identity(),
                             bool measurementIsWithBias = true,
                             const Vector2 & initBias = Vector2::Constant(0),
                             const Vector2 & initBiasuncertainty = Vector2::Constant(defaultBiasUncertainty));

  /// @brief Resets the estimator with first measurements
  /// @details Use this when initializing with an available DCM (biased / or not) measurement
  ///
  /// @param measuredDcm            the the measured position of the DCM in the world frame
  /// @param measuredZMP            the the measured position of the ZMP in the world frame
  /// @param yaw                    the initial yaw angle
  /// @param measurementIsWithBias  sets if yes or no the first measurement is biased
  /// @param biasLimit              the X and Y (expressed in local frame) largest accepted absolute values of the bias
  ///                               (zero means no limit)
  /// @param initBias               the initial value of the drift
  /// @param initBiasuncertainty    the uncertainty in the bias initial value in meters
  void resetWithMeasurements(const Vector2 & measuredDcm,
                             const Vector2 & measuredZMP,
                             double yaw = 0,
                             bool measurementIsWithBias = true,
                             const Vector2 & initBias = Vector2::Constant(0),
                             const Vector2 & initBiasuncertainty = Vector2::Constant(defaultBiasUncertainty))
  {
    resetWithMeasurements(measuredDcm, measuredZMP, Rotation2D(yaw).toRotationMatrix(), measurementIsWithBias, initBias,
                          initBiasuncertainty);
  }

  ///@brief Destroy the Lipm Dcm Bias Estimator object
  ~LipmDcmEstimator();

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
  void setSamplingTime(double dt);

  ///@brief Set the Bias object from a guess
  ///
  ///@param bias guess
  void setBias(const Vector2 & bias);

  ///@copydoc setBias(double bias)
  ///
  ///@param the uncertainty you have in this guess in meters
  void setBias(const Vector2 & bias, const Vector2 & uncertainty);

  /// @brief Set the Bias Drift Per Second
  ///
  /// @param driftPerSecond the standard deviation of the drift (m/s)
  void setBiasDriftPerSecond(double driftPerSecond);

  /// @brief Set the Bias Limit
  ///
  /// @param biasLimit the X and Y (expressed in local frame) largest accepted
  ///                   absolute values of the bias (zero means no limit)
  void setBiasLimit(const Vector2 & biasLimit);

  /// @brief set the unbiased DCM position from a guess
  ///
  /// @param dcm guess
  void setUnbiasedDCM(const Vector2 & dcm);

  /// @copydoc setUnbiasedDCM(double dcm)
  ///
  /// @param uncertainty the uncertainty in this guess
  void setUnbiasedDCM(const Vector2 & dcm, const Vector2 & uncertainty);

  /// @brief Set the Zmp Measurement Error Stamdard devbiation
  ///
  void setZmpMeasureErrorStd(double);

  /// @brief Set the Dcm Measurement Error Standard
  ///
  void setDcmMeasureErrorStd(double);

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
  void setInputs(const Vector2 & dcm, const Vector2 & zmp, const Matrix2 & R = Matrix2::Identity());

  /// @brief Runs the estimation. Needs to be called every timestep
  ///
  /// @return Vector2
  void update();

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

  Vector2 biasLimit_;

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
  inline static Matrix2 Vec2ToSqDiag(const Vector2 & v)
  {
    return Vector2(v.array().square()).asDiagonal();
  }

  /// @brief builds a constant 2x2 diagonal from a double
  inline static Matrix2 dblToDiag(const double & d)
  {
    return Vector2::Constant(d).asDiagonal();
  }

  /// @brief builds a constant 2x2 diagonal from a square of a double
  inline static Matrix2 dblToSqDiag(const double & d)
  {
    return dblToDiag(d * d);
  }

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // namespace stateObservation

#endif /// LIPMDCMBIASESTIMATOR_HPP
