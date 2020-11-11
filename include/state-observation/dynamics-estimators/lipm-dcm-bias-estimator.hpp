
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
#include <state-observation/observer/linear-kalman-filter.hpp>

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
  /// @param omega_0             the natural frequency of the DCM (rad/s)
  /// @param dt                  the sampling time in seconds
  /// @param biasDriftStd        the standard deviation of the drift (m/s)
  /// @param zmpMeasureErrorStd  the standard deviaiton of the zmp estimation error (m)
  /// @param dcmMeasureErrorStd  the standard deviation of the dcm estimation error, NOT including the bias (m)
  /// @param initDCM             the initial value of the DCM
  /// @param initBias            the initial value of the bias
  /// @param tinitDcmUncertainty the uncertainty in the DCM initial value in meters
  /// @param initBiasUncertainty the uncertainty in the bias initial value in meters
  LipmDcmBiasEstimator(double omega_0,
                       double dt = defaultDt_,
                       double biasDriftStd = defaultBiasDriftSecond_,
                       double zmpMeasureErrorStd = defaultZmpErrorStd_,
                       double dcmMeasureErrorStd = defaultDcmErrorStd_,
                       double initDcm = 0,
                       double initBias = 0,
                       double initDcmUncertainty = defaultDCMUncertainty,
                       double initBiasUncertainty = defaultBiasUncertainty);

  /// @brief Construct a new Lipm Dcm Bias Estimator object
  /// @details Use this when initializing with an available DCM (biased) measurement
  ///
  /// @param measurementIsWithBias  sets if yes or no the measurement is with the bias. It is better set to true than
  /// trying to remove it beforehand
  /// @param measuredDcm            the the measured position of the DCM
  /// @param omega_0                the natural frequency of the DCM (rad/s)
  /// @param dt                     the sampling time in seconds
  /// @param biasDriftStd           the standard deviation of the drift (m/s)
  /// @param zmpMeasureErrorStd     the standard deviaiton of the zmp estimation error (m)
  /// @param scmMeasureErrorStd     the standard deviation of the dcm estimation error, NOT including the bias (m)
  /// @param initBias               the initial value of the drift
  /// @param initBiasUncertainty the uncertainty in the bias initial value in meters
  LipmDcmBiasEstimator(bool measurementIsWithBias,
                       double measuredDcm,
                       double omega_0,
                       double dt = defaultDt_,
                       double biasDriftStd = defaultBiasDriftSecond_,
                       double zmpMeasureErrorStd = defaultZmpErrorStd_,
                       double dcmMeasureErrorStd = defaultDcmErrorStd_,
                       double initBias = 0,
                       double initBiasuncertainty = defaultBiasUncertainty);

  ///@brief Destroy the Lipm Dcm Bias Estimator object
  ~LipmDcmBiasEstimator();

  ///@brief Set the Lipm Natural Frequency
  ///
  ///@param omega_0  is the sampling time in seconds
  void setLipmNaturalFrequency(double omega_0);

  ///@brief Set the Sampling Time
  ///
  ///@param dt sampling time
  void setSamplingTime(double dt);

  ///@brief Set the Bias object from a guess
  ///
  ///@param bias guess
  void setBias(double bias);

  ///@copydoc setBias(double bias)
  ///
  ///@param the uncertainty you have in this guess in meters
  void setBias(double bias, double uncertainty);

  /// @brief Set the Bias Drift Per Second
  ///
  /// @param driftPerSecond the standard deviation of the drift (m/s)
  void setBiasDriftPerSecond(double driftPerSecond);

  /// @brief set the real DCM position from a guess
  ///
  /// @param dcm guess
  void setDCM(double dcm);

  /// @copydoc setDCM(double dcm)
  ///
  /// @param dcm
  /// @param uncertainty the uncertainty in this guess
  void setDCM(double dcm, double uncertainty);

  /// @brief Set the Zmp Measurement Error Stamdard devbiation
  ///
  void setZmpMeasureErrorStd(double);

  /// @brief Set the Dcm Measurement Error Standard
  ///
  void setDcmMeasureErrorStd(double);

  /// @brief Set the Inputs of the estimator
  ///
  /// @param dcm
  /// @param zmp
  void setInputs(double dcm, double zmp);

  /// @brief Runs the estimation. Needs to be called every timestep
  ///
  /// @return Vector2
  inline Vector2 update()
  {
    return filter_.getEstimatedState(filter_.getMeasurementTime());
  }

  /// @brief Get the Unbiased DCM filtered by the estimator
  ///
  /// @detailt This is the recommended output to take
  /// @return double
  double getUnbiasedDCM() const;

  /// @brief Get the estimated Bias
  ///
  /// @return double
  double getBias() const;

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

private:
  /// @brief set Matrices: A, B, Q
  void updateMatricesABQ_();

  double omega0_;
  double dt_;
  double biasDriftStd_;
  double zmpErrorStd_;

  LinearKalmanFilter filter_;
  Matrix2 A_;
  Vector2 B_;
  /// this needs to be transposed
  Vector2 C_;
  /// measurement noise
  Matrix1 R_;

  /// process noise
  Matrix2 Q_;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // namespace stateObservation

#endif /// LIPMDCMBIASESTIMATOR_HPP
