#include <state-observation/dynamics-estimators/lipm-dcm-bias-estimator.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>

namespace stateObservation
{
using namespace tools;
LipmDcmBiasEstimator::LipmDcmBiasEstimator(double omega_0,
                                           double dt,
                                           double biasDriftStd,
                                           double initZMP,
                                           double initDcm,
                                           double initBias,
                                           double zmpMeasureErrorStd,
                                           double dcmMeasureErrorStd,
                                           double initDcmUncertainty,
                                           double initBiasUncertainty)
: omega0_(omega_0), dt_(dt), biasDriftStd_(biasDriftStd), zmpErrorStd_(zmpMeasureErrorStd), previousZmp_(initZMP),
  filter_(2, 1, 1)
{
  updateMatricesABQ_();
  C_ << 1., 1.;
  R_ << square(dcmMeasureErrorStd);
  filter_.setC(C_);
  filter_.setMeasurementCovariance(R_);
  Vector2 x;
  x << initDcm, initBias;
  filter_.setState(x, 0);
  Matrix2 P;
  // clang-format off
  P<< square(initDcmUncertainty),  0,
      0,                           square(initBiasUncertainty);
  // clang-format on
  filter_.setStateCovariance(P);
}

LipmDcmBiasEstimator::LipmDcmBiasEstimator(bool measurementIsWithBias,
                                           double measuredDcm,
                                           double measuredZMP,
                                           double omega_0,
                                           double dt,
                                           double biasDriftStd,
                                           double zmpMeasureErrorStd,
                                           double dcmMeasureErrorStd,
                                           double initBias,
                                           double initBiasuncertainty)
: omega0_(omega_0), dt_(dt), biasDriftStd_(biasDriftStd), zmpErrorStd_(zmpMeasureErrorStd), previousZmp_(measuredZMP),
  filter_(2, 1, 1)
{
  updateMatricesABQ_();
  C_ << 1., 1.;
  R_ << square(dcmMeasureErrorStd);
  filter_.setC(C_.transpose());
  filter_.setMeasurementCovariance(R_);
  Vector2 x;

  /// initialize the state using the measurement
  if(measurementIsWithBias)
  {
    x << measuredDcm - initBias, initBias;
  }
  else
  {
    x << measuredDcm, initBias;
  }

  filter_.setState(x, 0);
  Matrix2 P;

  if(measurementIsWithBias)
  {
    /// The state and the
    // clang-format off
    P<< square(dcmMeasureErrorStd) + square(initBiasuncertainty),  -square(initBiasuncertainty),
        -square(initBiasuncertainty),                               square(initBiasuncertainty);
    // clang-format on
  }
  else
  {
    // clang-format off
    P<< square(dcmMeasureErrorStd),  0,
        0,                           square(initBiasuncertainty);
    // clang-format on
  }

  filter_.setStateCovariance(P);
}

LipmDcmBiasEstimator::~LipmDcmBiasEstimator() {}

void LipmDcmBiasEstimator::setLipmNaturalFrequency(double omega_0)
{
  omega0_ = omega_0;
  updateMatricesABQ_();
}

void LipmDcmBiasEstimator::setSamplingTime(double dt)
{
  dt_ = dt;
  updateMatricesABQ_();
}

void LipmDcmBiasEstimator::setBias(double bias)
{
  Vector2 x = filter_.getCurrentEstimatedState();
  /// update the bias
  x(1) = bias;
  filter_.setCurrentState(x);
}

void LipmDcmBiasEstimator::setBias(double bias, double uncertainty)
{
  setBias(bias);
  Matrix2 P = filter_.getStateCovariance();
  /// resetting the non diagonal parts
  P(0, 1) = P(1, 0) = 0;
  P(1, 1) = square(uncertainty);
  filter_.setStateCovariance(P);
}

void LipmDcmBiasEstimator::setBiasDriftPerSecond(double driftPerSecond)
{
  Matrix2 Q = filter_.getProcessCovariance();
  /// update the corresponding part in the process noise matrix
  Q(1, 1) = square(driftPerSecond);
  filter_.setProcessCovariance(Q);
}

void LipmDcmBiasEstimator::setDCM(double dcm)
{
  Vector2 x = filter_.getCurrentEstimatedState();
  /// update the bias
  x(0) = dcm;
  filter_.setCurrentState(x);
}

void LipmDcmBiasEstimator::setDCM(double dcm, double uncertainty)
{
  setDCM(dcm);
  Matrix2 P = filter_.getStateCovariance();
  /// resetting the non diagonal parts
  P(0, 1) = P(1, 0) = 0;
  P(0, 0) = square(uncertainty);
  filter_.setStateCovariance(P);
}

void LipmDcmBiasEstimator::setZmpMeasureErrorStd(double std)
{
  zmpErrorStd_ = std;
  updateMatricesABQ_();
}

void LipmDcmBiasEstimator::setDcmMeasureErrorStd(double std)
{
  Matrix1 R;
  R(0, 0) = square(std);
}

void LipmDcmBiasEstimator::setInputs(double dcm, double zmp)
{
  Vector1 u;
  Vector1 y;

  y(0) = dcm;

  /// The prediction of the state depends on the previous value of the ZMP
  u(0) = previousZmp_;
  previousZmp_ = zmp;

  filter_.pushMeasurement(y);
  filter_.pushInput(u);
}

double LipmDcmBiasEstimator::getUnbiasedDCM() const
{
  return filter_.getCurrentEstimatedState()(0);
}

double LipmDcmBiasEstimator::getBias() const
{
  return filter_.getCurrentEstimatedState()(1);
}

void LipmDcmBiasEstimator::LipmDcmBiasEstimator::updateMatricesABQ_()
{
  // clang-format off
  A_ << 1 + omega0_ * dt_, 0,
        0,                 1;

  B_ << -omega0_ * dt_,
        0;

  Q_ << square(omega0_* dt_ * zmpErrorStd_), 0,
        0,                                   square(biasDriftStd_*dt_);
  // clang-format on

  filter_.setA(A_);
  filter_.setB(B_);
  filter_.setProcessCovariance(Q_);
}

} // namespace stateObservation
