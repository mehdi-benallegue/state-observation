#include <state-observation/dynamics-estimators/lipm-dcm-bias-estimator.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>

namespace stateObservation
{

constexpr double LipmDcmBiasEstimator::defaultDCMUncertainty;
constexpr double LipmDcmBiasEstimator::defaultBiasUncertainty;

constexpr double LipmDcmBiasEstimator::defaultDt_;

/// default expected drift of the bias every second
constexpr double LipmDcmBiasEstimator::defaultBiasDriftSecond_;

/// default error in the estimation of the sensors
constexpr double LipmDcmBiasEstimator::defaultZmpErrorStd_;
constexpr double LipmDcmBiasEstimator::defaultDcmErrorStd_;

using namespace tools;
LipmDcmBiasEstimator::LipmDcmBiasEstimator(double dt,
                                           double omega_0,
                                           double biasDriftStd,
                                           const Vector2 & initZMP,
                                           const Vector2 & initDcm,
                                           const Vector2 & initBias,
                                           double dcmMeasureErrorStd,
                                           double zmpMeasureErrorStd,
                                           const Vector2 & initDcmUncertainty,
                                           const Vector2 & initBiasUncertainty)
: omega0_(omega_0), dt_(dt), biasDriftStd_(biasDriftStd), zmpErrorStd_(zmpMeasureErrorStd), previousZmp_(initZMP),
  filter_(4, 2, 2), A_(Matrix4::Identity()), previousOrientation_(Matrix2::Identity())
{
  updateMatricesABQ_();
  C_ << Matrix2::Identity(), Matrix2::Identity();
  R_ = dblToSqDiag(dcmMeasureErrorStd);
  filter_.setC(C_);
  filter_.setMeasurementCovariance(R_);
  Vector4 x;
  x << initDcm, initBias;
  filter_.setState(x, 0);
  Matrix4 P;
  // clang-format off
  P<< Vec2ToSqDiag(initDcmUncertainty),  Matrix2::Zero(),
      Matrix2::Zero(),                   Vec2ToSqDiag(initBiasUncertainty);
  // clang-format on
  filter_.setStateCovariance(P);
}

void LipmDcmBiasEstimator::resetWithMeasurements(const Vector2 & measuredDcm,
                                                 const Vector2 & measuredZMP,
                                                 const Matrix2 & yaw,
                                                 bool measurementIsWithBias,
                                                 double biasDriftPerSecondStd,
                                                 double dcmMeasureErrorStd,
                                                 double zmpMeasureErrorStd,
                                                 const Vector2 & initBias,
                                                 const Vector2 & initBiasuncertainty)

{
  filter_.reset();
  biasDriftStd_ = biasDriftPerSecondStd;
  zmpErrorStd_ = zmpMeasureErrorStd;
  previousZmp_ = measuredZMP;
  previousOrientation_ = yaw;
  A_.setIdentity();

  updateMatricesABQ_();
  C_ << Matrix2::Identity(), Matrix2::Identity();
  R_ = dblToSqDiag(dcmMeasureErrorStd);
  filter_.setC(C_);
  filter_.setMeasurementCovariance(R_);
  Vector4 x;

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
  Matrix4 P;

  if(measurementIsWithBias)
  {
    Matrix2 initBiasCov = Vec2ToSqDiag(initBiasuncertainty);
    /// The state and the
    // clang-format off
    P<< initBiasCov+dblToSqDiag(dcmMeasureErrorStd), -initBiasCov,
       -initBiasCov,                                  initBiasCov;
    // clang-format on
  }
  else
  {
    // clang-format off
    P<< dblToSqDiag(dcmMeasureErrorStd),  Matrix2::Zero(),
        Matrix2::Zero(),                  Vec2ToSqDiag(initBiasuncertainty);
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

void LipmDcmBiasEstimator::setBias(const Vector2 & bias)
{
  Vector4 x = filter_.getCurrentEstimatedState();
  /// update the bias
  x.tail<2>() = bias;
  filter_.setCurrentState(x);
}

void LipmDcmBiasEstimator::setBias(const Vector2 & bias, const Vector2 & uncertainty)
{
  setBias(bias);
  Matrix4 P = filter_.getStateCovariance();
  /// resetting the non diagonal parts
  P.topRightCorner<2, 2>().setZero();
  P.bottomLeftCorner<2, 2>().setZero();
  P.bottomRightCorner<2, 2>() = Vec2ToSqDiag(uncertainty);
  filter_.setStateCovariance(P);
}

void LipmDcmBiasEstimator::setBiasDriftPerSecond(double driftPerSecond)
{
  /// update the corresponding part in the process noise matrix
  Q_.bottomRightCorner<2, 2>() = dblToSqDiag(driftPerSecond);
  filter_.setProcessCovariance(Q_);
}

void LipmDcmBiasEstimator::setDCM(const Vector2 & dcm)
{
  Vector4 x = filter_.getCurrentEstimatedState();
  /// update the bias
  x.head<2>() = dcm;
  filter_.setCurrentState(x);
}

void LipmDcmBiasEstimator::setDCM(const Vector2 & dcm, const Vector2 & uncertainty)
{
  setDCM(dcm);
  Matrix4 P = filter_.getStateCovariance();
  /// resetting the non diagonal parts
  P.topRightCorner<2, 2>().setZero();
  P.bottomLeftCorner<2, 2>().setZero();
  P.topLeftCorner<2, 2>() = Vec2ToSqDiag(uncertainty);
  filter_.setStateCovariance(P);
}

void LipmDcmBiasEstimator::setZmpMeasureErrorStd(double std)
{
  zmpErrorStd_ = std;
  updateMatricesABQ_();
}

void LipmDcmBiasEstimator::setDcmMeasureErrorStd(double std)
{
  Matrix2 R;
  R = dblToSqDiag(std);
}

void LipmDcmBiasEstimator::setInputs(const Vector2 & dcm, const Vector2 & zmp)
{
  Vector2 u;
  Vector2 y;

  y = dcm;

  /// The prediction of the state depends on the previous value of the ZMP
  u = previousZmp_;
  previousZmp_ = zmp;

  filter_.pushMeasurement(y);
  filter_.pushInput(u);
}

void LipmDcmBiasEstimator::setInputs(const Vector2 & dcm, const Vector2 & zmp, const Matrix2 & orientation)
{
  setInputs(dcm, zmp);
  A_.bottomRightCorner<2, 2>() =
      orientation * previousOrientation_.transpose() * A_.bottomRightCorner<2, 2>(); /// add a rotation
  previousOrientation_ = orientation;
  filter_.setA(A_);
}

Vector2 LipmDcmBiasEstimator::getUnbiasedDCM() const
{
  return filter_.getCurrentEstimatedState().head<2>();
}

Vector2 LipmDcmBiasEstimator::getBias() const
{
  return filter_.getCurrentEstimatedState().tail<2>();
}

void LipmDcmBiasEstimator::LipmDcmBiasEstimator::updateMatricesABQ_()
{
  // clang-format off

  ///We only modify a corner to avoid resetting the orientation
  A_.topLeftCorner<2,2>() = dblToDiag(1 + omega0_ * dt_);

  B_ << dblToDiag(-omega0_ * dt_),
        Matrix2::Zero();

  Q_ << dblToSqDiag(omega0_* dt_ * zmpErrorStd_), Matrix2::Zero(),
         Matrix2::Zero(),                         dblToSqDiag(biasDriftStd_*dt_);
  // clang-format on

  filter_.setA(A_);
  filter_.setB(B_);
  filter_.setProcessCovariance(Q_);
}

} // namespace stateObservation
