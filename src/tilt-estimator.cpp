#include <state-observation/observer/tilt-estimator.hpp>

namespace stateObservation
{

  TiltEstimator::TiltEstimator(double alpha, double beta, double gamma)
    : ZeroDelayObserver(9, 6), alpha_(alpha), beta_(beta), gamma_(gamma),
      p_S_C(Vector3::Zero()), R_S_C(Matrix3::Identity()), v_S_C(Vector3::Zero()), w_S_C(Vector3::Zero()), v_C(Vector3::Zero())
  {
  }
  
  void TiltEstimator::setMeasurement(const Vector3 ya_k, const Vector3 yg_k, TimeIndex k)
  {
    ObserverBase::MeasureVector y_k(6);
    y_k << ya_k, yg_k;

    ZeroDelayObserver::setMeasurement(y_k, k);
  }

  ObserverBase::StateVector TiltEstimator::oneStepEstimation_()
  {
    TimeIndex k = this->x_.getTime();
    
    BOOST_ASSERT(this->y_.size() > 0 && this->y_.checkIndex(k+1) && "ERROR: The measurement vector is not set");

    Vector3 ya = getMeasurement(k+1).head<3>();
    Vector3 yg = getMeasurement(k+1).tail<3>();
    
    Vector3 x1 = R_S_C.transpose() * (v_C + v_S_C) + (yg - R_S_C.transpose() * w_S_C).cross(R_S_C.transpose() * p_S_C);

    ObserverBase::StateVector x_hat = getEstimatedState(k);
    Vector3 x1_hat = x_hat.segment<3>(0);
    Vector3 x2_hat_prime = x_hat.segment<3>(3);
    Vector3 x2_hat = x_hat.segment<3>(6);
    
    Vector dx_hat(9);
    dx_hat.segment<3>(0) = x1_hat.cross(yg) - cst::gravityConstant * x2_hat_prime + ya + alpha_ * (x1 - x1_hat);
    dx_hat.segment<3>(3) = x2_hat_prime.cross(yg) - beta_ * (x1 - x1_hat);
    dx_hat.segment<3>(6) = x2_hat.cross(yg - gamma_ * x2_hat.cross(x2_hat_prime));
    
    x_hat += dx_hat * dt_;

    x_hat.tail<3>() /= x_hat.tail<3>().norm();
    
    setState(x_hat, k+1);
    
    return x_hat;
  }

  
}
