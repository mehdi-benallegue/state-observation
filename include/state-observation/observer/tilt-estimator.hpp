/**
 * \file      tilt-estimator.hpp
 * \author    Rafael Cisneros, Mehdi Benallegue
 * \date       2018
 * \brief      Defines the class for the tilt estimator.
 *
 * \details
 *
 *
 */



#ifndef TILTESTIMATORHPP
#define TILTESTIMATORHPP

#include <state-observation/observer/zero-delay-observer.hpp>


namespace stateObservation
{

/**
  * \class  TiltEstimator
  * \brief
  *         Description is pending
  *
  *
  *
  */
  class TiltEstimator: public ZeroDelayObserver
  {
  public:

    /// The constructor
    ///  \li alpha : parameter related to the convergence of the linear velocity
    ///              of the IMU expressed in the control frame
    ///  \li beta  : parameter related to the fast convergence of the tilt
    ///  \li gamma : parameter related to the orthogonality
    TiltEstimator(double alpha, double beta, double gamma);

    void setAlpha(const double alpha) { alpha_ = alpha; }
    double getAlpha() const { return alpha_; }
    
    void setBeta(const double beta) { beta_ = beta; }
    double getBeta() const { return beta_; }

    void setGamma(const double gamma) { gamma_ = gamma; }
    double getGamma() const { return gamma_; }
    
    void setSamplingTime(const double dt) { dt_ = dt; }
    double getSamplingTime() const { return dt_; }

    void setSensorPositionInC(const Vector3& p) { p_S_C = p; }
    Vector3 getSensorPositionInC() { return p_S_C; }
    
    void setSensorOrientationInC(const Matrix3& R) { R_S_C = R; }
    Matrix3 getSensorOrientationInC() { return R_S_C; }
    
    void setSensorLinearVelocityInC(const Vector3& v) { v_S_C = v; }
    Vector3 getSensorLinearVelocityInC() { return v_S_C; }
    
    void setSensorAngularVelocityInC(const Vector3& w) { w_S_C = w; }
    Vector3 getSensorAngularVelocityInC() { return w_S_C; }

    void setControlOriginVelocityInW(const Vector3& v) { v_C = v; }
    Vector3 getControlOriginVelocityInW() { return v_C; }
    
    void setMeasurement(const Vector3 ya_k, const Vector3 yg_k, TimeIndex k);
    
  protected:

    /// The parameters of the estimator
    double alpha_, beta_, gamma_;

    /// Sampling time
    double dt_;
    
    /// Position of the IMU in the control frame
    Vector3 p_S_C;

    /// Orientation of the IMU in the control frame
    Matrix3 R_S_C;

    /// Linear velocity of the IMU in the control frame
    Vector3 v_S_C;

    /// Angular velocity of the IMU in the control frame
    Vector3 w_S_C;

    /// Linear velocity of the control frame
    Vector3 v_C;
    
    /// The tilt estimator loop
    StateVector oneStepEstimation_();
  };
  
}

#endif //TILTESTIMATORHPP
