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
  *         use getEstimatedState to obtain the state vector
  *         the tilt R.transpose()*e_z is constituted
  *         with the last three components of the state vector.
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

    ///set the gain of x1_hat variable
    void setAlpha(const double alpha) { alpha_ = alpha; }
    double getAlpha() const { return alpha_; }
    
    ///set the gain of x2prime_hat variable
    void setBeta(const double beta) { beta_ = beta; }
    double getBeta() const { return beta_; }

    ///set the gain of x2_hat variable
    void setGamma(const double gamma) { gamma_ = gamma; }
    double getGamma() const { return gamma_; }
    
    ///set the sampling time of the measurements
    void setSamplingTime(const double dt) { dt_ = dt; }
    double getSamplingTime() const { return dt_; }

    /// sets the position of the IMU sensor in the control frame
    void setSensorPositionInC(const Vector3& p) { p_S_C_ = p; }
    Vector3 getSensorPositionInC() { return p_S_C_; }
    
    /// sets the oriantation of the IMU sensor in the control frame
    void setSensorOrientationInC(const Matrix3& R) { R_S_C_ = R; }
    Matrix3 getSensorOrientationInC() { return R_S_C_; }
    
    /// sets teh linear velocity of the IMU sensor in the control frame
    void setSensorLinearVelocityInC(const Vector3& v) { v_S_C_ = v; }
    Vector3 getSensorLinearVelocityInC() { return v_S_C_; }
    
    /// sets the angular velocity of the IMU sensor in the control frame
    void setSensorAngularVelocityInC(const Vector3& w) { w_S_C_ = w; }
    Vector3 getSensorAngularVelocityInC() { return w_S_C_; }

    /// sets the velocity of the control origin in the world frame
    void setControlOriginVelocityInW(const Vector3& v) { v_C_ = v; }
    Vector3 getControlOriginVelocityInW() { return v_C_; }
    
    /// sets ths measurement (accelero and gyro stacked in one vector)
    void setMeasurement(const Vector3 ya_k, const Vector3 yg_k, TimeIndex k);

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
  protected:

    /// The parameters of the estimator
    double alpha_, beta_, gamma_;

    /// Sampling time
    double dt_;
    
    /// Position of the IMU in the control frame
    Vector3 p_S_C_;

    /// Orientation of the IMU in the control frame
    Matrix3 R_S_C_;

    /// Linear velocity of the IMU in the control frame
    Vector3 v_S_C_;

    /// Angular velocity of the IMU in the control frame
    Vector3 w_S_C_;

    /// Linear velocity of the control frame
    Vector3 v_C_;


    ///variables used for the computation
    Vector3 x1_;
    Vector3 x1_hat_;
    Vector3 x2_hat_prime_;
    Vector3 x2_hat_;
    Vector3 dx1_hat;
    
    

    
    /// The tilt estimator loop
    StateVector oneStepEstimation_();
  };
  
}

#endif //TILTESTIMATORHPP
