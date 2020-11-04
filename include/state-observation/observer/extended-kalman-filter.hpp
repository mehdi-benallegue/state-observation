/**
 * \file      extended-kalman-filter-base.hpp
 * \author    Mehdi Benallegue
 * \date       2012
 * \brief      Defined the class to intanciate to use an extended Kalman filter.
 *
 *             x_{k+1}=f(x_k,u_k)+v_k
 *             y_k=h(x_k,u_k)+w_k
 *
 * \details
 *
 *
 */

#ifndef STATEOBSERVER_EXTENDEDKALMANFILTERHPP
#define STATEOBSERVER_EXTENDEDKALMANFILTERHPP

#include <state-observation/api.h>
#include <state-observation/dynamical-system/dynamical-system-functor-base.hpp>
#include <state-observation/observer/kalman-filter-base.hpp>

namespace stateObservation
{
/**
 * \class  ExtendedKalmanFilter
 * \brief
 *
 *        The class to intanciate to use an extended Kalman filter.
 *        To use this class, one needs to provide a pointer on a functor
 *        that describes the state dynamics and the measurement dynamics.
 *        The functor type needs to be a derived class from the class
 *        DynamicsFunctorBase.
 *
 *        x_{k+1}=f(x_k,u_k)+v_k
 *
 *        y_k=h(x_k,u_k)+w_k
 *
 *
 *
 */

class STATE_OBSERVATION_DLLAPI ExtendedKalmanFilter : public KalmanFilterBase
{

public:
  /// The constructor.
  ///  \li n : size of the state vector
  ///  \li m : size of the measurements vector

  ExtendedKalmanFilter(Index n, Index m);

  /// The constructor.
  ///  \li n : size of the state vector
  ///  \li m : size of the measurements vector
  ///  \li p : size of the input vector
  ///  \li The parameter directInputOutputFeedthrough defines whether (true) or not (false) the measurement y_k requires
  ///  the input u_k \li The parameter directInputStateProcessFeedthrough defines whether (true) or not (false) the
  ///  state x_{k+1} requires the input u_k

  ExtendedKalmanFilter(Index n,
                       Index m,
                       Index p,
                       bool directInputOutputFeedthrough = true,
                       bool directInputStateProcessFeedthrough = true);

  /// The constructor.
  ///  \li n : size of the state vector
  ///  \li nt: size of the state tengent vector representation (usually nt<=n)
  ///  \li m : size of the measurements vector
  ///  \li p : size of the input vector
  ///  \li The parameter directInputOutputFeedthrough defines whether (true) or not (false) the measurement y_k requires
  ///  the input u_k \li The parameter directInputStateProcessFeedthrough defines whether (true) or not (false) the
  ///  state x_{k+1} requires the input u_k

  ExtendedKalmanFilter(Index n,
                       Index nt,
                       Index m,
                       Index mt,
                       Index p,
                       bool directInputOutputFeedthrough,
                       bool directInputStateProcessFeedthrough);

  /// Set a pointer to the functor that defines the dynamics of the states
  /// and the measurement the user is responsible for the validity of the
  /// pointer during the execution of the kalman filter
  void setFunctor(DynamicalSystemFunctorBase * f);
  DynamicalSystemFunctorBase * getFunctor(void) const;

  /// Gets a pointer to the functor
  DynamicalSystemFunctorBase * functor() const;

  /// Clear the value of the functor
  /// Does not destroy the pointed object
  void clearFunctor();

  /// Precise whether (true) or not (false) the measurement y_k requires
  /// the input u_k
  void setDirectInputOutputFeedthrough(bool b = true);

  /// Precise whether (true) or not (false) the estimation of the state x_{k+1} requires
  /// the input u_k
  void setDirectInputStateFeedthrough(bool b = true);

  /// Give an estimation of A matrix using
  /// finite difference method (the forward difference method)
  /// the parameter dx is the step vector for derivation
  /// it must have the size nt (tangent vector)
  virtual Amatrix getAMatrixFD(const Vector & dx);

  /// Give an estimation of C matrix using
  /// finite difference method (the forward difference method)
  /// the parameter dx is the step vector for derivation
  /// it must have the size nt (tangent vector)
  virtual Cmatrix getCMatrixFD(const Vector & dx);

  /// Reset the extended kalman filter (call also the reset function of the dynamics functor)
  virtual void reset();

protected:
  /// simulate the dynamics of the state using the functor
  virtual StateVector prediction_(TimeIndex k);

  /// simulate the dynamic of the measurement using the functor
  virtual MeasureVector simulateSensor_(const StateVector & x, TimeIndex k);

  /// predicts the measurement using the functor, assumed that the predicted state is up-to-date
  virtual MeasureVector predictSensor_(TimeIndex k);

  /// boolean that provides if theris a need of not for input for the masurement
  bool directInputOutputFeedthrough_;

  /// boolean that provides if theris a need of not for input for the state dynamics
  bool directInputStateProcessFeedthrough_;

  /// pointer on the dynamics functor
  DynamicalSystemFunctorBase * f_;

  // optimization
  struct Optimization
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    ObserverBase::InputVector u_;
    KalmanFilterBase::Amatrix a_;
    KalmanFilterBase::Cmatrix c_;
    ObserverBase::StateVector x_;
    ObserverBase::StateVector dx_;
    ObserverBase::StateVector xp_;
    ObserverBase::MeasureVector yp_;

  } opt;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // namespace stateObservation

#endif // EXTENDEDKALMANFILTERHPP
