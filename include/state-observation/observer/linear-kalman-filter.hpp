/**
 * \file      linear-kalman-filter.hpp
 * \author    Mehdi Benallegue
 * \date       2012
 * \brief      Defines the class of a Linear Kalman filter
 *
 *              It implements the Kalman filter for linear systems (LTI-LTV).
 *              This is the class to instanciate when you want to use Kalman filtering
 *              for linear systems
 *
 *             x_{k+1}=A_k x_k+ B_k u_k + v_k
 *
 *             y_k=C_k x_k + D_k u_k + w_k
 *
 * \details
 *
 *
 */

#ifndef STATEOBSERVER_LINEARKALMANFILTERHPP
#define STATEOBSERVER_LINEARKALMANFILTERHPP

#include <state-observation/api.h>
#include <state-observation/observer/kalman-filter-base.hpp>

namespace stateObservation
{

/**
 * \class  LinearKalmanFilter
 * \brief
 *        The class of a Linear Kalman filter
 *
 *        It implements the Kalman filter for linear systems (LTI-LTV).
 *        This is the class to instanciate when you want to use Kalman filtering
 *        for linear systems. To use this class, one needs to provide the
 *        matrices that describe the dynamics of the state and the measurement.
 *
 *             x_{k+1}=A_k x_k+ B_k u_k + v_k
 *
 *             y_k=C_k x_k + D_k u_k + w_k
 *
 *
 *
 */

class STATE_OBSERVATION_DLLAPI LinearKalmanFilter : public KalmanFilterBase
{
public:
  /// The constructor
  ///  \li n : size of the state vector
  ///  \li m : size of the measurements vector
  ///  \li p : size of the input vector
  LinearKalmanFilter(Index n, Index m, Index p = 0) : KalmanFilterBase(n, m, p) {}

  /// Default constructor
  LinearKalmanFilter() {}

  /// The type of the matrix linking the input to the state
  typedef Matrix Bmatrix;

  /// The type of the matrix linking the input to the measurement
  typedef Matrix Dmatrix;

  /// Set the value of the input-state matrix
  virtual void setB(const Bmatrix & B);

  /// Clear the value of the input-state Matrix
  virtual void clearB();

  /// Set the value of the input-measurement matrix
  virtual void setD(const Dmatrix & D);

  /// Clear the value of the input-measurement matrix
  virtual void clearD();

  /// Get a matrix having the size of the B matrix having "c" values
  Bmatrix getBmatrixConstant(double c) const;

  /// Get a matrix having the size of the B matrix having random values
  Bmatrix getBmatrixRandom() const;

  /// Get a matrix having the size of the B matrix having zero values
  Bmatrix getBmatrixZero() const;

  /// checks whether or not a matrix has the dimensions of the B matrix
  bool checkBmatrix(const Bmatrix &) const;

  /// Get a matrix having the size of the D matrix having "c" values
  Dmatrix getDmatrixConstant(double c) const;

  /// Get a matrix having the size of the D matrix having random values
  Dmatrix getDmatrixRandom() const;

  /// Get a matrix having the size of the D matrix having zero values
  Dmatrix getDmatrixZero() const;

  /// checks whether or not a matrix has the dimensions of the D matrix
  bool checkDmatrix(const Dmatrix &) const;

  /// changes the dimension of the state vector:
  /// resets the internal container for the state vector and
  /// the containers for the matrices A, B, C, Q, P
  virtual void setStateSize(Index n);

  /// changes the dimension of the measurement vector:
  /// resets the internal container for the measurement vectors and
  /// the containers for the matrices C, D, R
  virtual void setMeasureSize(Index m);

  /// changes the dimension of the input vector:
  /// resets the internal container for the input vectors and
  /// the containers for the matrices B, D
  virtual void setInputSize(Index p);

protected:
  /// The implementation of the (linear) prediction (state dynamics)
  virtual StateVector prediction_(TimeIndex k);

  /// The implementation of the (linear) measurement (state dynamics)
  virtual MeasureVector simulateSensor_(const StateVector & x, TimeIndex k);

  /// The container of the Input-State matrix
  Matrix d_;

  /// The container of the Input-Measurement matrix
  Matrix b_;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // namespace stateObservation

#endif // STATEOBSERVER_LINEARKALMANFILTERHPP
