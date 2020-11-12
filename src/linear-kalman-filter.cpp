#include <state-observation/observer/linear-kalman-filter.hpp>

namespace stateObservation
{

void LinearKalmanFilter::setB(const Bmatrix & B)
{
  BOOST_ASSERT(checkBmatrix(B) && "ERROR: The B matrix size is incorrect");
  b_ = B;
}

void LinearKalmanFilter::clearB()
{
  b_.resize(0, 0);
}

void LinearKalmanFilter::setD(const Dmatrix & D)
{
  BOOST_ASSERT(checkDmatrix(D) && "ERROR: The D matrix size is incorrect");

  d_ = D;
}

void LinearKalmanFilter::clearD()
{
  d_.resize(0, 0);
}

ObserverBase::StateVector LinearKalmanFilter::prediction_(TimeIndex k)
{
  (void)k; // unused

  BOOST_ASSERT(checkAmatrix(a_) && "ERROR: The A matrix is not initialized");
  BOOST_ASSERT(checkBmatrix(b_) && "ERROR: The B matrix is not initialized");
  BOOST_ASSERT(checkCmatrix(c_) && "ERROR: The C matrix is not initialized");

  xbar_.set(a_ * x_(), k);

  if(p_ > 0 && b_ != getBmatrixZero())
  {
    BOOST_ASSERT(u_.checkIndex(k - 1)
                 && "ERROR: The input feedthrough of the state dynamics is not set "
                    "(the state at time k+1 needs the input at time k which was not given) "
                    "if you don't need the input in the computation of state, you "
                    "must set B matrix to zero");
    xbar_().noalias() += b_ * this->u_[k - 1];
  }

  return xbar_();
}

ObserverBase::MeasureVector LinearKalmanFilter::simulateSensor_(const StateVector & x, TimeIndex k)
{

  BOOST_ASSERT(checkCmatrix(c_) && "ERROR: The C matrix is not initialized");

  if(p_ > 0 && checkDmatrix(d_) && d_ != getDmatrixZero())
  {
    BOOST_ASSERT(u_.checkIndex(k)
                 && "ERROR: The input feedthrough of the measurements is not set "
                    "(the measurement at time k needs the input at time k which was not given) "
                    "if you don't need the input in the computation of measurement, you "
                    "must set D matrix to zero");
    ybar_.set(c_ * x + d_ * u_[k], k);
  }
  else
  {
    ybar_.set(c_ * x, k);
  }
  return ybar_();
}

void LinearKalmanFilter::reset()
{
  KalmanFilterBase::reset();

  clearB();
  clearD();
}

LinearKalmanFilter::Bmatrix LinearKalmanFilter::getBmatrixConstant(double c) const
{
  return Bmatrix::Constant(n_, p_, c);
}

LinearKalmanFilter::Bmatrix LinearKalmanFilter::getBmatrixRandom() const
{
  return Bmatrix::Random(n_, p_);
}

LinearKalmanFilter::Bmatrix LinearKalmanFilter::getBmatrixZero() const
{
  return Bmatrix::Zero(n_, p_);
}

bool LinearKalmanFilter::checkBmatrix(const Bmatrix & a) const
{
  return (a.rows() == n_ && a.cols() == p_);
}

LinearKalmanFilter::Dmatrix LinearKalmanFilter::getDmatrixConstant(double c) const
{
  return Dmatrix::Constant(m_, p_, c);
}

LinearKalmanFilter::Dmatrix LinearKalmanFilter::getDmatrixRandom() const
{
  return Dmatrix::Random(m_, p_);
}

LinearKalmanFilter::Dmatrix LinearKalmanFilter::getDmatrixZero() const
{
  return Dmatrix::Zero(m_, p_);
}

bool LinearKalmanFilter::checkDmatrix(const Dmatrix & a) const
{
  return (a.rows() == m_ && a.cols() == p_);
}

void LinearKalmanFilter::setStateSize(Index n)
{
  if(n != n_)
  {
    KalmanFilterBase::setStateSize(n);
    clearB();
  }
}

void LinearKalmanFilter::setMeasureSize(Index m)
{
  if(m != m_)
  {
    KalmanFilterBase::setMeasureSize(m);
    clearD();
  }
}

void LinearKalmanFilter::setInputSize(Index p)
{
  if(p != p_)
  {
    KalmanFilterBase::setInputSize(p);
    clearB();
    clearD();
  }
}

} // namespace stateObservation
