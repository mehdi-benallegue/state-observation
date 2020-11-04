#include <state-observation/dynamical-system/bidim-elastic-inv-pendulum-dyn-sys.hpp>

namespace stateObservation
{

using stateObservation::cst::gravityConstant;

BidimElasticInvPendulum::BidimElasticInvPendulum() : m_(1), h_(1), processNoise_(0x0), dt_(1)
{
#ifdef STATEOBSERVATION_VERBOUS_CONSTRUCTORS
  std::cout << std::endl << "IMUFixedContactDynamicalSystem Constructor" << std::endl;
#endif // STATEOBSERVATION_VERBOUS_CONSTRUCTOR
       // ctor
}

BidimElasticInvPendulum::~BidimElasticInvPendulum()
{
  // dtor
}

Vector BidimElasticInvPendulum::stateDynamics(const Vector & x, const Vector & u, TimeIndex)
{
  assertStateVector_(x);
  assertInputVector_(u);

  // x_{k+1}
  Vector xk1 = Vector::Zero(stateSize_, 1);

  double p = x[0];
  double pdot = x[2];
  double pdotdot = u[0];
  double theta = x[1];
  double thetadot = x[3];
  double thetadotdot = (1 / (m_ * (p * p + h_ * h_)))
                       * (-k_ * theta - m_ * gravityConstant * (cos(theta) * p - sin(theta) * h_)
                          + m_ * (h_ * pdotdot - 2 * thetadot * (h_ * pdot)));

  xk1[0] = p + dt_ * pdot;
  xk1[1] = theta + dt_ * thetadot;
  xk1[2] = pdot + dt_ * pdotdot;
  xk1[3] = thetadot + dt_ * thetadotdot;

  if(processNoise_ != 0x0)
    return processNoise_->getNoisy(xk1);
  else
    return xk1;
}

/// set the height of the com of the pendulum
void BidimElasticInvPendulum::setHeight(const double & h)
{
  h_ = h;
}

/// set the mass of the pendulum
void BidimElasticInvPendulum::setMass(const double & m)
{
  m_ = m;
}

/// set the elasticity of the pendulum
void BidimElasticInvPendulum::setElasticity(const double & k)
{
  k_ = k;
}

Vector BidimElasticInvPendulum::measureDynamics(const Vector &, const Vector &, TimeIndex)
{
  /// There is no measurements
  return Vector::Zero(0);
}

void BidimElasticInvPendulum::setProcessNoise(NoiseBase * n)
{
  processNoise_ = n;
}

void BidimElasticInvPendulum::resetProcessNoise()
{
  processNoise_ = 0x0;
}

void BidimElasticInvPendulum::setMeasurementNoise(NoiseBase *)
{
  // there is no measurement
}

void BidimElasticInvPendulum::resetMeasurementNoise()
{
  // there is no measurement
}

void BidimElasticInvPendulum::setSamplingPeriod(double dt)
{
  dt_ = dt;
}

Index BidimElasticInvPendulum::getStateSize() const
{
  return stateSize_;
}

Index BidimElasticInvPendulum::getInputSize() const
{
  return inputSize_;
}

Index BidimElasticInvPendulum::getMeasurementSize() const
{
  return measurementSize_;
}

NoiseBase * BidimElasticInvPendulum::getProcessNoise() const
{
  return processNoise_;
}

NoiseBase * BidimElasticInvPendulum::getMeasurementNoise() const
{
  return 0x0;
}
} // namespace stateObservation
