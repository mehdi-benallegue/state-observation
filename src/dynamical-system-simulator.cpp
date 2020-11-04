#include <state-observation/dynamical-system/dynamical-system-simulator.hpp>

namespace stateObservation
{
DynamicalSystemSimulator::DynamicalSystemSimulator() : f_(0x0)
{
  // ctor
}

DynamicalSystemSimulator::~DynamicalSystemSimulator()
{
  // dtor
}

void DynamicalSystemSimulator::setDynamicsFunctor(DynamicalSystemFunctorBase * f)
{
  f_ = f;
}

void DynamicalSystemSimulator::setState(const Vector & x, TimeIndex k)
{
  BOOST_ASSERT((x_.size() == 0 || (x_.getFirstIndex() <= k && x_.getNextIndex() >= k))
               && "ERROR: Only consecutive states can be set. If you want to restart a new dynamics please call "
                  "resetDynamics before");
  x_.truncateAfter(k - 1);
  x_.setValue(x, k);
}

void DynamicalSystemSimulator::setInput(const Vector & u, TimeIndex k)
{
  u_.insert(std::pair<TimeIndex, Vector>(k, u));
}

Vector DynamicalSystemSimulator::getMeasurement(TimeIndex k)
{
  BOOST_ASSERT((y_.size() == 0 || y_.getFirstIndex() < k) && "ERROR: Only future measurements can be obtained");
  if(y_.getNextIndex() < k) simulateDynamicsTo(k);

  return y_[k];
}

Vector DynamicalSystemSimulator::getState(TimeIndex k)
{
  BOOST_ASSERT((x_.size() == 0 || x_.getFirstIndex() <= k) && "ERROR: Only future measurements can be obtained");
  if(x_.getLastIndex() < TimeIndex(k)) simulateDynamicsTo(k - 1);

  return x_[k];
}

Vector DynamicalSystemSimulator::getCurrentState() const
{
  return x_[TimeIndex(x_.getLastIndex())];
}

TimeIndex DynamicalSystemSimulator::getCurrentTime() const
{
  return x_.getLastIndex();
}

Vector DynamicalSystemSimulator::getInput(TimeIndex k) const
{
  BOOST_ASSERT(u_.size() > 0 && "ERROR: the input is not set");
  std::map<TimeIndex, Vector>::const_iterator i = u_.upper_bound(k);

  --i;

  BOOST_ASSERT(i->first <= k && "ERROR: the input is not set for the current time");

  return i->second;
}

void DynamicalSystemSimulator::simulateDynamics()
{
  BOOST_ASSERT(f_ != 0x0 && "ERROR: A dynamics functor must be set");

  TimeIndex k(x_.getLastIndex());
  Vector u = getInput(k);
  y_.setValue(f_->measureDynamics(x_[k], u, k), k);
  x_.setValue(f_->stateDynamics(x_[k], u, k), k + 1);
}

void DynamicalSystemSimulator::simulateDynamicsTo(TimeIndex k)
{
  for(TimeIndex i = x_.getLastIndex(); i < k; ++i)
  {
    simulateDynamics();
  }
}

IndexedVectorArray DynamicalSystemSimulator::getMeasurementArray(TimeIndex startingTime, TimeSize duration)
{
  BOOST_ASSERT(startingTime > y_.getFirstIndex() && "ERROR: The starting time is too early, try later starting time");
  IndexedVectorArray a;

  for(TimeIndex i = startingTime; i < startingTime + TimeIndex(duration); ++i)
  {
    a.setValue(getMeasurement(i), i);
  }
  return a;
}

IndexedVectorArray DynamicalSystemSimulator::getStateArray(TimeIndex startingTime, TimeSize duration)
{
  BOOST_ASSERT(startingTime > x_.getFirstIndex() && "ERROR: The starting time is too early, try later starting time");

  IndexedVectorArray a;

  for(TimeIndex i = startingTime; i < startingTime + TimeIndex(duration); ++i)
  {
    a.setValue(getState(i), i);
  }
  return a;
}

void DynamicalSystemSimulator::resetDynamics()
{
  x_.reset();
  y_.reset();
  u_.clear();
}

void DynamicalSystemSimulator::resetSimulator()
{
  resetDynamics();
  f_ = 0x0;
}
} // namespace stateObservation
