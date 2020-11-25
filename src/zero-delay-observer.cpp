#include <state-observation/observer/zero-delay-observer.hpp>

namespace stateObservation
{

void ZeroDelayObserver::setState(const ObserverBase::StateVector & x_k, TimeIndex k)
{
  BOOST_ASSERT(checkStateVector(x_k) && "The size of the state vector is incorrect");

  x_.set(x_k, k);
  while(y_.size() > 0 && y_.getFirstIndex() <= k)
  {
    y_.popFront();
  }

  if(p_ > 0)
    while(u_.size() > 0 && u_.getFirstIndex() < k)
    {
      u_.popFront();
    }
}

void ZeroDelayObserver::setCurrentState(const ObserverBase::StateVector & x_k)
{
  BOOST_ASSERT(x_.isSet() && "The state vector has not been set");
  BOOST_ASSERT(checkStateVector(x_k) && "The size of the state vector is incorrect");

  x_() = x_k;
}

void ZeroDelayObserver::clearStates()
{
  x_.reset();
}

void ZeroDelayObserver::setMeasurement(const ObserverBase::MeasureVector & y_k, TimeIndex k)
{

  BOOST_ASSERT(checkMeasureVector(y_k) && "The size of the measure vector is incorrect");
  if(y_.size() > 0)
    BOOST_ASSERT((y_.getNextIndex() == k || y_.checkIndex(k)) && "ERROR: The time is set incorrectly for \
                                the measurements (order or gap)");
  else
    BOOST_ASSERT((!x_.isSet() || x_.getTime() == k - 1) && "ERROR: The time is set incorrectly for the measurements \
                                (must be [current_time+1])");

  y_.setValue(y_k, k);
}

void ZeroDelayObserver::pushMeasurement(const ObserverBase::MeasureVector & y_k)
{

  BOOST_ASSERT(checkMeasureVector(y_k) && "The size of the measure vector is incorrect");
  if(y_.size() > 0)
  {
    y_.pushBack(y_k);
  }
  else
  {
    BOOST_ASSERT(x_.isSet()
                 && "Unable to initialize measurement without time index, the state vector has not been set.");
    /// we need the measurement of the next state to correct for the prediction
    y_.setValue(y_k, x_.getTime() + 1);
  }
}

void ZeroDelayObserver::clearMeasurements()
{
  y_.reset();
}

void ZeroDelayObserver::setInput(const ObserverBase::InputVector & u_k, TimeIndex k)
{
  if(p_ > 0)
  {
    BOOST_ASSERT(checkInputVector(u_k) && "The size of the input vector is incorrect");

    if(u_.size() > 0)
      BOOST_ASSERT((u_.getNextIndex() == k || u_.checkIndex(k)) && "ERROR: The time is set incorrectly \
                                for the inputs (order or gap)");
    else
    {
      BOOST_ASSERT((!x_.isSet() || x_.getTime() == k || x_.getTime() == k - 1)
                   && "ERROR: The time is set incorrectly for the \
                          inputs (must be [current_time] or [current_time+1])");
    }

    u_.setValue(u_k, k);
  }
}

void ZeroDelayObserver::pushInput(const ObserverBase::InputVector & u_k)
{
  if(p_ > 0)
  {
    BOOST_ASSERT(checkInputVector(u_k) && "The size of the input vector is incorrect");

    if(u_.size() > 0)
    {
      u_.pushBack(u_k);
    }
    else
    {
      BOOST_ASSERT(x_.isSet() && "Unable to initialize input without time index, the state vector has not been set.");
      /// we need the input at the time of the state vector to predict the next one
      u_.setValue(u_k, x_.getTime());
    }
  }
}

void ZeroDelayObserver::clearInputs()
{
  u_.reset();
}

void ZeroDelayObserver::clearInputsAndMeasurements()
{
  u_.reset();
  y_.reset();
}

TimeIndex ZeroDelayObserver::estimateState()
{
  getEstimatedState(getMeasurementTime());
  return getCurrentTime();
}

ObserverBase::StateVector ZeroDelayObserver::getEstimatedState(TimeIndex k)
{
  BOOST_ASSERT(x_.isSet() && "The state vector has not been set");

  TimeIndex k0 = getCurrentTime();

  BOOST_ASSERT(k0 <= k && "ERROR: The observer cannot estimate previous states");

  for(TimeIndex i = k0; i < k; ++i)
  {
    oneStepEstimation_();
    if(y_.getFirstIndex() < k) y_.popFront();

    if(p_ > 0)
      if(u_.getFirstIndex() < k) u_.popFront();
  }

  return x_();
}

ObserverBase::StateVector ZeroDelayObserver::getCurrentEstimatedState() const
{
  BOOST_ASSERT(x_.isSet() && "The state vector has not been set");
  return x_();
}

TimeIndex ZeroDelayObserver::getCurrentTime() const
{
  BOOST_ASSERT(x_.isSet() && "The state vector has not been set");
  return x_.getTime();
}

Vector ZeroDelayObserver::getInput(TimeIndex k) const
{
  return u_[k];
}

TimeSize ZeroDelayObserver::getInputsNumber() const
{
  return u_.size();
}

TimeIndex ZeroDelayObserver::getInputTime() const
{

  BOOST_ASSERT(y_.size() > 0 && "ERROR: There is no measurements registered (past measurements are erased)");
  return u_.getLastIndex();
}

Vector ZeroDelayObserver::getMeasurement(TimeIndex k) const
{
  return y_[k];
}

TimeIndex ZeroDelayObserver::getMeasurementTime() const
{
  BOOST_ASSERT(y_.size() > 0 && "ERROR: There is no measurements registered (past measurements are erased)");
  return y_.getLastIndex();
}

TimeSize ZeroDelayObserver::getMeasurementsNumber() const
{
  return TimeSize(y_.size());
}

void ZeroDelayObserver::setStateSize(Index n)
{
  if(n != n_)
  {
    ObserverBase::setStateSize(n);
    clearStates();
  }
}

void ZeroDelayObserver::setMeasureSize(Index m)
{
  if(m != m_)
  {
    ObserverBase::setMeasureSize(m);
    clearMeasurements();
  }
}

void ZeroDelayObserver::setInputSize(Index p)
{
  if(p != p_)
  {
    ObserverBase::setInputSize(p);
    clearInputs();
  }
}
} // namespace stateObservation
