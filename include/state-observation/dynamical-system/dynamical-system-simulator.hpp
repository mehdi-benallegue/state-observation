/**
 * \file      dynamical-system-simulator.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief      Provides an interface to simulate the dynamics provided by a
 *             dynamics functor.
 *
 * \details
 *
 *
 */

#ifndef STATEOBSERVATIONDYNAMICALSYSTEMSIMULATOR_H
#define STATEOBSERVATIONDYNAMICALSYSTEMSIMULATOR_H

#include <map>

#include <state-observation/api.h>
#include <state-observation/dynamical-system/dynamical-system-functor-base.hpp>

namespace stateObservation
{
/**
 * \class  DynamicalSystemSimulator
 * \brief
 *        The class gives a small encapsulation of the dynamics functor,
 *       which enables the simulation of the dynamics and the storage of
 *       states, inputs and measurements.
 *
 */
class STATE_OBSERVATION_DLLAPI DynamicalSystemSimulator
{
public:
  /// Constructor
  DynamicalSystemSimulator();

  /// Virtual destructor
  virtual ~DynamicalSystemSimulator();

  /// Sets a pointer to the dynamics Functor (the class does not destroy
  /// the Functor when the class is destroyed.)
  virtual void setDynamicsFunctor(DynamicalSystemFunctorBase *);

  /// Sets the state value at instant k
  virtual void setState(const Vector & x, TimeIndex k);

  /// Set the value of the input at instant k
  /// If no value of the input is available when simulating the system
  /// the previous input value is automatically set.
  virtual void setInput(const Vector & u, TimeIndex k);

  /// Gets the state of the current time
  virtual Vector getCurrentState() const;

  /// Gets the current time
  virtual TimeIndex getCurrentTime() const;

  /// Runs one loop of the dynamics simulation
  virtual void simulateDynamics();

  /// Runs the simulation until a given time index
  virtual void simulateDynamicsTo(TimeIndex k);

  /// Gets the input of a given time index,
  /// if the Input is not available the previous input is provided
  virtual Vector getInput(TimeIndex k) const;

  /// Gets the measurement value at a given time index
  /// if the value is not available, the dynamics is simulated
  /// to the time k
  virtual Vector getMeasurement(TimeIndex k);

  /// Gets the state value at a given time index
  /// if the value is not available, the dynamics is simulated
  /// to the time k
  virtual Vector getState(TimeIndex k);

  /// Gives a IndexedVectorArray of the measurements starting at startingTime
  /// and having the given  duration
  virtual IndexedVectorArray getMeasurementArray(TimeIndex startingTime, TimeSize duration);

  /// Gives a IndexedVectorArray of the states starting at startingTime
  /// and having the given  duration
  virtual IndexedVectorArray getStateArray(TimeIndex startingTime, TimeSize duration);

  /// resets all the states, the measurements and the inputs
  virtual void resetDynamics();

  /// resets all the simulator with even the dynamics functor
  virtual void resetSimulator();

protected:
  DynamicalSystemFunctorBase * f_;

  IndexedVectorArray x_;

  IndexedVectorArray y_;

  std::map<TimeIndex, Vector> u_;
};
} // namespace stateObservation
#endif // STATEOBSERVATIONDYNAMICALSYSTEMSIMULATOR_H
