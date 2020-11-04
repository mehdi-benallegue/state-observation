#include <fstream>
#include <iostream>

#include <state-observation/dynamical-system/dynamical-system-simulator.hpp>
#include <state-observation/examples/offline-ekf-flexibility-estimation.hpp>
#include <state-observation/noise/gaussian-white-noise.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>

using namespace stateObservation;

int test()
{
  /// The number of samples
  const Index kmax = 2000;

  /// sampling period
  const double dt = 5e-3;

  /// Sizes of the states for the state, the measurement, and the input vector
  const unsigned stateSize = 18;
  const unsigned measurementSize = 6;
  const unsigned inputSize = 15;
  (void)measurementSize;
  (void)inputSize;

  for(Index i = 0; i < 1; ++i)
  {
    /// The array containing all the states, the measurements and the inputs
    IndexedVectorArray x;
    IndexedVectorArray u;

    /// simulation of the signal
    /// the IMU dynamical system functor
    flexibilityEstimation::IMUFixedContactDynamicalSystem imu(dt);

    /// the simulator initalization
    DynamicalSystemSimulator sim;
    sim.setDynamicsFunctor(&imu); // choice of the dynamical system

    Vector x0 = (Vector::Zero(stateSize, 1));

    /// initialization of the state vector
    for(Index j = 0; j < stateSize; ++j) x0[j] = (double(rand()) / RAND_MAX - 0.5) * 2;

    sim.setState(x0, 0);

    Vector uk = Vector::Zero(imu.getInputSize(), 1);

    for(Index k = 0; k < kmax; ++k)
    {
      u.setValue(uk, k);

      /// give the input to the simulator
      /// we only need to give one value and the
      /// simulator takes automatically the appropriate value
      sim.setInput(uk, k);
    }

    /// set the sampling perdiod to the functor
    imu.setSamplingPeriod(dt);

    sim.simulateDynamicsTo(kmax);

    x = sim.getStateArray(1, kmax);

    x.writeInFile("state.dat");
    // std::cout <<x[kmax].norm ()<< " "<< x[kmax].transpose() << std::endl;
  }

  return 0;
}

int main()
{

  return test();
}
