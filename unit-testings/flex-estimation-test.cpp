#include <fstream>
#include <iostream>

#include <state-observation/dynamical-system/dynamical-system-simulator.hpp>
#include <state-observation/examples/offline-ekf-flexibility-estimation.hpp>
#include <state-observation/noise/gaussian-white-noise.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>

using namespace stateObservation;
typedef kine::indexes<kine::rotationVector> indexes;

int testConstant()
{
  /// The number of samples
  const Index kmax = 5000;

  /// sampling period
  const double dt = 1e-3;

  /// Sizes of the states for the state, the measurement, and the input vector
  const unsigned stateSize = 18;
  const unsigned measurementSize = 6;
  const unsigned inputSize = 15;

  /// The array containing all the states, the measurements and the inputs
  IndexedVectorArray y;
  IndexedVectorArray u;

  Vector3 contact(Vector3::Zero());

  contact[0] = 0.31123600000000001;
  contact[1] = -0.35577300000000001;
  contact[2] = 1.1000000000000001;

  Vector yConstant = Vector::Zero(measurementSize, 1);
  yConstant[1] = 9.8;

  Vector uConstant = Vector::Zero(inputSize, 1);
  uConstant[2] = 1.8;

  {
    for(Index i = 1; i < kmax + 1; ++i)
    {
      y.setValue(yConstant, i);
      u.setValue(uConstant, i);
    }
  }

  u.pushBack(uConstant);

  /// the initalization of a random estimation of the initial state
  Vector xh0 = Vector::Zero(stateSize, 1);

  std::vector<Vector3, Eigen::aligned_allocator<Vector3>> contactPositions;

  contactPositions.push_back(contact);

  stateObservation::IndexedVectorArray xh =
      stateObservation::examples::offlineEKFFlexibilityEstimation(y, u, xh0, 1, contactPositions, dt);

  /// file of output
  std::ofstream f;
  f.open("trajectory.dat");

  /// the reconstruction of the state
  for(TimeIndex i = xh.getFirstIndex(); i < xh.getNextIndex(); ++i)
  {
    f << i << " " << xh[i].transpose() << std::endl;

    Vector v2(Matrix::Zero(6, 1));

    v2.head(3) = Vector(xh[i]).segment(indexes::pos, 3);
    v2.tail(3) = Vector(xh[i]).segment(indexes::ori, 3);

    Matrix4 m = kine::vector6ToHomogeneousMatrix(v2);
    Vector4 v;
    v.head(3) = contact;
    v[3] = 1;

    std::cout << m - kine::invertHomoMatrix(kine::invertHomoMatrix(m)) << std::endl << std::endl;
  }

  return 0;
}

int test()
{
  /// The number of samples
  const Index kmax = 3000;

  /// sampling period
  const double dt = 1e-3;

  /// Sizes of the states for the state, the measurement, and the input vector
  const unsigned stateSize = 18;
  const unsigned measurementSize = 6;
  // const unsigned inputSize=15;

  /// The array containing all the states, the measurements and the inputs
  IndexedVectorArray x;
  IndexedVectorArray y;
  IndexedVectorArray u;

  /// The covariance matrix of the process noise and the measurement noise
  Matrix q;
  Matrix r;

  Vector3 contact(Vector3::Random());

  {
    /// simulation of the signal
    /// the IMU dynamical system functor
    flexibilityEstimation::IMUFixedContactDynamicalSystem imu(dt);

    /// The process noise initialization
    Matrix q1 = Matrix::Zero(stateSize, stateSize);
    q1(indexes::angAcc, indexes::angAcc) = q1(indexes::angAcc + 1, indexes::angAcc + 1) =
        q1(indexes::angAcc + 2, indexes::angAcc + 2) = 0.000001;
    q1(indexes::angVel, indexes::angVel) = q1(indexes::angVel + 1, indexes::angVel + 1) =
        q1(indexes::angVel + 2, indexes::angVel + 2) = 0.0001;
    q1(indexes::ori, indexes::ori) = q1(indexes::ori + 1, indexes::ori + 1) = q1(indexes::ori + 2, indexes::ori + 2) =
        0.001;

    GaussianWhiteNoise processNoise(imu.getStateSize());
    processNoise.setStandardDeviation(q1);
    imu.setProcessNoise(&processNoise);

    /// The measurement noise initialization
    Matrix r1 = Matrix::Identity(measurementSize, measurementSize) * 0.0001;
    GaussianWhiteNoise MeasurementNoise(imu.getMeasurementSize());
    MeasurementNoise.setStandardDeviation(r1);
    imu.setMeasurementNoise(&MeasurementNoise);

    /// the simulator initalization
    DynamicalSystemSimulator sim;
    sim.setDynamicsFunctor(&imu);

    /// initialization of the state vector
    Vector x0 = Vector::Zero(stateSize, 1);

    x0[indexes::angAcc] = x0[indexes::angAcc + 1] = x0[indexes::angAcc + 2] = 0.00001;

    x0[indexes::angVel] = x0[indexes::angVel + 1] = x0[indexes::angVel + 2] = 0.0001;

    x0[indexes::ori] = x0[indexes::ori + 1] = x0[indexes::ori + 2] = 0.3;

    // x0=x0*100;

    sim.setState(x0, 0);

    Vector uk = Vector::Zero(imu.getInputSize(), 1);

    int i;
    /// construction of the input
    /// the input is constant over 10 time samples
    for(i = 0; i < kmax / 10.0; ++i)
    {
      uk[indexes::pos] = 0.4 * sin(M_PI / 10 * i);
      uk[indexes::pos + 1] = 0.6 * sin(M_PI / 12 * i);
      uk[indexes::pos + 2] = 0.2 * sin(M_PI / 5 * i);

      uk[indexes::linVel] = 0.1 * sin(M_PI / 12 * i);
      uk[indexes::linVel + 1] = 0.07 * sin(M_PI / 15 * i);
      uk[indexes::linVel + 2] = 0.05 * sin(M_PI / 5 * i);

      uk[indexes::linAcc] = 1 * sin(M_PI / 12 * i);
      uk[indexes::linAcc + 1] = 0.07 * sin(M_PI / 15 * i);
      uk[indexes::linAcc + 2] = 0.05 * sin(M_PI / 10 * i);

      uk[indexes::ori] = 2 * sin(M_PI / 12 * i);
      uk[indexes::ori + 1] = 1.5 * sin(M_PI / 18 * i);
      uk[indexes::ori + 2] = 0.8 * sin(M_PI / 6 * i);

      uk[indexes::angVel] = 0.2 * sin(M_PI / 12 * i);
      uk[indexes::angVel + 1] = 0.07 * sin(M_PI / 12 * i);
      uk[indexes::angVel + 2] = 0.05 * sin(M_PI / 5 * i);

      /// filling the 10 time samples of the constant input
      for(int j = 0; j < 10; ++j)
      {
        u.setValue(uk, i * 10 + j);
      }

      /// give the input to the simulator
      /// we only need to give one value and the
      /// simulator takes automatically the appropriate value
      sim.setInput(uk, 10 * i);
    }

    /// Last sample needed
    u.setValue(uk, i * 10);

    /// set the sampling perdiod to the functor
    imu.setSamplingPeriod(dt);

    /// launched the simulation to the time kmax+1
    for(Index i = 0; i < kmax + 1; ++i)
    {
      Vector x = sim.getState(i);

      Vector3 orientationFlexV = x.segment(indexes::ori, 3);
      Vector3 angularVelocityFlex = x.segment(indexes::angVel, 3);
      Vector3 angularAccelerationFlex = x.segment(indexes::angAcc, 3);

      Matrix3 orientationFlexR = kine::rotationVectorToAngleAxis(orientationFlexV).matrix();

      Vector3 positionFlex;
      Vector3 velocityFlex;
      Vector3 accelerationFlex;

      kine::fixedPointRotationToTranslation(orientationFlexR, angularVelocityFlex, angularAccelerationFlex, contact,
                                            positionFlex, velocityFlex, accelerationFlex);

      x.segment(indexes::pos, 3) = positionFlex;
      x.segment(indexes::linVel, 3) = velocityFlex;
      x.segment(indexes::linAcc, 3) = accelerationFlex;

      sim.setState(x, i);

      sim.simulateDynamics();
    }

    /// extract the array of measurements and states
    y = sim.getMeasurementArray(1, kmax);
    x = sim.getStateArray(1, kmax);
  }

  /// the initalization of a random estimation of the initial state
  Vector xh0 = Vector::Zero(stateSize, 1);

  std::vector<Vector3, Eigen::aligned_allocator<Vector3>> contactPositions;

  contactPositions.push_back(contact);

  stateObservation::IndexedVectorArray xh =
      stateObservation::examples::offlineEKFFlexibilityEstimation(y, u, xh0, 1, contactPositions, dt);

  /// file of output
  std::ofstream f;
  f.open("trajectory.dat");

  double error = 0;

  /// the reconstruction of the state
  for(TimeIndex i = y.getFirstIndex(); i < y.getNextIndex(); ++i)
  {
    /// display part, useless
    Vector3 g;
    {
      Matrix3 R;
      Vector3 orientationV = Vector(x[i]).segment(indexes::ori, 3);
      double angle = orientationV.norm();
      if(angle > cst::epsilonAngle)
        R = AngleAxis(angle, orientationV / angle).toRotationMatrix();
      else
        R = Matrix3::Identity();
      g = R.transpose() * Vector3::UnitZ();
      g.normalize();
    }

    Vector3 gh;
    {
      Matrix3 Rh;

      Vector3 orientationV = Vector(xh[i]).segment(indexes::ori, 3);
      double angle = orientationV.norm();
      if(angle > cst::epsilonAngle)
        Rh = AngleAxis(angle, orientationV / angle).toRotationMatrix();
      else
        Rh = Matrix3::Identity();
      gh = Rh.transpose() * Vector3::UnitZ();
      gh.normalize();
    }

    error = acos(double(g.transpose() * gh)) * 180 / M_PI;

    f << i << " \t " << error << " \t\t\t " << g.transpose() << " \t\t\t " << gh.transpose() << " \t\t\t "
      << x[i].transpose() << " \t\t\t%%%%%%\t\t\t " << xh[i].transpose() << std::endl;
  }

  std::cout << "Error " << error << ", test: ";

  if(error > 2)
  {
    std::cout << "FAILED !!!!!!!";
    return 1;
  }
  else
  {
    std::cout << "SUCCEEDED !!!!!!!";
    return 0;
  }
}

int main()
{

  return test();
}
