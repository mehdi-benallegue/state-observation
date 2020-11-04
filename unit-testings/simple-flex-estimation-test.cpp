#include <fstream>
#include <iostream>

#include <state-observation/dynamical-system/dynamical-system-simulator.hpp>
#include <state-observation/examples/offline-ekf-flexibility-estimation.hpp>
#include <state-observation/noise/gaussian-white-noise.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>

using namespace stateObservation;

typedef kine::indexes<kine::rotationVector> indexes;

int test()
{
  /// The number of samples
  const Index kmax = 1000;

  /// sampling period
  const double dt = 5e-3;

  /// Sizes of the states for the state, the measurement, and the input vector
  const unsigned stateSize = 18;
  // const unsigned measurementSize=6;
  const unsigned inputSize = 15;

  /// The array containing all the states, the measurements and the inputs
  IndexedVectorArray x;
  IndexedVectorArray y;
  IndexedVectorArray u;
  IndexedVectorArray z;

  IndexedVectorArray ino;
  IndexedVectorArray prediMea;

  /// Contact vector
  Vector3 contact(-1, 0, 0);

  /// Generation
  {
    Quaternion q(Quaternion::Identity());
    Vector3 odoti(Vector3::Zero());
    Vector3 oi(Vector3::Zero());
    Vector3 pos;
    Vector3 vel;
    Vector3 acc;
    kine::fixedPointRotationToTranslation(q.matrix(), oi, odoti, contact, pos, vel, acc);
    AngleAxis aa(q);

    Vector Xi(Vector::Zero(stateSize, 1));
    Xi.segment(indexes::pos, 3) = pos;
    Xi.segment(indexes::ori, 3) = aa.angle() * aa.axis();
    Xi.segment(indexes::linVel, 3) = vel;
    Xi.segment(indexes::angVel, 3) = oi;
    Xi.segment(indexes::linAcc, 3) = acc;
    Xi.segment(indexes::angAcc, 3) = odoti;
    x.setValue(Xi, 0);

    Quaternion qCtrl(Quaternion::Identity());
    Vector3 odotCtrl(Vector3::Zero());
    Vector3 oCtrl(Vector3::Zero());
    Vector3 posCtrl(contact);
    Vector3 velCtrl(Vector3::Zero());
    Vector3 accCtrl(Vector3::Zero());

    Quaternion qImu(q * qCtrl);
    Vector3 oImu(Vector3::Zero());
    Vector3 posImu(q.matrix() * posCtrl + pos);
    Vector3 velImu(Vector3::Zero());
    Vector3 accImu(Vector3::Zero());

    AccelerometerGyrometer imu;

    for(Index i = 1; i < kmax; ++i)
    {
      q = kine::rotationVectorToAngleAxis(oi * dt) * q;
      aa = q;
      oi += odoti * dt;
      double id = double(i);
      odoti << 0.1 * sin(0.007 * id), 0.2 * sin(0.03 * id), 0.25 * sin(0.02 * id);

      kine::fixedPointRotationToTranslation(q.matrix(), oi, odoti, contact, pos, vel, acc);

      Xi.segment(indexes::pos, 3) = pos;
      Xi.segment(indexes::ori, 3) = aa.angle() * aa.axis();
      Xi.segment(indexes::linVel, 3) = vel;
      Xi.segment(indexes::angVel, 3) = oi;
      Xi.segment(indexes::linAcc, 3) = acc;
      Xi.segment(indexes::angAcc, 3) = odoti;

      x.setValue(Xi, i);

      qCtrl = kine::rotationVectorToAngleAxis(oCtrl * dt) * qCtrl;
      AngleAxis aaCtrl(qCtrl);
      oCtrl += odotCtrl * dt;

      odotCtrl << 0.15 * sin(0.008 * id), 0.1 * sin(0.023 * id), 0.2 * sin(0.025 * id);
      posCtrl += velCtrl * dt;
      velCtrl += accCtrl * dt;
      accCtrl << 0.12 * sin(0.018 * id), 0.08 * sin(0.035 * id), 0.3 * sin(0.027 * id);

      Vector Ui(Vector::Zero(inputSize, 1));
      Ui.segment(indexes::pos, 3) = posCtrl;
      Ui.segment(indexes::ori, 3) = aaCtrl.angle() * aaCtrl.axis();
      Ui.segment(indexes::linVel, 3) = velCtrl;
      Ui.segment(indexes::angVel, 3) = oCtrl;
      Ui.segment(indexes::linAcc, 3) = accCtrl;
      u.setValue(Ui, i);

      Quaternion newqImu(q * qCtrl);
      Vector3 newPosImu(q.matrix() * posCtrl + pos);
      Vector3 newVelImu(tools::derivate(posImu, newPosImu, dt));

      accImu = tools::derivate(velImu, newVelImu, dt);
      velImu = newVelImu;
      posImu = newPosImu;

      oImu = kine::derivateRotationFD(qImu, newqImu, dt);
      qImu = newqImu;

      Vector Ximu(Vector::Zero(10, 1));
      Ximu.head<4>() = qImu.coeffs();
      Ximu.segment(4, 3) = accImu;
      Ximu.segment(7, 3) = oImu;

      imu.setState(Ximu, i);
      z.setValue(Ximu, i);
      y.setValue(imu.getMeasurements(), i);
    }
  }

  /// the initalization of an estimation of the initial state
  Vector xh0 = Vector::Zero(stateSize, 1);

  std::vector<Vector3, Eigen::aligned_allocator<Vector3>> contactPositions;

  contactPositions.push_back(contact);

  stateObservation::IndexedVectorArray xh =
      stateObservation::examples::offlineEKFFlexibilityEstimation(y, u, xh0, 1, contactPositions, dt, &ino, &prediMea);

  double error = 0;

  /// the reconstruction of the state
  for(TimeIndex i = y.getFirstIndex(); i < y.getNextIndex(); ++i)
  {
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
  }

  std::cout << "Error " << error << ", test: ";

  if(error > 2.)
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
