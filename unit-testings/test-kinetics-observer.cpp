#include <iostream>

#include <state-observation/dynamics-estimators/kinetics-observer.hpp>
#include <state-observation/tools/probability-law-simulation.hpp>

using namespace stateObservation;
using namespace kine;

int testKineticsObserverCodeAccessor(int errorcode)
{
  double error = 0;
  double dt = 0.001;

  double mass = 100;
  KineticsObserver o(4, 2);

  Vector x0(o.getStateSize());
  x0.setZero();
  Vector xf(x0);
  Vector xs(x0);

  o.setSamplingTime(dt);

  Kinematics stateKine;

  stateKine.position.set() << 0.1, 0, 0.7;
  stateKine.orientation = Vector3(0, 0, 0);
  stateKine.linVel.set().setZero();
  stateKine.angVel.set().setZero();

  o.setStateVector(x0);

  o.setStateKinematics(stateKine);
  o.setGyroBias(Vector3(1, 2, 3));

  Vector6 wrench;
  wrench << 4, 5, 6, 7, 8, 9;

  o.setStateUnmodeledWrench(wrench);

  Vector x = o.getCurrentStateVector();
  stateObservation::TimeIndex index = o.getStateVectorTimeIndex();

  Kinematics contactKine;
  contactKine.position.set() << 0, 0.1, 0;
  contactKine.orientation.setZeroRotation();

  o.addContact(contactKine, 0);

  Matrix3 linStiffness, angStiffness, linDamping, angDamping;
  linStiffness.setZero();
  linStiffness.diagonal().setConstant(50000);

  angStiffness.setZero();
  angStiffness.diagonal().setConstant(400);

  linDamping.setZero();
  linDamping.diagonal().setConstant(500);

  angDamping.setZero();
  angDamping.diagonal().setConstant(20);

  contactKine.position.set() << 0, -0.1, 0;
  o.addContact(contactKine, 3, linStiffness, linDamping, angStiffness, angDamping);

  Matrix12 initialCov, processCov;

  initialCov.setZero();
  initialCov.diagonal().setConstant(0.01);
  processCov.setZero();
  processCov.diagonal().setConstant(0.0001);

  contactKine.position.set() << 1, 0.1, 0;
  int i = o.addContact(contactKine, initialCov, processCov);

  (void)i; /// avoid warning in release mode
  assert(i == 1);

  contactKine.position.set() << 1, -0.1, 0;
  o.addContact(contactKine, initialCov, processCov, 2, linDamping, linStiffness, angStiffness, angDamping);

  std::cout << index << " " << x.transpose() << std::endl;

  o.update();

  Kinematics k = o.getKinematics();

  std::cout << k;

  Kinematics l = o.getKinematicsOf(k);

  std::cout << l;

  std::cout << o.kineIndex() << " " << o.posIndex() << " " << o.oriIndex() << " " << o.linVelIndex() << " "
            << o.angVelIndex() << " " << o.gyroBiasIndex(0) << " " << o.gyroBiasIndex(1) << " "
            << o.unmodeledWrenchIndex() << " " << o.unmodeledForceIndex() << " " << o.unmodeledTorqueIndex() << " "
            << o.contactsIndex() << " " << o.contactIndex(0) << " " << o.contactKineIndex(0) << " "
            << o.contactPosIndex(0) << " " << o.contactOriIndex(0) << " " << o.contactForceIndex(0) << " "
            << o.contactTorqueIndex(0) << " " << o.contactWrenchIndex(0) << " " <<

      o.contactIndex(1) << " " << o.contactKineIndex(1) << " " << o.contactPosIndex(1) << " " << o.contactOriIndex(1)
            << " " << o.contactForceIndex(1) << " " << o.contactTorqueIndex(1) << " " << o.contactWrenchIndex(1) << " "
            <<

      o.contactIndex(2) << " " << o.contactKineIndex(2) << " " << o.contactPosIndex(2) << " " << o.contactOriIndex(2)
            << " " << o.contactForceIndex(2) << " " << o.contactTorqueIndex(2) << " " << o.contactWrenchIndex(2) << " "
            <<

      o.contactIndex(3) << " " << o.contactKineIndex(3) << " " << o.contactPosIndex(3) << " " << o.contactOriIndex(3)
            << " " << o.contactForceIndex(3) << " " << o.contactTorqueIndex(3) << " " << o.contactWrenchIndex(3) << " "
            << std::endl;

  std::cout << o.kineIndexTangent() << " " << o.posIndexTangent() << " " << o.oriIndexTangent() << " "
            << o.linVelIndexTangent() << " " << o.angVelIndexTangent() << " " << o.gyroBiasIndexTangent(0) << " "
            << o.gyroBiasIndexTangent(1) << " " << o.unmodeledWrenchIndexTangent() << " "
            << o.unmodeledForceIndexTangent() << " " << o.unmodeledTorqueIndexTangent() << " "
            << o.contactsIndexTangent() << " " <<

      o.contactIndexTangent(0) << " " << o.contactKineIndexTangent(0) << " " << o.contactPosIndexTangent(0) << " "
            << o.contactOriIndexTangent(0) << " " << o.contactForceIndexTangent(0) << " "
            << o.contactTorqueIndexTangent(0) << " " << o.contactWrenchIndexTangent(0) << " " <<

      o.contactIndexTangent(1) << " " << o.contactKineIndexTangent(1) << " " << o.contactPosIndexTangent(1) << " "
            << o.contactOriIndexTangent(1) << " " << o.contactForceIndexTangent(1) << " "
            << o.contactTorqueIndexTangent(1) << " " << o.contactWrenchIndexTangent(1) << " " <<

      o.contactIndexTangent(2) << " " << o.contactKineIndexTangent(2) << " " << o.contactPosIndexTangent(2) << " "
            << o.contactOriIndexTangent(2) << " " << o.contactForceIndexTangent(2) << " "
            << o.contactTorqueIndexTangent(2) << " " << o.contactWrenchIndexTangent(2) << " " <<

      o.contactIndexTangent(3) << " " << o.contactKineIndexTangent(3) << " " << o.contactPosIndexTangent(3) << " "
            << o.contactOriIndexTangent(3) << " " << o.contactForceIndexTangent(3) << " "
            << o.contactTorqueIndexTangent(3) << " " << o.contactWrenchIndexTangent(3) << " " << std::endl;

  o.setWithUnmodeledWrench(true);
  o.setWithAccelerationEstimation(true);
  o.setWithGyroBias(true);

  Matrix3 acceleroCov, gyroCov;

  acceleroCov = Matrix3::Identity() * 1e-4;
  gyroCov = Matrix3::Identity() * 1e-8;

  o.setIMUDefaultCovarianceMatrix(acceleroCov, gyroCov);

  Matrix6 wrenchCov;

  wrenchCov << Matrix3::Identity() * 1e-0, Matrix3::Zero(), Matrix3::Zero(), Matrix3::Identity() * 1e-4;

  o.setContactWrenchSensorDefaultCovarianceMatrix(wrenchCov);

  o.setMass(mass);

  Vector state1 = o.getEKF().stateVectorRandom();
  Vector state2 = o.getEKF().stateVectorRandom();
  Vector statediff;
  o.stateDifference(state1, state2, statediff);
  Vector state3;
  o.stateSum(state2, statediff, state3);

  std::cout << state1.transpose() << std::endl;
  std::cout << state3.transpose() << std::endl;

  Matrix statecomp(state1.size(), 2);

  statecomp << state1, state3;

  std::cout << statecomp << std::endl;

  std::cout << "Sum error" << (error = o.stateDifference(state1, state3).norm()) << std::endl;

  state2 = o.getEKF().stateVectorRandom();
  statediff = o.getEKF().stateTangentVectorRandom();

  Vector statediff_bis;
  o.stateSum(state2, statediff, state1);
  o.stateDifference(state1, state2, statediff_bis);

  std::cout << statediff.transpose() << std::endl;
  std::cout << statediff_bis.transpose() << std::endl;

  Matrix statecompdiff(statediff.size(), 2);

  statecompdiff << statediff, statediff_bis;

  std::cout << statecompdiff << std::endl;

  std::cout << "DIff error" << (error += (statediff - statediff_bis).norm()) << std::endl;

  if(error > 1e-8)
  {
    return errorcode;
  }

  o.clearContacts();

  return 0;
}

int main()
{
  int returnVal;

  if((returnVal = testKineticsObserverCodeAccessor(3)))
  {
    std::cout << "Kinetics Observer test failed, code : 3" << std::endl;
    return returnVal;
  }
  else
  {
    std::cout << "Kinetics Observer test succeeded" << std::endl;
  }

  std::cout << "test succeeded" << std::endl;
  return 0;
}
