#include <fstream>
#include <iostream>

#include <state-observation/tools/definitions.hpp>
#include <state-observation/tools/rigid-body-kinematics.hpp>

using namespace stateObservation;

const double dt = 5e-3;

int homoMatrixDerivationTestFromFile(char * homo, char * vel, const std::string & prefix)
{
  IndexedMatrixArray velocities;
  velocities.readFromFile(vel, 6);
  IndexedMatrixArray homoMatrices;
  homoMatrices.readFromFile(homo, 4, 4);

  IndexedMatrixArray computedVelocities;

  std::ofstream f, f1, f2;
  f.open((prefix + "output.dat").c_str());
  f1.open((prefix + "computed.dat").c_str());
  f2.open((prefix + "real.dat").c_str());

  for(TimeIndex i = homoMatrices.getFirstIndex(); i < homoMatrices.getLastIndex(); ++i)
  {
    computedVelocities.setValue(kine::derivateHomogeneousMatrixFD(homoMatrices[i], homoMatrices[i + 1], dt), i);

    f << i << " \t ";
    for(long j = 0; j < computedVelocities[i].rows(); ++j)
    {
      f << " " << computedVelocities[i](j) - velocities[i](j);
    }

    f << std::endl;

    f1 << i << " \t " << computedVelocities[i].transpose() << std::endl;
    f2 << i << " \t " << velocities[i].transpose() << std::endl;
  }

  return 0;
}

int invertMatrixTest() // tested OK
{
  stateObservation::Quaternion q(stateObservation::Vector4::Random());

  q.normalize();

  Matrix3 R1(q.matrix());

  Matrix4 h(Matrix4::Identity());

  h.block(0, 0, 3, 3) = R1;
  h.block(0, 3, 3, 1) = Vector3::Random();

  Matrix4 hi = kine::invertHomoMatrix(h);

  std::cout << h << std::endl << std::endl;
  std::cout << hi << std::endl << std::endl;
  std::cout << h.inverse() << std::endl << std::endl;
  std::cout << hi * h << std::endl << std::endl;
  return 0;
}

int transformationTest() // tested ok
{

  Vector6 v2 = Vector6::Random();
  Matrix4 m = kine::vector6ToHomogeneousMatrix(v2);

  std::cout << v2.transpose() << std::endl << std::endl;
  std::cout << m << std::endl << std::endl;
  std::cout << kine::homogeneousMatrixToVector6(m) << std::endl << std::endl;
  std::cout << kine::vector6ToHomogeneousMatrix(kine::homogeneousMatrixToVector6(m)) << std::endl << std::endl;

  return 0;
}

int testHomoDerivation()
{
  Vector6 p = Vector6::Random();
  Matrix4 m = kine::vector6ToHomogeneousMatrix(p);

  Vector6 v = Vector6::Random();

  Vector6 p2;

  p2.head(3) = p.head(3) + dt * v.head(3);

  AngleAxis aa(Quaternion(kine::rotationVectorToAngleAxis(v.tail(3) * dt)) * m.block(0, 0, 3, 3));

  p2.tail(3) = aa.axis() * aa.angle();

  Matrix4 m2 = kine::vector6ToHomogeneousMatrix(p2);

  std::cout << kine::derivateHomogeneousMatrixFD(m, m2, dt).transpose() << std::endl << std::endl;

  std::cout << v.transpose() << std::endl << std::endl;
  return 0;
}

int testVector6Derivation()
{
  Vector6 p = Vector6::Random();
  Matrix4 m = kine::vector6ToHomogeneousMatrix(p);

  Vector6 v = Vector6::Random();

  Vector6 p2;

  p2.head(3) = p.head(3) + dt * v.head(3);

  AngleAxis aa(Quaternion(kine::rotationVectorToAngleAxis(v.tail(3) * dt)) * m.block(0, 0, 3, 3));

  p2.tail(3) = aa.axis() * aa.angle();

  std::cout << kine::derivatePoseThetaUFD(p, p2, dt).transpose() << std::endl << std::endl;

  std::cout << v.transpose() << std::endl << std::endl;

  return 0;
}

int main()
{
  homoMatrixDerivationTestFromFile(const_cast<char *>("/tmp/featurecompensateR_ref-position.dat"),
                                   const_cast<char *>("/tmp/featurecompensateR_ref-velocity.dat"),
                                   const_cast<char *>("reference"));

  homoMatrixDerivationTestFromFile(const_cast<char *>("/tmp/tranformation_right-gMl.dat"),
                                   const_cast<char *>("/tmp/tranformation_right-gVl.dat"),
                                   const_cast<char *>("flexInverse"));

  homoMatrixDerivationTestFromFile(const_cast<char *>("/tmp/flextimator-flexTransformationMatrix.dat"),
                                   const_cast<char *>("/tmp/flextimator-flexVelocityVector.dat"),
                                   const_cast<char *>("flex"));

  // invertMatrixTest();

  // transformationTest();

  testHomoDerivation();

  testVector6Derivation();

  return 0;
}
