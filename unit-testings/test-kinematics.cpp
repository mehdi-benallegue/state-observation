#include <iostream>

#include <state-observation/tools/probability-law-simulation.hpp>
#include <state-observation/tools/rigid-body-kinematics.hpp>

using namespace stateObservation;
using namespace kine;

/// @brief test rotationMatrix2Angle
///
/// @param errorCode
/// @return int
int testRotationOperations(int errorCode)
{
  unsigned numberOfTests = 1000;
  for(unsigned currentTest = 0; currentTest < numberOfTests; ++currentTest)
  {
    { /// test isRotation

      Matrix3 R = Matrix3::Random();
      //(M.isUnitary() && isApprox(M.topLeftCorner<2, 2>().determinant(), M(2, 2))
      //&& isApprox(M.bottomLeftCorner<2, 2>().determinant(), M(0, 2)));
      if(isRotationMatrix(R, cst::epsilon1 * 10))
      {
        std::cout << "Test number " << currentTest << "isRotationTest failed: false positive" << std::endl;
        return errorCode;
      }
      R = randomRotationQuaternion().toRotationMatrix();
      if(!isRotationMatrix(R, cst::epsilon1 * 10))
      {
        std::cout << R.isUnitary() << " " << R.topLeftCorner<2, 2>().determinant() << " " << R(2, 2) << " "
                  << R.bottomLeftCorner<2, 2>().determinant() << " " << R(0, 2) << " "
                  << R.topLeftCorner<2, 2>().determinant() - R(2, 2) << " "
                  << R.bottomLeftCorner<2, 2>().determinant() - R(0, 2) << std::endl;
        std::cout << "Test number " << currentTest << "isRotationTest failed: false negative" << std::endl;
        return errorCode;
      }
      Matrix3 R2;
      R2 << R.col(1), R.col(0), R.col(2); /// Non right handed orthogonal matrix
      if(isRotationMatrix(R2, cst::epsilon1 * 10))
      {
        std::cout << "isRotationTest failed: false positive (right-handedness)" << std::endl;
        return errorCode;
      }
    }
    {
      /// test pure yaw
      if(!isPureYaw(AngleAxis(randomAngle(), Vector3::UnitZ()).matrix())
         || isPureYaw(randomRotationQuaternion().toRotationMatrix()))
      {
        std::cout << "Test number " << currentTest << "Pure yaw detection failed " << std::endl;
        return errorCode;
      }
    }
    {
      /// Test the angle made by one vector by a rotation
      Vector3 axis = Vector3::Random().normalized();
      Vector3 v = Vector3::Random().cross(axis).normalized();
      double angle = randomAngle();
      Matrix3 m = (AngleAxis(angle, axis)).matrix() * AngleAxis(randomAngle(), v).matrix();

      double error = fabs(angle - kine::rotationMatrixToAngle(m, axis, v));

      if(error > cst::epsilon1 * 2 * M_PI)
      {
        std::cout << "Test number " << currentTest << "Test vector angle failed. Angle error " << error << std::endl;
        return errorCode;
      }
    }
    {
      /// Test yaw extraction with custom vector
      double angle = randomAngle();

      Vector2 v = Vector2::Random().normalized();
      Vector3 v3;
      v3 << v, 0;
      Matrix3 m = AngleAxis(angle, Vector3::UnitZ()).matrix() * AngleAxis(randomAngle(), v3).matrix();

      double error = fabs(angle - kine::rotationMatrixToYaw(m, v));

      if(error > cst::epsilon1 * 2 * M_PI)
      {
        std::cout << "Test number " << currentTest << "Test Yaw extraction with custom vector failed. Angle error "
                  << error << std::endl;
        return errorCode;
      }
    }
    {
      /// Test yaw extraction from roll pitch yaw using the x axis alignment
      double rollangle = randomAngle();
      double pitchangle = randomAngle() / 2; /// constrain the pitch to be smaller than pi/2
      double yawangle = randomAngle();

      Matrix3 m = rollPitchYawToRotationMatrix(rollangle, pitchangle, yawangle);

      double error = fabs(yawangle - kine::rotationMatrixToYaw(m));

      if(error > cst::epsilon1 * 1e5 * M_PI) /// this function is really not precise
      {
        std::cout << "Test number " << currentTest << " axis-based failed. Angle" << yawangle << " Angle error "
                  << error << std::endl;
        return errorCode;
      }

      /// Test yaw extraction from roll pitch yaw using the traditional conversion
      Vector3 rpy = kine::rotationMatrixToRollPitchYaw(m);
      error = fabs(yawangle - rpy(2));

      if(error > cst::epsilon1 * 1e5 * M_PI)
      {
        std::cout << "Test number " << currentTest << "eigen-based failed. Angle" << yawangle << " Angle error "
                  << error << std::endl;
        return errorCode;
      }
    }

    {
      /// Test the automatic detection of the vector to use to extract yaw from a rotation matrix
      double angle = randomAngle(); /// random value
      Vector2 v = Vector2::Random().normalized();
      Vector3 v3;
      v3 << v, 0;

      Matrix3 m = AngleAxis(angle, Vector3::UnitZ()).matrix() * AngleAxis(randomAngle(), v3).matrix();

      double error = fabs(angle - kine::rotationMatrixToYawAxisAgnostic(m));

      if(error > cst::epsilon1 * 2 * M_PI)
      {
        std::cout << " Test rotationMatrixToYawAxisAgnostic failed. Angle error " << error << std::endl;

        return errorCode;
      }
    }
    {
      Vector3 horizAxis1;
      horizAxis1(2) = 0;
      horizAxis1.head<2>() = Vector2::Random().normalized();
      double tiltAngle1 = randomAngle();
      double yawAngle = randomAngle();

      Matrix3 yaw = AngleAxis(yawAngle, Vector3::UnitZ()).matrix();

      /// we should get back this matrix
      Matrix3 initialMatrix = yaw * AngleAxis(tiltAngle1, horizAxis1).matrix();

      Vector3 tilt = initialMatrix.transpose() * Vector3::UnitZ();

      Vector3 horizAxis2;
      horizAxis2(2) = 0;
      horizAxis2.head<2>() = Vector2::Random().normalized();

      Vector3 ml = initialMatrix.transpose() * horizAxis2;

      Vector3 newTilt = Vector3::Random();
      Matrix3 mTemp1, mTemp2;

      mTemp1 << horizAxis2, horizAxis2.cross(Vector3::UnitZ().cross(horizAxis2)).normalized(),
          horizAxis2.cross(Vector3::UnitZ()).normalized();

      mTemp2 << ml, ml.cross(newTilt.cross(ml)).normalized(), ml.cross(newTilt).normalized();

      Matrix3 M2 = mTemp1 * mTemp2.transpose();

      Matrix3 estimatedMatrix = kine::mergeTiltWithYawAxisAgnostic(tilt, M2);
      std::cout << " " << isRotationMatrix(M2, 0.001) << " " << isRotationMatrix(mTemp1, 0.0001) << " "
                << isRotationMatrix(mTemp2, 0.0001) << std::endl;

      if(!isRotationMatrix(estimatedMatrix, cst::epsilon1 * 10))
      {
        std::cout << "Test mergeTiltWithYawAxisAgnostic failed. Reconstructed matrix is not a Rotation Matrix"
                  << std::endl;
        return errorCode;
      }

      double error = AngleAxis(initialMatrix * estimatedMatrix.transpose()).angle();

      if(error > cst::epsilon1 * 1000 * M_PI)
      {
        std::cout << "Test mergeTiltWithYawAxisAgnostic failed. Reconstructed matrix is wrong" << std::endl;
        return errorCode;
      }
    }
  } /// end of for
  return 0;
}

int testOrientation(int errcode)
{

  Vector4 q_i;

  typedef tools::ProbabilityLawSimulation ran;

  /// random orientation
  q_i = ran::getGaussianVector(Matrix4::Identity(), Vector4::Zero(), 4);
  q_i.normalize();

  /// several representations of the orientation
  Quaternion q(q_i);
  AngleAxis aa(q);
  Matrix3 M(q.toRotationMatrix());

  double err = 0.;

  std::cout << "Orientation test started " << err << std::endl;

  {
    kine::Orientation ori1;
    ori1 = q;
    Matrix3 M1 = ori1;

    err +=
        (Quaternion(ori1.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm() + (M1 - M).norm();
    std::cout << "Assignment operaton test 1 (Quaternion) done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori2;
    ori2 = M;

    err += (Quaternion(ori2.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm()
           + (ori2.toMatrix3() - M).norm();
    std::cout << "Assignment operaton test 2 (Matrix3) done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori3;
    ori3 = aa;

    err += (Quaternion(ori3.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm()
           + (ori3.toMatrix3() - M).norm();
    std::cout << "Assignment operaton test 3 (AngleAxis) done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori4;
    ori4 = Vector3(aa.angle() * aa.axis());

    err += (Quaternion(ori4.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm()
           + (ori4.toMatrix3() - M).norm();
    std::cout << "Assignment operaton test 3 (Vector3) done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori1(q);
    Matrix3 M1 = ori1;

    err +=
        (Quaternion(ori1.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm() + (M1 - M).norm();
    std::cout << "Cast operaton test 1 (Matrix3) done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori2(M);

    err += (Quaternion(ori2.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm()
           + (ori2.toMatrix3() - M).norm();
    std::cout << "copy constructor test 1 (Matrix3) done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori3(aa);

    err += (Quaternion(ori3.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm()
           + (ori3.toMatrix3() - M).norm();
    std::cout << "copy constructor test 2 (AngleAxis) done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori4(Vector3(aa.angle() * aa.axis()));

    err += (Quaternion(ori4.toVector4()).toRotationMatrix() - M).norm() + (ori4.toMatrix3() - M).norm();
    std::cout << "copy constructor test 3 (AngleAxis) done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori4(q, M);

    err += (Quaternion(ori4.toVector4()).toRotationMatrix() - M).norm() + (ori4.toMatrix3() - M).norm();
    std::cout << "Constructor test 1 (Quaternion, Matrix3) done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori4;
    ori4.setRandom();
    ori4 = M;

    err += (Quaternion(ori4.toVector4()).toRotationMatrix() - M).norm() + (ori4.toMatrix3() - M).norm();
    std::cout << "Update Assignement operator test 1 (Matrix3) done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori4;
    ori4.setRandom();
    ori4 = q;

    err += (Quaternion(ori4.toVector4()).toRotationMatrix() - M).norm() + (ori4.toMatrix3() - M).norm();
    std::cout << "Update Assignement operator test 2 (Quaternion) done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori1(q);

    kine::Orientation ori2 = ori1.inverse();
    kine::Orientation ori3 = ori2.inverse();

    Matrix3 M3 = ori1;

    err += (Quaternion(ori3.toVector4()).toRotationMatrix() - M).norm() + (M3 - M).norm();
    std::cout << "Inverse operator test 1 done. Current error " << err << std::endl;
  }

  {

    kine::Orientation ori0(q);
    kine::Orientation ori1(M);
    ori1 = ori1.inverse();
    const kine::Orientation & ori2 = ori0;
    const kine::Orientation & ori3 = ori1;

    kine::Orientation ori00 = ori0 * ori1;
    ori0 = q;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori01 = ori1 * ori0;
    ori0 = q;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori02 = ori0 * ori3;
    ori0 = q;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori03 = ori3 * ori0;
    ori0 = q;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori04 = ori2 * ori1;
    ori0 = q;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori05 = ori1 * ori2;
    ori0 = q;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori06 = ori2 * ori3;
    kine::Orientation ori07 = ori3 * ori2;

    kine::Orientation ori10(ori0, ori1);
    ori0 = q;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori11(ori1, ori0);
    ori0 = q;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori12(ori0, ori3);
    ori0 = q;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori13(ori3, ori0);
    ori0 = q;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori14(ori2, ori1);
    ori0 = q;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori15(ori1, ori2);
    ori0 = q;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori16(ori2, ori3);
    kine::Orientation ori17(ori3, ori2);

    err += (ori00.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  1 done. Current error " << err << std::endl;
    err += (ori01.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  2 done. Current error " << err << std::endl;
    err += (ori02.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  3 done. Current error " << err << std::endl;
    err += (ori03.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  4 done. Current error " << err << std::endl;
    err += (ori04.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  5 done. Current error " << err << std::endl;
    err += (ori05.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  6 done. Current error " << err << std::endl;
    err += (ori06.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  7 done. Current error " << err << std::endl;
    err += (ori07.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  8 done. Current error " << err << std::endl;
    err += (ori10.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  9 done. Current error " << err << std::endl;
    err += (ori11.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  10 done. Current error " << err << std::endl;
    err += (ori12.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  11 done. Current error " << err << std::endl;
    err += (ori13.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  12 done. Current error " << err << std::endl;
    err += (ori14.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  13 done. Current error " << err << std::endl;
    err += (ori15.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  14 done. Current error " << err << std::endl;
    err += (ori16.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  15 done. Current error " << err << std::endl;
    err += (ori17.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  16 done. Current error " << err << std::endl;
  }

  {

    kine::Orientation ori0(q);
    kine::Orientation ori1(q);
    ori1 = ori1.inverse();
    const kine::Orientation & ori2 = ori0;
    const kine::Orientation & ori3 = ori1;

    kine::Orientation ori00 = ori0 * ori1;
    ori0 = q;
    ori1 = q;
    ori1 = ori1.inverse();
    kine::Orientation ori01 = ori1 * ori0;
    ori0 = q;
    ori1 = q;
    ori1 = ori1.inverse();
    kine::Orientation ori02 = ori0 * ori3;
    ori0 = q;
    ori1 = q;
    ori1 = ori1.inverse();
    kine::Orientation ori03 = ori3 * ori0;
    ori0 = q;
    ori1 = q;
    ori1 = ori1.inverse();
    kine::Orientation ori04 = ori2 * ori1;
    ori0 = q;
    ori1 = q;
    ori1 = ori1.inverse();
    kine::Orientation ori05 = ori1 * ori2;
    ori0 = q;
    ori1 = q;
    ori1 = ori1.inverse();
    kine::Orientation ori06 = ori2 * ori3;
    kine::Orientation ori07 = ori3 * ori2;

    kine::Orientation ori10(ori0, ori1);
    ori0 = q;
    ori1 = q;
    ori1 = ori1.inverse();
    kine::Orientation ori11(ori1, ori0);
    ori0 = q;
    ori1 = q;
    ori1 = ori1.inverse();
    kine::Orientation ori12(ori0, ori3);
    ori0 = q;
    ori1 = q;
    ori1 = ori1.inverse();
    kine::Orientation ori13(ori3, ori0);
    ori0 = q;
    ori1 = q;
    ori1 = ori1.inverse();
    kine::Orientation ori14(ori2, ori1);
    ori0 = q;
    ori1 = q;
    ori1 = ori1.inverse();
    kine::Orientation ori15(ori1, ori2);
    ori0 = q;
    ori1 = q;
    ori1 = ori1.inverse();
    kine::Orientation ori16(ori2, ori3);
    kine::Orientation ori17(ori3, ori2);

    err += (ori00.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  17 done. Current error " << err << std::endl;
    err += (ori01.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  18 done. Current error " << err << std::endl;
    err += (ori02.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  19 done. Current error " << err << std::endl;
    err += (ori03.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  20 done. Current error " << err << std::endl;
    err += (ori04.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  21 done. Current error " << err << std::endl;
    err += (ori05.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  22 done. Current error " << err << std::endl;
    err += (ori06.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  23 done. Current error " << err << std::endl;
    err += (ori07.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  24 done. Current error " << err << std::endl;
    err += (ori10.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  25 done. Current error " << err << std::endl;
    err += (ori11.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  26 done. Current error " << err << std::endl;
    err += (ori12.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  27 done. Current error " << err << std::endl;
    err += (ori13.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  28 done. Current error " << err << std::endl;
    err += (ori14.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  29 done. Current error " << err << std::endl;
    err += (ori15.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  30 done. Current error " << err << std::endl;
    err += (ori16.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  31 done. Current error " << err << std::endl;
    err += (ori17.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  32 done. Current error " << err << std::endl;
  }

  {

    kine::Orientation ori0(M);
    kine::Orientation ori1(M);
    ori1 = ori1.inverse();
    const kine::Orientation & ori2 = ori0;
    const kine::Orientation & ori3 = ori1;

    ori0 = M;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori00 = ori0 * ori1;
    ori0 = M;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori01 = ori1 * ori0;
    ori0 = M;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori02 = ori0 * ori3;
    ori0 = M;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori03 = ori3 * ori0;
    ori0 = M;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori04 = ori2 * ori1;
    ori0 = M;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori05 = ori1 * ori2;
    ori0 = M;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori06 = ori2 * ori3;
    kine::Orientation ori07 = ori3 * ori2;

    kine::Orientation ori10(ori0, ori1);
    ori0 = M;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori11(ori1, ori0);
    ori0 = M;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori12(ori0, ori3);
    ori0 = M;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori13(ori3, ori0);
    ori0 = M;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori14(ori2, ori1);
    ori0 = M;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori15(ori1, ori2);
    ori0 = M;
    ori1 = M;
    ori1 = ori1.inverse();
    kine::Orientation ori16(ori2, ori3);
    kine::Orientation ori17(ori3, ori2);

    err += (ori00.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  33 done. Current error " << err << std::endl;
    err += (ori01.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  34 done. Current error " << err << std::endl;
    err += (ori02.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  35 done. Current error " << err << std::endl;
    err += (ori03.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  36 done. Current error " << err << std::endl;
    err += (ori04.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  37 done. Current error " << err << std::endl;
    err += (ori05.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  38 done. Current error " << err << std::endl;
    err += (ori06.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  39 done. Current error " << err << std::endl;
    err += (ori07.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  40 done. Current error " << err << std::endl;
    err += (ori10.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  41 done. Current error " << err << std::endl;
    err += (ori11.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  42 done. Current error " << err << std::endl;
    err += (ori12.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  43 done. Current error " << err << std::endl;
    err += (ori13.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  44 done. Current error " << err << std::endl;
    err += (ori14.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  45 done. Current error " << err << std::endl;
    err += (ori15.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  46 done. Current error " << err << std::endl;
    err += (ori16.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  47 done. Current error " << err << std::endl;
    err += (ori17.toMatrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  48 done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori1(q);
    kine::Orientation ori2(M);

    kine::Orientation ori3 = ori2.inverse() * ori1;

    err += (Quaternion(ori3.toVector4()).toRotationMatrix() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  49 done. Current error " << err << std::endl;
  }

  {
    Vector3 v1 = Vector3::Random() * 0.1;
    kine::Orientation ori1(M);
    kine::Orientation ori2(ori1);
    ori2.integrate(v1);

    Vector3 v2 = ori1.differentiate(ori2);

    err += (v1 - v2).norm();
    std::cout << "Integration/differentiation test 1 done. Current error " << err << std::endl;
  }

  {
    Vector3 v1;
    kine::Orientation ori1;
    kine::Orientation ori2;

    ori1.setRandom();
    ori2.setRandom();

    v1 = ori2.differentiate(ori1);
    kine::Orientation oris = ori2.integrate(v1);

    err += oris.differentiate(ori1).norm();
    std::cout << "Integration/differentiation test 2 done. Current error " << err << std::endl;
  }

  if(err > 1e-13)
  {
    return errcode;
  }
  return 0;
}

int testKinematics(int errcode)
{

  std::cout << "Kinematics test started" << std::endl;
  typedef kine::Kinematics::Flags Flags;

  kine::Kinematics k0;
  kine::Kinematics k1;

  Flags::Byte flag0 = BOOST_BINARY(000000);
  Flags::Byte flag1 = BOOST_BINARY(000000);

  kine::Kinematics k, l, k2;

  int count = int(pow(2, 6) * pow(2, 6));
  double err = 0;
  double threshold = 1e-30 * count;

  for(int i = 0; i < count; i++)
  {
    Vector3 pos0 = Vector3::Random();
    kine::Orientation ori0 = kine::Orientation::randomRotation();
    Vector3 linvel0 = Vector3::Random();
    Vector3 angvel0 = Vector3::Random();
    Vector3 linacc0 = Vector3::Random();
    Vector3 angacc0 = Vector3::Random();

    Vector3 pos1 = Vector3::Random();
    kine::Orientation ori1 = kine::Orientation::randomRotation();
    Vector3 linvel1 = Vector3::Random();
    Vector3 angvel1 = Vector3::Random();
    Vector3 linacc1 = Vector3::Random();
    Vector3 angacc1 = Vector3::Random();

    k0.reset();
    k1.reset();

    if(flag0 & Flags::position)
    {
      k0.position = pos0;
    }
    if(true) /// the orientation has to be set
    {
      k0.orientation = ori0;
    }
    if(flag0 & Flags::linVel)
    {
      k0.linVel = linvel0;
    }
    if(flag0 & Flags::angVel)
    {
      k0.angVel = angvel0;
    }
    if(flag0 & Flags::linAcc)
    {
      k0.linAcc = linacc0;
    }
    if(flag0 & Flags::angAcc)
    {
      k0.angAcc = angacc0;
    }

    if(flag1 & Flags::position)
    {
      k1.position = pos1;
    }
    if(!k1.position.isSet() || flag1 & Flags::orientation) /// the orientation has to be set
    {
      k1.orientation = ori1;
    }
    if(flag1 & Flags::linVel)
    {
      k1.linVel = linvel1;
    }
    if(flag1 & Flags::angVel)
    {
      k1.angVel = angvel1;
    }
    if(flag1 & Flags::linAcc)
    {
      k1.linAcc = linacc1;
    }
    if(flag1 & Flags::angAcc)
    {
      k1.angAcc = angacc1;
    }

    if((flag0 = (flag0 + 1) & Flags::all) == 0) /// update the flags to span all the possibilties
      flag1 = (flag1 + 1) & Flags::all;

    k2 = k1;
    if(!k1.orientation.isSet())
    {
      k2.orientation.setZeroRotation();
    }

    k = ((k0 * k2) * k2.getInverse()) * k0.getInverse();

    if(k.position.isSet())
    {
      err += k.position().squaredNorm();
    }
    if(k.orientation.isSet())
    {
      err += k.orientation.toRotationVector().squaredNorm();
    }
    if(k.linVel.isSet())
    {
      err += k.linVel().squaredNorm();
    }
    if(k.angVel.isSet())
    {
      err += k.angVel().squaredNorm();
    }
    if(k.linAcc.isSet())
    {
      err += k.linAcc().squaredNorm();
    }
    if(k.angAcc.isSet())
    {
      err += k.angAcc().squaredNorm();
    }
  }

  std::cout << "Error 1 : " << err << std::endl;

  if(err > threshold)
  {
    std::cout << "Error too large : " << err << std::endl;
    return errcode;
  }

  err = 0;

  threshold = 1e-19 * count;

  flag0 = BOOST_BINARY(000000);
  flag1 = BOOST_BINARY(000000);

  for(int i = 0; i < count; i++)
  {
    Vector3 pos0 = Vector3::Random();
    kine::Orientation ori0 = kine::Orientation::randomRotation();
    Vector3 linvel0 = Vector3::Random();
    Vector3 angvel0 = Vector3::Random();
    Vector3 linacc0 = Vector3::Random();
    Vector3 angacc0 = Vector3::Random();

    k.reset();
    k.position = pos0;
    k.orientation = ori0;
    k.linVel = linvel0;
    k.angVel = angvel0;
    k.linAcc = linacc0;
    k.angAcc = angacc0;

    l = k;

    double dt = 0.001;
    k.integrate(dt);

    k0.reset();
    k1.reset();

    if(flag0 & Flags::position)
    {
      k0.position = l.position;
    }
    if(flag0 & Flags::orientation)
    {
      k0.orientation = l.orientation;
    }
    if(flag0 & Flags::linVel)
    {
      k0.linVel = l.linVel;
    }
    if(flag0 & Flags::angVel)
    {
      k0.angVel = l.angVel;
    }
    if(flag0 & Flags::linAcc)
    {
      k0.linAcc = l.linAcc;
    }
    if(flag0 & Flags::angAcc)
    {
      k0.angAcc = l.angAcc;
    }

    if(flag1 & Flags::position)
    {
      k1.position = k.position;
    }
    if(flag1 & Flags::orientation)
    {
      k1.orientation = k.orientation;
    }
    if(flag1 & Flags::linVel)
    {
      k1.linVel = k.linVel;
    }
    if(flag1 & Flags::angVel)
    {
      k1.angVel = k.angVel;
    }
    if(flag1 & Flags::linAcc)
    {
      k1.linAcc = k.linAcc;
    }
    if(flag1 & Flags::angAcc)
    {
      k1.angAcc = k.angAcc;
    }

    if((flag0 = (flag0 + 1) & Flags::all) == 0) /// update the flags to span all the possibilties
      flag1 = (flag1 + 1) & Flags::all;

    kine::Kinematics::Flags::Byte flag = BOOST_BINARY(000000);

    if(k1.position.isSet() || (k0.position.isSet() && k0.linVel.isSet() && k0.linAcc.isSet()))
    {
      flag = flag | kine::Kinematics::Flags::position;
    }
    if(k1.orientation.isSet() || (k0.orientation.isSet() && k0.angVel.isSet() && k0.angAcc.isSet()))
    {
      flag = flag | kine::Kinematics::Flags::orientation;
    }
    if(k1.linVel.isSet() || (k0.position.isSet() && k1.position.isSet()) || (k0.linAcc.isSet() && k0.linVel.isSet()))
    {
      flag = flag | kine::Kinematics::Flags::linVel;
    }
    if(k1.angVel.isSet() || (k0.orientation.isSet() && k1.orientation.isSet())
       || (k0.angVel.isSet() && k0.angAcc.isSet()))
    {
      flag = flag | kine::Kinematics::Flags::angVel;
    }
    if(k1.linAcc.isSet() || (k0.linVel.isSet() && k1.linVel.isSet())
       || (k0.linVel.isSet() && k0.position.isSet() && k1.position.isSet()))
    {
      flag = flag | kine::Kinematics::Flags::linAcc;
    }
    if(k1.angAcc.isSet() || (k0.angVel.isSet() && k1.angVel.isSet())
       || (k0.angVel.isSet() && k0.orientation.isSet() && k1.orientation.isSet()))
    {
      flag = flag | kine::Kinematics::Flags::angAcc;
    }

    k0.update(k1, dt, flag);

    if(k0.position.isSet())
    {
      err += (k.position() - k0.position()).squaredNorm();
    }
    if(k0.orientation.isSet())
    {
      err += (k.orientation.differentiate(k0.orientation)).squaredNorm();
    }
    if(k0.linVel.isSet())
    {
      if((k.linVel() - k0.linVel()).squaredNorm() < 1e-10)
      {
        err += (k.linVel() - k0.linVel()).squaredNorm();
      }
      else
      {
        err += ((k.position() - l.position()) / dt - k0.linVel()).squaredNorm();
      }
    }
    if(k0.angVel.isSet())
    {
      if((k.angVel() - k0.angVel()).squaredNorm() < 1e-10)
      {
        err += (k.angVel() - k0.angVel()).squaredNorm();
      }
      else
      {
        err += (l.orientation.differentiate(k.orientation) / dt - k0.angVel()).squaredNorm();
      }
    }
    if(k0.linAcc.isSet())
    {
      if((k.linAcc() - k0.linAcc()).squaredNorm() < 1e-10)
      {
        err += (k.linAcc() - k0.linAcc()).squaredNorm();
      }
      else
      {
        err += (k.linAcc() - 2 * k0.linAcc()).squaredNorm();
      }
    }
    if(k0.angAcc.isSet())
    {
      if((k.angAcc() - k0.angAcc()).squaredNorm() < 1e-10)
      {
        err += (k.angAcc() - k0.angAcc()).squaredNorm();
      }
      else
      {
        err += (k.angAcc() - 2 * k0.angAcc()).squaredNorm();
      }
    }

    //        std::cout<< i<<" "<<err << std::endl;

    if(err > threshold)
    {
      break;
    }
  }

  std::cout << "Error 2 : " << err << std::endl;

  if(err > threshold)
  {
    std::cout << "Error too large !" << std::endl;
    return errcode;
  }

  return 0;
}

int main()
{
  int returnVal;
  int errorcode = 0;

  if((returnVal = testRotationOperations(++errorcode)))
  {
    std::cout << "testRotationOperations Failed, error code: " << returnVal << std::endl;
    return returnVal;
  }
  else
  {
    std::cout << "testRotationOperations succeeded" << std::endl;
  }

  if((returnVal = testOrientation(++errorcode))) /// it is not an equality check
  {
    std::cout << "Orientation test failed, code : 1" << std::endl;
    return returnVal;
  }
  else
  {
    std::cout << "Orientation test succeeded" << std::endl;
  }

  if((returnVal = testKinematics(++errorcode))) /// it is not an equality check
  {
    std::cout << "Kinematics test failed, code : 2" << std::endl;
    return returnVal;
  }
  else
  {
    std::cout << "Kinematics test succeeded" << std::endl;
  }

  std::cout << "test succeeded" << std::endl;
  return 0;
}