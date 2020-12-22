#include <iostream>

#include <state-observation/dynamics-estimators/kinetics-observer.hpp>
#include <state-observation/tools/probability-law-simulation.hpp>

using namespace stateObservation;
using namespace kine;

int testOrientation(int errcode)
{

  Vector4 q_i;

  tools::ProbabilityLawSimulation s;

  /// random orientation
  q_i = s.getGaussianVector(Matrix4::Identity(), Vector4::Zero(), 4);
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
           + (ori2.matrix3() - M).norm();
    std::cout << "Assignment operaton test 2 (Matrix3) done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori3;
    ori3 = aa;

    err += (Quaternion(ori3.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm()
           + (ori3.matrix3() - M).norm();
    std::cout << "Assignment operaton test 3 (AngleAxis) done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori4;
    ori4 = Vector3(aa.angle() * aa.axis());

    err += (Quaternion(ori4.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm()
           + (ori4.matrix3() - M).norm();
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
           + (ori2.matrix3() - M).norm();
    std::cout << "copy constructor test 1 (Matrix3) done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori3(aa);

    err += (Quaternion(ori3.toVector4()).toRotationMatrix() - Quaternion(q_i).toRotationMatrix()).norm()
           + (ori3.matrix3() - M).norm();
    std::cout << "copy constructor test 2 (AngleAxis) done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori4(Vector3(aa.angle() * aa.axis()));

    err += (Quaternion(ori4.toVector4()).toRotationMatrix() - M).norm() + (ori4.matrix3() - M).norm();
    std::cout << "copy constructor test 3 (AngleAxis) done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori4(q, M);

    err += (Quaternion(ori4.toVector4()).toRotationMatrix() - M).norm() + (ori4.matrix3() - M).norm();
    std::cout << "Constructor test 1 (Quaternion, Matrix3) done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori4;
    ori4.setRandom();
    ori4 = M;

    err += (Quaternion(ori4.toVector4()).toRotationMatrix() - M).norm() + (ori4.matrix3() - M).norm();
    std::cout << "Update Assignement operator test 1 (Matrix3) done. Current error " << err << std::endl;
  }

  {
    kine::Orientation ori4;
    ori4.setRandom();
    ori4 = q;

    err += (Quaternion(ori4.toVector4()).toRotationMatrix() - M).norm() + (ori4.matrix3() - M).norm();
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

    err += (ori00.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  1 done. Current error " << err << std::endl;
    err += (ori01.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  2 done. Current error " << err << std::endl;
    err += (ori02.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  3 done. Current error " << err << std::endl;
    err += (ori03.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  4 done. Current error " << err << std::endl;
    err += (ori04.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  5 done. Current error " << err << std::endl;
    err += (ori05.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  6 done. Current error " << err << std::endl;
    err += (ori06.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  7 done. Current error " << err << std::endl;
    err += (ori07.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  8 done. Current error " << err << std::endl;
    err += (ori10.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  9 done. Current error " << err << std::endl;
    err += (ori11.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  10 done. Current error " << err << std::endl;
    err += (ori12.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  11 done. Current error " << err << std::endl;
    err += (ori13.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  12 done. Current error " << err << std::endl;
    err += (ori14.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  13 done. Current error " << err << std::endl;
    err += (ori15.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  14 done. Current error " << err << std::endl;
    err += (ori16.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  15 done. Current error " << err << std::endl;
    err += (ori17.matrix3() - Matrix3::Identity()).norm();
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

    err += (ori00.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  17 done. Current error " << err << std::endl;
    err += (ori01.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  18 done. Current error " << err << std::endl;
    err += (ori02.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  19 done. Current error " << err << std::endl;
    err += (ori03.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  20 done. Current error " << err << std::endl;
    err += (ori04.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  21 done. Current error " << err << std::endl;
    err += (ori05.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  22 done. Current error " << err << std::endl;
    err += (ori06.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  23 done. Current error " << err << std::endl;
    err += (ori07.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  24 done. Current error " << err << std::endl;
    err += (ori10.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  25 done. Current error " << err << std::endl;
    err += (ori11.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  26 done. Current error " << err << std::endl;
    err += (ori12.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  27 done. Current error " << err << std::endl;
    err += (ori13.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  28 done. Current error " << err << std::endl;
    err += (ori14.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  29 done. Current error " << err << std::endl;
    err += (ori15.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  30 done. Current error " << err << std::endl;
    err += (ori16.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  31 done. Current error " << err << std::endl;
    err += (ori17.matrix3() - Matrix3::Identity()).norm();
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

    err += (ori00.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  33 done. Current error " << err << std::endl;
    err += (ori01.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  34 done. Current error " << err << std::endl;
    err += (ori02.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  35 done. Current error " << err << std::endl;
    err += (ori03.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  36 done. Current error " << err << std::endl;
    err += (ori04.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  37 done. Current error " << err << std::endl;
    err += (ori05.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  38 done. Current error " << err << std::endl;
    err += (ori06.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  39 done. Current error " << err << std::endl;
    err += (ori07.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  40 done. Current error " << err << std::endl;
    err += (ori10.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  41 done. Current error " << err << std::endl;
    err += (ori11.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  42 done. Current error " << err << std::endl;
    err += (ori12.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  43 done. Current error " << err << std::endl;
    err += (ori13.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  44 done. Current error " << err << std::endl;
    err += (ori14.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  45 done. Current error " << err << std::endl;
    err += (ori15.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  46 done. Current error " << err << std::endl;
    err += (ori16.matrix3() - Matrix3::Identity()).norm();
    std::cout << "Multiplication test  47 done. Current error " << err << std::endl;
    err += (ori17.matrix3() - Matrix3::Identity()).norm();
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

  Vector x = o.getStateVector();
  stateObservation::TimeIndex index = o.getStateVectorSampleTime();

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
  o.addContact(contactKine, linStiffness, linDamping, angStiffness, angDamping, 3);

  Matrix12 initialCov, processCov;

  initialCov.setZero();
  initialCov.diagonal().setConstant(0.01);
  processCov.setZero();
  processCov.diagonal().setConstant(0.0001);

  contactKine.position.set() << 1, 0.1, 0;
  int i = o.addContact(contactKine, initialCov, processCov, -1);

  (void)i; /// avoid warning in release mode
  assert(i == 1);

  contactKine.position.set() << 1, -0.1, 0;
  o.addContact(contactKine, initialCov, processCov, linDamping, linStiffness, angStiffness, angDamping, 2);

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

  if((returnVal = testOrientation(1))) /// it is not an equality check
  {
    std::cout << "Orientation test failed, code : 1" << std::endl;
    return returnVal;
  }
  else
  {
    std::cout << "Orientation test succeeded" << std::endl;
  }

  if((returnVal = testKinematics(2))) /// it is not an equality check
  {
    std::cout << "Kinematics test failed, code : 2" << std::endl;
    return returnVal;
  }
  else
  {
    std::cout << "Kinematics test succeeded" << std::endl;
  }

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
