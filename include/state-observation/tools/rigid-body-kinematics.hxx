/*
 * Copyright (c) 2019-2020
 * @author Mehdi BENALLEGUE
 *
 * National Institute of Advanced Industrial Science and Technology (AIST)
 */

#include <cmath>

namespace stateObservation
{
namespace kine

{

inline void integrateKinematics(Vector3 & position, const Vector3 & velocity, double dt)
{
  position.noalias() += dt * velocity;
}

inline void integrateKinematics(Vector3 & position, Vector3 & velocity, const Vector3 & acceleration, double dt)
{
  position.noalias() += dt * velocity + 0.5 * dt * dt * acceleration;
  velocity.noalias() += dt * acceleration;
}

inline void integrateKinematics(Matrix3 & orientation, const Vector3 & rotationVelocityVector, double dt)
{
  orientation = kine::rotationVectorToRotationMatrix(rotationVelocityVector * dt) * orientation;
}

inline void integrateKinematics(Matrix3 & orientation,
                                Vector3 & rotationVelocityVector,
                                const Vector3 & rotationVelocityVectorRate,
                                double dt)
{
  orientation =
      kine::rotationVectorToRotationMatrix(rotationVelocityVector * dt + 0.5 * dt * dt * rotationVelocityVectorRate)
      * orientation;

  rotationVelocityVector += dt * rotationVelocityVectorRate;
}

inline void integrateKinematics(Quaternion & orientation, const Vector3 & rotationVelocityVector, double dt)
{
  orientation = kine::rotationVectorToQuaternion(rotationVelocityVector * dt) * orientation;
}

inline void integrateKinematics(Quaternion & orientation,
                                Vector3 & rotationVelocityVector,
                                const Vector3 & rotationVelocityVectorRate,
                                double dt)
{
  orientation =
      kine::rotationVectorToQuaternion(rotationVelocityVector * dt + 0.5 * dt * dt * rotationVelocityVectorRate)
      * orientation;

  rotationVelocityVector.noalias() += dt * rotationVelocityVectorRate;
}

inline void integrateKinematics(Vector3 & position,
                                Vector3 & velocity,
                                const Vector3 & acceleration,
                                Quaternion & orientation,
                                Vector3 & rotationVelocityVector,
                                const Vector3 & rotationVelocityVectorRate,
                                double dt)
{
  integrateKinematics(position, velocity, acceleration, dt);
  integrateKinematics(orientation, rotationVelocityVector, rotationVelocityVectorRate, dt);
}

inline void integrateKinematics(Vector3 & position,
                                Vector3 & velocity,
                                const Vector3 & acceleration,
                                Matrix3 & orientation,
                                Vector3 & rotationVelocityVector,
                                const Vector3 & rotationVelocityVectorRate,
                                double dt)
{
  integrateKinematics(position, velocity, acceleration, dt);
  integrateKinematics(orientation, rotationVelocityVector, rotationVelocityVectorRate, dt);
}

inline void integrateKinematics(Vector3 & position,
                                const Vector3 & velocity,
                                Matrix3 & orientation,
                                const Vector3 & rotationVelocity,
                                double dt)
{
  integrateKinematics(position, velocity, dt);
  integrateKinematics(orientation, rotationVelocity, dt);
}

inline void integrateKinematics(Vector3 & position,
                                const Vector3 & velocity,
                                Quaternion & orientation,
                                const Vector3 & rotationVelocity,
                                double dt)
{
  integrateKinematics(position, velocity, dt);
  integrateKinematics(orientation, rotationVelocity, dt);
}

/// Puts the orientation vector norm between 0 and Pi if it
/// gets close to 2pi
inline Vector regulateOrientationVector(const Vector3 & v)
{
  double n2 = v.squaredNorm();
  if(n2 > std::pow(3. / 2. * M_PI, 2))
  {
    double n = sqrt(n2);
    unsigned k = unsigned(ceil((n - M_PI) / (2 * M_PI)));
    return (v / n) * (n - k * 2 * M_PI);
  }
  else
  {
    return v;
  }
}

/// Transform the rotation vector into angle axis
inline AngleAxis rotationVectorToAngleAxis(const Vector3 & v)
{
  double angle(v.squaredNorm());
  if(angle > cst::epsilonAngle * cst::epsilonAngle)
  {
    angle = sqrt(angle);
    return AngleAxis(angle, v / angle);
  }
  else
  {
    return AngleAxis(0.0, Vector3::UnitZ());
  }
}

/// Tranbsform the rotation vector into rotation matrix
inline Matrix3 rotationVectorToRotationMatrix(const Vector3 & v)
{
  return (rotationVectorToAngleAxis(Vector3(v))).toRotationMatrix();
}

/// Tranbsform the rotation vector into rotation matrix
inline Quaternion rotationVectorToQuaternion(const Vector3 & v)
{
  return Quaternion(rotationVectorToAngleAxis(Vector3(v)));
}

/// Tranbsform the rotation matrix into rotation vector
inline Vector3 rotationMatrixToRotationVector(const Matrix3 & R)
{
  AngleAxis a(R);
  Vector3 v(a.axis());
  v.noalias() = a.angle() * v;
  return a.angle() * a.axis();
}

/// Tranbsform a quaternion into rotation vector
inline Vector3 quaternionToRotationVector(const Quaternion & q)
{
  AngleAxis aa(q);

  return aa.angle() * aa.axis();
}

/// Tranbsform a quaternion into rotation vector
inline Vector3 quaternionToRotationVector(const Vector4 & v)
{
  Quaternion q(v);
  AngleAxis aa(q);

  return aa.angle() * aa.axis();
}

/// Transform the rotation matrix into roll pitch yaw
///(decompose R into Ry*Rp*Rr)
inline Vector3 rotationMatrixToRollPitchYaw(const Matrix3 & R, Vector3 & v)
{
  /// source http://planning.cs.uiuc.edu/node102.html
  /// and http://planning.cs.uiuc.edu/node103.html

  //      v<<atan2(R(2,1),R(2,2)),
  //      atan2(-R(2,0),sqrt(tools::square(R(2,1))+tools::square(R(2,2)))),
  //      atan2(R(1,0),R(0,0));
  //      return v;

  return v = rotationMatrixToRollPitchYaw(R);
}

inline Vector3 rotationMatrixToRollPitchYaw(const Matrix3 & R)
{
  //      v<<atan2(R(2,1),R(2,2)),
  //      atan2(-R(2,0),sqrt(tools::square(R(2,1))+tools::square(R(2,2)))),
  //      atan2(R(1,0),R(0,0));
  //      return v;
  return R.eulerAngles(2, 1, 0).reverse();
}

/// Transform the roll pitch yaw into rotation matrix
///( R = Ry*Rp*Rr)
inline Matrix3 rollPitchYawToRotationMatrix(double roll, double pitch, double yaw)
{
  return rollPitchYawToQuaternion(roll, pitch, yaw).toRotationMatrix();
}

inline Matrix3 rollPitchYawToRotationMatrix(const Vector3 & rpy)
{
  return rollPitchYawToRotationMatrix(rpy[0], rpy[1], rpy[2]);
}

/// Transform the roll pitch yaw into rotation matrix
///( R = Ry*Rp*Rr)
inline Quaternion rollPitchYawToQuaternion(double roll, double pitch, double yaw)
{
  AngleAxis rollAngle(roll, Eigen::Vector3d::UnitX());
  AngleAxis pitchAngle(pitch, Eigen::Vector3d::UnitY());
  AngleAxis yawAngle(yaw, Eigen::Vector3d::UnitZ());

  Quaternion q(yawAngle);

  q = q * pitchAngle;
  q = q * rollAngle;

  return q;
}

inline Quaternion rollPitchYawToQuaternion(const Vector3 & rpy)
{
  return rollPitchYawToQuaternion(rpy[0], rpy[1], rpy[2]);
}

/// scalar component of a quaternion
inline double scalarComponent(const Quaternion & q)
{
  return q.coeffs()(3);
}

/// vector part of the quaternion
inline Vector3 vectorComponent(const Quaternion & q)
{
  return q.coeffs().segment<3>(0);
}

/// Projects the Matrix to so(3)
inline Matrix3 orthogonalizeRotationMatrix(const Matrix3 & M)
{
  Eigen::JacobiSVD<Matrix3> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  return svd.matrixU() * svd.matrixV().transpose();
}

/// transform a 3d vector into a skew symmetric 3x3 matrix
inline Matrix3 skewSymmetric(const Vector3 & v, Matrix3 & R)
{
  // R <<     0, -v[2],  v[1],
  //      v[2],     0, -v[0],
  //     -v[1],  v[0],     0;

  R(0, 0) = R(1, 1) = R(2, 2) = 0.;
  R(0, 1) = -(R(1, 0) = v[2]);
  R(2, 0) = -(R(0, 2) = v[1]);
  R(1, 2) = -(R(2, 1) = v[0]);

  return R;
}

/// transform a 3d vector into a skew symmetric 3x3 matrix
inline Matrix3 skewSymmetric(const Vector3 & v)
{
  Matrix3 R;

  return skewSymmetric(v, R);
}

/// transform a 3d vector into a squared skew symmetric 3x3 matrix
inline Matrix3 skewSymmetric2(const Vector3 & v, Matrix3 & R)
{
  R.noalias() = v * v.transpose();

  double n = R.trace();

  R(0, 0) -= n;
  R(1, 1) -= n;
  R(2, 2) -= n;

  return R;
}

/// transform a 3d vector into a squared skew symmetric 3x3 matrix
inline Matrix3 skewSymmetric2(const Vector3 & v)
{
  Matrix3 R;

  return skewSymmetric2(v, R);
}

/// transforms a homogeneous matrix into 6d vector (position theta mu)
inline Vector6 homogeneousMatrixToVector6(const Matrix4 & M)
{
  Vector6 v;
  AngleAxis a(AngleAxis(Matrix3(M.block(0, 0, 3, 3))));

  v.head(3) = M.block(0, 3, 3, 1);
  v.tail(3) = a.angle() * a.axis();

  return v;
}

/// transforms a 6d vector (position theta mu) into a homogeneous matrix
inline Matrix4 vector6ToHomogeneousMatrix(const Vector6 & v)
{
  Matrix4 M(Matrix4::Identity());
  M.block(0, 0, 3, 3) = rotationVectorToAngleAxis(Vector3(v.tail(3))).toRotationMatrix();
  M.block(0, 3, 3, 1) = v.head(3);
  return M;
}

inline Matrix3 twoVectorsToRotationMatrix(const Vector3 & v1, const Vector3 Rv1)
{

  BOOST_ASSERT(v1.isUnitary(cst::epsilon1) && Rv1.isUnitary(cst::epsilon1)
               && "The vectors v1 and Rv1 need to be normalized");

  double v1dotRv1 = v1.dot(Rv1);
  if(1 - v1dotRv1 < cst::epsilon1)
  {
    /// the angle is zero;
    return Matrix3::Identity();
  }
  else if(1 + v1dotRv1 < cst::epsilon1)
  {
    ///it is a singularity
    throw std::invalid_argument("twoVectorsToRotationMatrix Vectors are opposite: there is an infinity of rotation "
                                "matrices with minimal angle between them");
  }
  else
  {
    Vector3 v1xRv1 = v1.cross(Rv1);
    return Matrix3::Identity() + skewSymmetric(v1xRv1)
           + ((1 - v1.dot(Rv1)) / v1xRv1.squaredNorm()) * skewSymmetric2(v1xRv1);
  }
}

inline bool isPureYaw(const Matrix3 & R)
{
  return R.block<1, 2>(2, 0).isZero(cst::epsilon1);
}

inline Vector3 getInvariantHorizontalVector(const Matrix3 & R)
{
  if(isPureYaw(R))
  { /// Pure yaw, any vector is a solution
    return Vector3::UnitX();
  }
  else
  {
    return Vector3(R(2, 1), -R(2, 0), 0);
  }
}

inline Vector3 getInvariantOrthogonalVector(const Matrix3 & Rhat, const Vector3 & Rtez)
{
  Vector3 Rhat_Rtez = Rhat * Rtez;
  if (Rhat_Rtez.head<2>().isZero(cst::epsilon1))
  { /// any vector is a solution
    return Vector3::UnitX();
  }
  else
  {
    return Vector3(Rhat_Rtez(1), -Rhat_Rtez(0), 0);
  }
}

inline Matrix3 mergeTiltWithYaw(const Vector3 & Rtez, const Matrix3 & R2, const Vector3 & m)
{
  /*
  R&=\left(\begin{array}{ccc}
  \frac{m\times e_{z}}{\left\Vert m\times e_{z}\right\Vert } & \frac{e_{z}\times m\times e_{z}}{\left\Vert m\times
  e_{z}\right\Vert } & e_{z}\end{array}\right)\left(\begin{array}{ccc} \frac{m_{l}\times v_{1}}{\left\Vert
  m_{l}\times v_{1}\right\Vert } & \frac{v_{1}\times m_{l}\times v_{1}}{\left\Vert m_{l}\times v_{1}\right\Vert } &
  v_{1}\end{array}\right)^{T}\\&v_{1}=R_{1}^{T}e_{z}\qquad m_{l}=R_{2}^{T}m
  */
  const Vector3 & ez = Vector3::UnitZ();

  const Vector3 & v1 = Rtez;

  Vector3 mlxv1 = (R2.transpose() * m).cross(v1);

  double n2 = mlxv1.squaredNorm();

  if(n2 > cst::epsilonAngle * cst::epsilonAngle)
  {
    /// ezxmxez = ez.cross(mxez);
    /// R_temp1 << mxez, ezxmxez, ez;

    Matrix3 R_temp1, R_temp2;
    R_temp1 << m.cross(ez), m, ez;

    mlxv1 /= sqrt(n2);

    // clang-format off
        R_temp2 << mlxv1.transpose(), 
                   v1.cross(mlxv1).transpose(),
                   v1.transpose();
    // clang-format on

    return R_temp1 * R_temp2;
  }
  else
  {
    throw std::invalid_argument(
        "Tilt is singular with regard to the provided axis (likely gimbal lock). Please use the "
        "angle-agnostic version");
  }
}

inline Matrix3 mergeRoll1Pitch1WithYaw2(const Matrix3 & R1, const Matrix3 & R2, const Vector3 v)
{
  return mergeTiltWithYaw(R1.transpose() * Vector3::UnitZ(), R2, v);
}

inline Matrix3 mergeTiltWithYawAxisAgnostic(const Vector3 & Rtez, const Matrix3 & R2)
{
  /*
 R&=\left(\begin{array}{ccc}
 \frac{m\times e_{z}}{\left\Vert m\times e_{z}\right\Vert } & \frac{e_{z}\times m\times e_{z}}{\left\Vert m\times
 e_{z}\right\Vert } & e_{z}\end{array}\right)\left(\begin{array}{ccc} \frac{m_{l}\times v_{1}}{\left\Vert
 m_{l}\times v_{1}\right\Vert } & \frac{v_{1}\times m_{l}\times v_{1}}{\left\Vert m_{l}\times v_{1}\right\Vert } &
 v_{1}\end{array}\right)^{T}\\&v_{1}=R_{1}^{T}e_{z}\qquad m_{l}=R_{2}^{T}m
 */
  const Vector3 & ez = Vector3::UnitZ();

  const Vector3 & v1 = Rtez;

  Vector3 m = getInvariantOrthogonalVector(R2, Rtez).normalized();

  Vector3 ml = R2.transpose() * m;

  Matrix3 R_temp1, R_temp2;
  R_temp1 << m.cross(ez), m, ez;

  // clang-format off
  R_temp2 << ml.cross(v1).transpose(),
             ml.transpose(),
             v1.transpose();
  // clang-format on

  return R_temp1 * R_temp2;
}

inline Matrix3 mergeRoll1Pitch1WithYaw2AxisAgnostic(const Matrix3 & R1, const Matrix3 & R2)
{
  return mergeTiltWithYawAxisAgnostic(R1.transpose() * Vector3::UnitZ(), R2);
}

inline double rotationMatrixToAngle(const Matrix3 & rotation, const Vector3 & axis, const Vector3 & v)
{
  Vector3 rotV_proj = axis.cross((rotation * v).cross(axis)); /// no need to normalize for atan2
  return atan2(v.cross(rotV_proj).dot(axis), v.dot(rotV_proj));
}

inline double rotationMatrixToYaw(const Matrix3 & rotation, const Vector2 & v)
{
  Vector2 rotV = (rotation.topLeftCorner<2, 2>() * v);
  return atan2(v.x() * rotV.y() - v.y() * rotV.x(), v.dot(rotV));
}

inline double rotationMatrixToYaw(const Matrix3 & rotation)
{
  return atan2(rotation(1, 0), rotation(0, 0));
}

inline double rotationMatrixToYawAxisAgnostic(const Matrix3 & rotation)
{

  if(isPureYaw(rotation))
  { /// this is the case of pure yaw rotationm, so simply extract yaw
    return rotationMatrixToYaw(rotation);
  }
  else
  { /// this is the case where there are non-yaw rotations (roll & pitch)
    /// v is the vector that stays horizontal by rotation
    Vector2 v = getInvariantHorizontalVector(rotation).head<2>();
    Vector2 rotV = (rotation.topLeftCorner<2, 2>() * v);
    return atan2(v.x() * rotV.y() - v.y() * rotV.x(), v.dot(rotV));
  }
}

inline Quaternion zeroRotationQuaternion()
{
  return Quaternion(1, 0, 0, 0);
}

inline Quaternion randomRotationQuaternion()
{
  return Quaternion(tools::ProbabilityLawSimulation::getGaussianVector<Vector4>(Matrix4::Identity(), Vector4::Zero(), 4).normalized());
}

inline double randomAngle()
{
  return tools::ProbabilityLawSimulation::getUniformScalar(-M_PI, M_PI);
}

inline bool isRotationMatrix(const Matrix3 & M, double precision)
{
  return (M.isUnitary()
          && (isApproxAbs(M.topLeftCorner<2, 2>().determinant(), M(2, 2) * precision)
              || isApproxAbs(M.bottomLeftCorner<2, 2>().determinant(), M(0, 2), precision)));
}

/// transforms a rotation into translation given a constraint of a fixed point
/// which means the global position of the fixed point is constantly at its
/// constrained position
inline void fixedPointRotationToTranslation(const Matrix3 & R,
                                            const Vector3 & rotationVelocity,
                                            const Vector3 & rotationAcceleration,
                                            const Vector3 & fixedPoint,
                                            Vector3 & outputTranslation,
                                            Vector3 & outputLinearVelocity,
                                            Vector3 & outputLinearAcceleration)
{
  Vector3 localFixed = R * fixedPoint;
  outputTranslation = fixedPoint - localFixed;
  outputLinearVelocity = -rotationVelocity.cross(localFixed);
  outputLinearAcceleration =
      -rotationVelocity.cross(rotationVelocity.cross(localFixed)) - rotationAcceleration.cross(localFixed);
}

/// Computes the multiplicative Jacobians for
/// Kalman filtering for examole
inline void derivateRotationMultiplicative(const Vector3 & deltaR, Matrix3 & dRdR, Matrix3 & dRddeltaR)
{
  dRdR = rotationVectorToRotationMatrix(deltaR);
  double d = deltaR.norm();

  if(d < cst::epsilonAngle)
  {
    dRddeltaR = Matrix3::Identity();
  }
  else
  {
    double cosd = cos(d);
    double sind = sin(d);

    double d2 = d * d;
    double d3 = d2 * d;

    dRddeltaR.noalias() = deltaR * deltaR.transpose() * (d - sind) / d3 + Matrix3::Identity() * sind / d
                          + skewSymmetric(deltaR) * (1 - cosd) / d2;
  }
}

/// Computes the multiplicative Jacobians for a vector expressed in
/// local frame
inline Matrix3 derivateRtvMultiplicative(const Matrix3 & R, const Vector3 & v)
{
  return -R.transpose() * skewSymmetric(v);
}

/// derivates a quaternion using finite difference to get a angular velocity vector
inline Vector3 derivateRotationFD(const Quaternion & q1, const Quaternion & q2, double dt)
{
  AngleAxis aa(q2 * q1.conjugate());

  return (aa.angle() / dt) * aa.axis();
}

/// derivates a quaternion using finite difference to get a angular velocity vector
inline Vector3 derivateRotationFD(const Vector3 & o1, const Vector3 & o2, double dt)
{
  Quaternion q1(rotationVectorToAngleAxis(o1));
  Quaternion q2(rotationVectorToAngleAxis(o2));

  return derivateRotationFD(q1, q2, dt);
}

inline Vector6 derivateHomogeneousMatrixFD(const Matrix4 & m1, const Matrix4 & m2, double dt)
{
  Vector6 out;

  Matrix3 r1 = m1.block(0, 0, 3, 3);
  Matrix3 r2 = m2.block(0, 0, 3, 3);

  AngleAxis aa(r2 * r1.transpose());

  out.tail(3) = (aa.angle() / dt) * aa.axis();

  out.head(3) = (m2.block(0, 3, 3, 1) - m1.block(0, 3, 3, 1)) / dt;

  return out;
}

inline Vector6 derivatePoseThetaUFD(const Vector6 & v1, const Vector6 & v2, double dt)
{
  Vector6 out;

  out.tail(3) = derivateRotationFD(v1.tail(3), v2.tail(3), dt);

  out.head(3) = (v2.head(3) - v1.head(3)) / dt;

  return out;
}

/// uses the derivation to reconstruct the velocities and accelerations given
/// trajectories in positions and orientations only
inline IndexedVectorArray reconstructStateTrajectory(const IndexedVectorArray & positionOrientation, double dt)
{
  typedef kine::indexes<kine::rotationVector> indexes;

  Vector r(Vector::Zero(18, 1));

  const IndexedVectorArray & po = positionOrientation;

  TimeIndex i0 = positionOrientation.getFirstIndex();
  TimeIndex i1 = positionOrientation.getNextIndex();

  IndexedVectorArray a;
  a.setValue(r, i0);
  a.resize(po.size(), r);

  for(TimeIndex i = i0; i < i1; ++i)
  {
    Vector poi = po[i];

    r.segment(indexes::pos, 3) = poi.head(3);
    r.segment(indexes::ori, 3) = poi.tail(3);
    a.setValue(r, i);
  }

  for(TimeIndex i = i0; i < i1 - 1; ++i)
  {
    r = a[i];

    Vector poi = po[i];
    Vector poi1 = po[i + 1];

    r.segment(indexes::linVel, 3) = tools::derivate(Vector3(poi.head(3)), Vector3(poi1.head(3)), dt);
    r.segment(indexes::angVel, 3) = derivateRotationFD(Quaternion(rotationVectorToAngleAxis(poi.tail(3))),
                                                       Quaternion(rotationVectorToAngleAxis(poi1.tail(3))), dt);

    a.setValue(r, i);
  }

  for(TimeIndex i = i0; i < i1 - 2; ++i)
  {
    r = a[i];
    Vector r2 = a[i + 1];

    r.segment(indexes::linAcc, 3) =
        tools::derivate(Vector3(r.segment(indexes::linVel, 3)), Vector3(r2.segment(indexes::linVel, 3)), dt);
    r.segment(indexes::angAcc, 3) =
        tools::derivate(Vector3(r.segment(indexes::angVel, 3)), Vector3(r2.segment(indexes::angVel, 3)), dt);

    a.setValue(r, i);
  }

  return a;
}

inline Vector invertState(const Vector & state)
{
  typedef kine::indexes<kine::rotationVector> indexes;
  Matrix3 r2 = (rotationVectorToAngleAxis(-state.segment(indexes::ori, 3))).toRotationMatrix(); // inverse
  Vector3 omega1 = state.segment(indexes::angVel, 3);
  Vector3 omega1dot = state.segment(indexes::angAcc, 3);
  Matrix3 omega1x = skewSymmetric(omega1);
  Matrix3 omega1dotx = skewSymmetric(omega1dot);

  Vector3 t1 = state.segment(indexes::pos, 3);
  Vector3 t1dot = state.segment(indexes::linVel, 3);
  Vector3 t1dotdot = state.segment(indexes::linAcc, 3);

  Vector state2(Vector::Zero(18, 1));
  state2.segment(indexes::pos, 3) = -r2 * t1; // t2
  state2.segment(indexes::linVel, 3) = r2 * (omega1x * t1 - t1dot); // t2dot
  state2.segment(indexes::linAcc, 3) =
      r2 * (omega1x * (2 * t1dot - omega1x * t1) - t1dotdot + omega1dotx * t1); // t2dotdot
  state2.segment(indexes::ori, 3) = -state.segment(indexes::ori, 3); // thetaU2
  state2.segment(indexes::angVel, 3) = -r2 * omega1; // omega2
  state2.segment(indexes::angAcc, 3) = r2 * (omega1x * omega1 - omega1dot); // omega2dot
  return state2;
}

inline Matrix4 invertHomoMatrix(const Matrix4 & m)
{
  Matrix4 m2(Matrix4::Identity());
  Matrix3 rt = m.block(0, 0, 3, 3).transpose();
  m2.block(0, 0, 3, 3) = rt;
  m2.block(0, 3, 3, 1) = -rt * m.block(0, 3, 3, 1);
  return m2;
}

inline Orientation::Orientation(bool b) : q_(b), m_(b) {}

inline Orientation::Orientation(const Vector3 & v) : q_(rotationVectorToQuaternion(v)) {}

inline Orientation::Orientation(const Quaternion & q) : q_(q) {}

inline Orientation::Orientation(const Matrix3 & m) : m_(m) {}

inline Orientation::Orientation(const AngleAxis & aa) : q_(Quaternion(aa)) {}

inline Orientation::Orientation(const double & roll, const double & pitch, const double & yaw)
: q_(kine::rollPitchYawToQuaternion(roll, pitch, yaw))
{
}

inline Orientation::Orientation(const Quaternion & q, const Matrix3 & m) : q_(q), m_(m) {}

inline Orientation::Orientation(const Orientation & operand1, const Orientation & operand2) : q_(false), m_(false)
{
  setToProductNoAlias(operand1, operand2);
}

inline Orientation & Orientation::operator=(const Vector3 & v)
{
  m_.reset();
  q_ = rotationVectorToQuaternion(v);
  return *this;
}

inline Orientation & Orientation::operator=(const Quaternion & q)
{
  m_.reset();
  q_ = q;
  return *this;
}

inline Orientation & Orientation::operator=(const Matrix3 & m)
{
  q_.reset();
  m_ = m;
  return *this;
}

inline Orientation & Orientation::operator=(const AngleAxis & aa)
{
  m_.reset();
  q_ = Quaternion(aa);
  return *this;
}

inline Orientation & Orientation::setValue(const Quaternion & q, const Matrix3 & m)
{
  q_ = q;
  m_ = m;
  return *this;
}

inline Orientation & Orientation::fromVector4(const Vector4 & v)
{
  return (*this) = Quaternion(v);
}

inline Orientation & Orientation::setRandom()
{
  return (*this) = randomRotationQuaternion();
}

template<>
inline Orientation & Orientation::setZeroRotation<Quaternion>()
{
  return (*this) = zeroRotationQuaternion();
}

template<>
inline Orientation & Orientation::setZeroRotation<Matrix3>()
{
  return (*this) = Matrix3(Matrix3::Identity());
}

inline Orientation & Orientation::setZeroRotation()
{
  return setZeroRotation<Quaternion>();
}

inline Vector4 Orientation::toVector4() const
{
  return toQuaternion().coeffs();
}

inline Vector3 Orientation::toRotationVector() const
{
  check_();
  if(isQuaternionSet())
  {
    return kine::quaternionToRotationVector(q_());
  }
  else
  {
    return kine::rotationMatrixToRotationVector(m_());
  }
}

inline Vector3 Orientation::toRollPitchYaw() const
{
  check_();
  if(!isMatrixSet())
  {
    quaternionToMatrix_();
  }
  return kine::rotationMatrixToRollPitchYaw(m_());
}

inline AngleAxis Orientation::toAngleAxis() const
{
  check_();
  if(isMatrixSet())
  {
    return AngleAxis(m_());
  }
  else
  {
    return AngleAxis(q_());
  }
}

inline const Matrix3 & Orientation::toMatrix3() const
{
  check_();
  if(!isMatrixSet())
  {
    quaternionToMatrix_();
  }
  return m_;
}

inline const Quaternion & Orientation::toQuaternion() const
{
  check_();
  if(!isQuaternionSet())
  {
    matrixToQuaternion_();
  }
  return q_;
}

inline Orientation::operator const Matrix3 &() const
{
  return toMatrix3();
}

inline Orientation::operator const Quaternion &() const
{
  return toQuaternion();
}

inline const Orientation & Orientation::setToProductNoAlias(const Orientation & R1, const Orientation & R2)
{
  R1.check_();
  R2.check_();
  if(R1.isQuaternionSet() && R2.isQuaternionSet())
  {
    if(R1.isMatrixSet() && R2.isMatrixSet())
    {
      m_.set().noalias() = R1.m_() * R2.m_();
    }
    else
    {
      m_.reset();
    }
    q_ = R1.q_() * R2.q_();
  }
  else
  {
    m_.set(true); /// we set the matrix as initialized before giving the value
    if(!R1.isMatrixSet())
    {
      m_().noalias() = R1.quaternionToMatrix_() * R2.m_();
    }
    else if(!R2.isMatrixSet())
    {
      m_().noalias() = R1.m_() * R2.quaternionToMatrix_();
    }
    else
    {
      m_().noalias() = R1.m_() * R2.m_();
    }
    q_.reset();
  }

  return *this;
}

inline Orientation Orientation::operator*(const Orientation & R2) const
{
  return Orientation(*this, R2);
}

inline Orientation Orientation::inverse() const
{
  check_();
  if(q_.isSet())
  {
    if(isMatrixSet())
    {
      return Orientation(q_().inverse(), m_().transpose());
    }
    else
    {
      return Orientation(q_().inverse());
    }
  }
  else
  {
    return Orientation(Matrix3(m_().transpose()));
  }
}

inline const Orientation & Orientation::integrate(Vector3 dt_x_omega)
{
  check_();
  if(q_.isSet())
  {
    if(isMatrixSet())
    {
      Quaternion q = kine::rotationVectorToQuaternion(dt_x_omega);
      q_ = q * q_();
      m_ = q.toRotationMatrix() * m_();
    }
    else
    {
      q_ = kine::rotationVectorToQuaternion(dt_x_omega) * q_();
    }
  }
  else
  {
    m_ = kine::rotationVectorToRotationMatrix(dt_x_omega) * m_();
  }
  return *this;
}

inline Vector3 Orientation::differentiate(Orientation R_k1) const
{
  check_();
  return (Orientation(R_k1, inverse())).toRotationVector();
}

inline bool Orientation::isSet() const
{
  return (isMatrixSet() || q_.isSet());
}

inline void Orientation::synchronize()
{
  check_();
  if(isMatrixSet())
  {
    if(!isQuaternionSet())
    {
      matrixToQuaternion_();
    }
  }
  else
  {
    if(isQuaternionSet())
    {
      quaternionToMatrix_();
    }
  }
}

inline void Orientation::reset()
{
  q_.reset();
  m_.reset();
}

inline bool Orientation::isMatrixSet() const
{
  return (m_.isSet());
}

inline bool Orientation::isQuaternionSet() const
{
  return (q_.isSet());
}

inline void Orientation::setMatrix(bool b)
{
  m_.set(b);
}

inline void Orientation::setQuaternion(bool b)
{
  q_.set(b);
}

inline Vector3 Orientation::operator*(const Vector3 & v) const
{
  check_();
  if(!isMatrixSet())
  {
    quaternionToMatrix_();
  }
  return m_() * v;
}

inline CheckedMatrix3 & Orientation::getMatrixRefUnsafe()
{
  return m_;
}

inline CheckedQuaternion & Orientation::getQuaternionRefUnsafe()
{
  return q_;
}

inline void Orientation::check_() const
{
  BOOST_ASSERT((isQuaternionSet() || isMatrixSet()) && "The orientation is not initialized");
}

inline const Matrix3 & Orientation::quaternionToMatrix_() const
{
  return m_ = q_().toRotationMatrix();
}

inline const Quaternion & Orientation::matrixToQuaternion_() const
{
  return q_ = Quaternion(m_());
}

inline Orientation Orientation::zeroRotation()
{
  return Orientation(Quaternion::Identity(), Matrix3::Identity());
}

inline Orientation Orientation::randomRotation()
{
  Orientation o;
  o.setRandom();
  return o;
}

///////////////////////////////////////////////////////////////////////
/// -------------------Kinematics structure implementation-------------
///////////////////////////////////////////////////////////////////////

inline Kinematics::Kinematics(const Vector & v, Kinematics::Flags::Byte flags)
{
  fromVector(v, flags);
}

inline Kinematics::Kinematics(const Kinematics & multiplier1, const Kinematics & multiplier2)
{
  setToProductNoAlias(multiplier1, multiplier2);
}

inline Kinematics & Kinematics::fromVector(const Vector & v, Kinematics::Flags::Byte flags)
{
  int index = 0;

  bool flagPos = flags & Flags::position;
  bool flagLinVel = flags & Flags::linVel;
  bool flagLinAcc = flags & Flags::linAcc;
  bool flagOri = flags & Flags::orientation;
  bool flagAngVel = flags & Flags::angVel;
  bool flagAngAcc = flags & Flags::angAcc;

  if(flagPos)
  {
    BOOST_ASSERT(v.size() >= index + 3 && "The kinematics vector size is incorrect (loading position)");
    if(v.size() >= index + 3)
    {
      position = v.segment<3>(index);
      index += 3;
    }
  }

  if(flagOri)
  {
    BOOST_ASSERT(v.size() >= index + 4 && "The kinematics vector size is incorrect (loading orientaTion)");
    if(v.size() >= index + 4)
    {
      orientation.fromVector4(v.segment<4>(index));
      index += 4;
    }
  }

  if(flagLinVel)
  {
    BOOST_ASSERT(v.size() >= index + 3 && "The kinematics vector size is incorrect (loading linear velocity)");
    if(v.size() >= index + 3)
    {
      linVel = v.segment<3>(index);
      index += 3;
    }
  }

  if(flagAngVel)
  {
    BOOST_ASSERT(v.size() >= index + 3 && "The kinematics vector size is incorrect (loading angular velocity)");
    if(v.size() >= index + 3)
    {
      angVel = v.segment<3>(index);
      index += 3;
    }
  }

  if(flagLinAcc)
  {
    BOOST_ASSERT(v.size() >= index + 3 && "The kinematics vector size is incorrect (loading linear acceleration)");
    if(v.size() >= index + 3)
    {
      linAcc = v.segment<3>(index);
      index += 3;
    }
  }

  if(flagAngAcc)
  {
    BOOST_ASSERT(v.size() >= index + 3 && "The kinematics vector size is incorrect (loading angular acceleration)");
    if(v.size() >= index + 3)
    {
      angAcc = v.segment<3>(index);
      // index+=3; ///useless
    }
  }

  return *this;
}

template<typename t>
inline Kinematics & Kinematics::setZero(Kinematics::Flags::Byte flags)
{

  bool flagPos = flags & Flags::position;
  bool flagLinVel = flags & Flags::linVel;
  bool flagLinAcc = flags & Flags::linAcc;
  bool flagOri = flags & Flags::orientation;
  bool flagAngVel = flags & Flags::angVel;
  bool flagAngAcc = flags & Flags::angAcc;

  if(flagPos)
  {
    position.set().setZero();
  }

  if(flagOri)
  {
    orientation.setZeroRotation<t>();
  }

  if(flagLinVel)
  {
    linVel.set().setZero();
  }

  if(flagAngVel)
  {
    angVel.set().setZero();
  }

  if(flagLinAcc)
  {
    linAcc.set().setZero();
  }

  if(flagAngAcc)
  {
    angAcc.set().setZero();
  }

  return *this;
}

inline Kinematics & Kinematics::setZero(Kinematics::Flags::Byte flags)
{
  return setZero<Quaternion>(flags);
}

inline const Kinematics & Kinematics::integrate(double dt)
{
  if(angVel.isSet())
  {
    if(angAcc.isSet())
    {
      if(orientation.isSet())
      {
        orientation.integrate(angVel() * dt + angAcc() * dt * dt / 2);
      }
      angVel() += angAcc() * dt;
    }
    else
    {
      if(orientation.isSet())
      {
        orientation.integrate(angVel() * dt);
      }
    }
  }

  if(linVel.isSet())
  {
    if(linAcc.isSet())
    {
      if(position.isSet())
      {
        position() += linVel() * dt + linAcc() * dt * dt / 2;
      }
      linVel() += linAcc() * dt;
    }
    else
    {
      if(position.isSet())
      {
        position() += linVel() * dt;
      }
    }
  }

  return *this;
}

inline const Kinematics & Kinematics::update(const Kinematics & newValue, double dt, Flags::Byte flags)
{
  {
    bool flagPos = flags & Flags::position;
    bool flagLinVel = flags & Flags::linVel;
    bool flagLinAcc = flags & Flags::linAcc;

    CheckedVector3 & thisPos = position;
    CheckedVector3 & thisVel = linVel;
    CheckedVector3 & thisAcc = linAcc;

    const CheckedVector3 & newPos = newValue.position;
    const CheckedVector3 & newVel = newValue.linVel;
    const CheckedVector3 & newAcc = newValue.linAcc;

    enum method
    {
      noUpdate,
      usePosition,
      useVelocity,
      useAcceleration,
      useVelAndAcc,
      usePosAndVel
    };

    method posMethod = noUpdate;
    method velMethod = noUpdate;
    method accMethod = noUpdate;

    if(flagPos)
    {
      if(newPos.isSet())
      {
        posMethod = usePosition;
      }
      else
      {
        BOOST_ASSERT(thisPos.isSet() && "The position cannot be updated without initial value");
        if(thisVel.isSet())
        {
          if(thisAcc.isSet())
          {
            posMethod = useVelAndAcc;
          }
          else
          {
            posMethod = useVelocity;
          }
        }
        else
        {
          if(thisAcc.isSet())
          {
            posMethod = useAcceleration;
          }
        }
      }
    }

    if(flagLinVel)
    {
      if(newVel.isSet())
      {
        velMethod = useVelocity;
      }
      else
      {
        if(thisPos.isSet() && newPos.isSet())
        {
          velMethod = usePosition;
        }
        else
        {
          BOOST_ASSERT(thisVel.isSet() && "The linear velocity cannot be updated without initial value");
          if(thisAcc.isSet())
          {
            velMethod = useAcceleration;
          }
        }
      }
    }

    if(flagLinAcc)
    {
      if(newAcc.isSet())
      {
        accMethod = useAcceleration;
      }
      else
      {
        if(thisVel.isSet() && newVel.isSet())
        {
          accMethod = useVelocity;
        }
        else
        {
          if(thisVel.isSet() && velMethod == usePosition)
          {
            accMethod = usePosAndVel;
          }

          else /// velocity is not available
          {
            if(thisPos.isSet() && newPos.isSet())
            {
              accMethod = usePosition;
            }
            else
            {
              BOOST_ASSERT(thisAcc.isSet() && "The linear accleration cannot be updated without initial value");
            }
          }
        }
      }
    }

    if(accMethod == useVelocity) // then velocity cannot use accelerations
    {
      thisAcc = (newVel() - thisVel()) / dt;
    }
    else
    {
      if(accMethod == usePosition) // then velocity cannot use accelerations
      {
        thisAcc = 2 * (newPos() - thisPos()) / (dt * dt);
      }
    }

    if(velMethod == usePosition) // then position cannot use velocity
    {
      if(accMethod == usePosAndVel) // use position and velocity to get the accelerations
      {
        thisAcc = -thisVel();
        thisVel = (newPos() - thisPos()) / dt;
        thisAcc() += thisVel();
        thisAcc() /= dt;
      }
      else
      {
        thisVel = (newPos() - thisPos()) / dt;
      }
    }

    if(posMethod == usePosition)
    {
      thisPos = newPos;
    }
    else
    {
      if(posMethod == useVelocity)
      {
        thisPos() += thisVel() * dt;
      }
      else
      {
        if(posMethod == useVelAndAcc)
        {
          thisPos() += thisVel() * dt + thisAcc() * dt * dt / 2;
        }
        else
        {
          if(posMethod == useAcceleration)
          {
            thisPos() += thisAcc() * dt * dt / 0.5;
          }
        }
      }
    }

    if(velMethod == useVelocity)
    {
      thisVel = newVel;
    }
    else
    {
      if(velMethod == useAcceleration)
      {
        thisVel() += thisAcc() * dt;
      }
    }

    if(accMethod == useAcceleration)
    {
      thisAcc = newAcc();
    }

    if(!flagPos)
    {
      thisPos.reset();
    }
    if(!flagLinVel)
    {
      thisVel.reset();
    }
    if(!flagLinAcc)
    {
      thisAcc.reset();
    }
  }

  {
    bool flagOri = flags & Flags::orientation;
    bool flagAngVel = flags & Flags::angVel;
    bool flagAngAcc = flags & Flags::angAcc;

    Orientation & thisOri = orientation;
    CheckedVector3 & thisVel = angVel;
    CheckedVector3 & thisAcc = angAcc;

    const Orientation & newOri = newValue.orientation;
    const CheckedVector3 & newVel = newValue.angVel;
    const CheckedVector3 & newAcc = newValue.angAcc;

    enum method
    {
      noUpdate,
      useOrientation,
      useVelocity,
      useAcceleration,
      useVelAndAcc,
      useOriAndVel
    };

    method posMethod = noUpdate;
    method velMethod = noUpdate;
    method accMethod = noUpdate;

    if(flagOri)
    {
      if(newOri.isSet())
      {
        posMethod = useOrientation;
      }
      else
      {
        BOOST_ASSERT(thisOri.isSet() && "The orientation is trying to be updated without initial value");
        if(thisVel.isSet())
        {
          if(thisAcc.isSet())
          {
            posMethod = useVelAndAcc;
          }
          else
          {
            posMethod = useVelocity;
          }
        }
        else
        {
          if(thisAcc.isSet())
          {
            posMethod = useAcceleration;
          }
        }
      }
    }

    if(flagAngVel)
    {
      if(newVel.isSet())
      {
        velMethod = useVelocity;
      }
      else
      {
        if(thisOri.isSet() && newOri.isSet())
        {
          velMethod = useOrientation;
        }
        else
        {
          BOOST_ASSERT(thisVel.isSet() && "The angular velocity is trying to be updated without initial value");

          if(thisAcc.isSet())
          {
            velMethod = useAcceleration;
          }
        }
      }
    }

    if(flagAngAcc)
    {
      if(newAcc.isSet())
      {
        accMethod = useAcceleration;
      }
      else
      {
        if(thisVel.isSet() && newVel.isSet())
        {
          accMethod = useVelocity;
        }
        else
        {
          if(thisVel.isSet() && velMethod == useOrientation)
          {
            accMethod = useOriAndVel;
          }
          else /// velocity is not available
          {
            if(thisOri.isSet() && newOri.isSet())
            {
              accMethod = useOrientation;
            }
            else
            {
              BOOST_ASSERT(thisAcc.isSet() && "The angular accleration is trying to be updated without initial value");
            }
          }
        }
      }
    }

    if(accMethod == useVelocity) // then velocity cannot use accelerations
    {
      thisAcc = (newVel() - thisVel()) / dt;
    }
    else
    {
      if(accMethod == useOrientation) // then velocity cannot use accelerations
      {
        thisAcc = 2 * newOri.differentiate(thisOri) / (dt * dt);
      }
    }

    if(velMethod == useOrientation) // then position cannot use velocity
    {
      if(accMethod == useVelAndAcc)
      {
        thisAcc = -thisVel();
        thisVel = thisOri.differentiate(newOri) / dt;
        thisAcc() += thisVel();
        thisAcc() /= dt;
      }
      else
      {
        thisVel = thisOri.differentiate(newOri) / dt;
      }
    }

    if(posMethod == useOrientation)
    {
      thisOri = newOri;
    }
    else
    {
      if(posMethod == useVelocity)
      {
        thisOri.integrate(thisVel() * dt);
      }
      else
      {
        if(posMethod == useVelAndAcc)
        {
          thisOri.integrate(thisVel() * dt + thisAcc() * dt * dt / 2);
        }
        else
        {
          if(posMethod == useAcceleration)
          {
            thisOri.integrate(thisAcc() * dt * dt / 0.5);
          }
        }
      }
    }

    if(velMethod == useVelocity)
    {
      thisVel = newVel;
    }
    else
    {
      if(velMethod == useAcceleration)
      {
        thisVel() += thisAcc() * dt;
      }
    }

    if(accMethod == useAcceleration)
    {
      thisAcc = newAcc();
    }

    if(!flagOri)
    {
      thisOri.reset();
    }
    if(!flagAngVel)
    {
      thisVel.reset();
    }
    if(!flagAngAcc)
    {
      thisAcc.reset();
    }
  }
  return *this;
}

inline Kinematics Kinematics::getInverse() const
{
  Kinematics inverted;

  BOOST_ASSERT(orientation.isSet() && "The orientation is not initialized, the kinematics cannot be inverted.");

  inverted.orientation = orientation.inverse();
  Orientation & r2 = inverted.orientation;

  if(angVel.isSet())
  {
    inverted.angVel = -(r2 * angVel()); // omega2

    if(angAcc.isSet())
    {
      inverted.angAcc = r2 * (angVel().cross(angVel()) - angAcc()); // omega2dot
    }
  }

  if(position.isSet())
  {
    inverted.position = -(r2 * position());

    if(linVel.isSet() && angVel.isSet())
    {

      Vector3 omegaxp = angVel().cross(position());
      inverted.linVel = r2 * (omegaxp - linVel()); // t2dot
      if(linAcc.isSet() && (angAcc.isSet()))
      {

        inverted.linAcc =
            r2 * (angVel().cross(2 * linVel - omegaxp) - linAcc() + angAcc().cross(position())); // t2dotdot
      }
    }
  }

  return inverted;
}

/// composition of transformation
inline Kinematics Kinematics::operator*(const Kinematics & multiplier) const
{
  return Kinematics(*this, multiplier);
}

inline Kinematics Kinematics::setToProductNoAlias(const Kinematics & multiplier1, const Kinematics & multiplier2)
{
  BOOST_ASSERT(multiplier1.orientation.isSet()
               && "The multiplier 1 orientation is not initialized, the multiplication is not possible.");

  BOOST_ASSERT((multiplier2.position.isSet() || multiplier2.orientation.isSet())
               && "The multiplier 2 kinematics is not initialized, the multiplication is not possible.");

  if(multiplier2.position.isSet() && multiplier1.position.isSet())
  {
    position.set(true);
    Vector3 & R1p2 = position(); /// reference ( Vector3&  )
    R1p2.noalias() = multiplier1.orientation * multiplier2.position();

    if(multiplier2.linVel.isSet() && multiplier1.linVel.isSet() && multiplier1.angVel.isSet())
    {
      Vector3 & R1p2d = tempVec_; /// reference
      R1p2d.noalias() = multiplier1.orientation * multiplier2.linVel();

      linVel.set(true);
      Vector3 & w1xR1p2 = linVel(); /// reference
      w1xR1p2.noalias() = multiplier1.angVel().cross(R1p2);

      Vector3 & w1xR1p2_R1p2d = w1xR1p2; ///  reference ( =linVel() )
      w1xR1p2_R1p2d += R1p2d;

      if(multiplier2.linAcc.isSet() && multiplier1.linAcc.isSet() && multiplier1.angAcc.isSet())
      {
        linAcc.set(true);
        linAcc().noalias() = multiplier1.orientation * multiplier2.linAcc();
        linAcc().noalias() += multiplier1.angAcc().cross(R1p2);
        linAcc().noalias() += multiplier1.angVel().cross(w1xR1p2_R1p2d + R1p2d);
        linAcc() += multiplier1.linAcc();
      }
      else
      {
        linAcc.reset();
      }

      linVel() += multiplier1.linVel();
    }
    else
    {
      linVel.reset();
      linAcc.reset();
    }

    position() += multiplier1.position();
  }
  else
  {
    position.reset();
    linVel.reset();
    linAcc.reset();
  }

  if(multiplier2.orientation.isSet())
  {
    orientation.setToProductNoAlias(multiplier1.orientation, multiplier2.orientation);

    if(multiplier2.angVel.isSet() && multiplier1.angVel.isSet())
    {
      angVel.set(true);
      Vector3 & R1w2 = angVel(); /// reference
      R1w2.noalias() = multiplier1.orientation * multiplier2.angVel();

      if(multiplier2.angAcc.isSet() && multiplier1.angAcc.isSet())
      {
        angAcc.set(true);
        angAcc().noalias() = multiplier1.orientation * multiplier2.angAcc();
        angAcc().noalias() += multiplier1.angVel().cross(R1w2);
        angAcc() += multiplier1.angAcc();
      }
      else
      {
        angAcc.reset();
      }

      angVel() += multiplier1.angVel();
    }
    else
    {
      angVel.reset();
      angAcc.reset();
    }
  }
  else
  {
    orientation.reset();
    angVel.reset();
    angAcc.reset();
  }

  return *this;
}

inline Vector Kinematics::toVector(Flags::Byte flags) const
{
  int size = 0;

  if(flags & Flags::position)
  {
    size += 3;
  }
  if(flags & Flags::orientation)
  {
    size += 4;
  }
  if(flags & Flags::linVel)
  {
    size += 3;
  }
  if(flags & Flags::angVel)
  {
    size += 3;
  }
  if(flags & Flags::linAcc)
  {
    size += 3;
  }
  if(flags & Flags::angAcc)
  {
    size += 3;
  }

  Vector output(size);

  int curIndex = 0;
  if((flags & Flags::position))
  {
    output.segment<3>(curIndex) = position();
    curIndex += 3;
  }
  if((flags & Flags::orientation))
  {
    output.segment<4>(curIndex) = orientation.toVector4();
    curIndex += 4;
  }
  if((flags & Flags::linVel))
  {
    output.segment<3>(curIndex) = linVel();
    curIndex += 3;
  }
  if((flags & Flags::angVel))
  {
    output.segment<3>(curIndex) = angVel();
    curIndex += 3;
  }
  if((flags & Flags::linAcc))
  {
    output.segment<3>(curIndex) = linAcc();
    curIndex += 3;
  }
  if((flags & Flags::angAcc))
  {
    output.segment<3>(curIndex) = angAcc();
  }

  return output;
}

inline Vector Kinematics::toVector() const
{
  int size = 0;
  if(position.isSet())
  {
    size += 3;
  }
  if(orientation.isSet())
  {
    size += 4;
  }
  if(linVel.isSet())
  {
    size += 3;
  }
  if(angVel.isSet())
  {
    size += 3;
  }
  if(linAcc.isSet())
  {
    size += 3;
  }
  if(angAcc.isSet())
  {
    size += 3;
  }

  Vector output(size);

  int curIndex = 0;
  if(position.isSet())
  {
    output.segment<3>(curIndex) = position();
    curIndex += 3;
  }
  if(orientation.isSet())
  {
    output.segment<4>(curIndex) = orientation.toQuaternion().coeffs();
    curIndex += 4;
  }
  if(linVel.isSet())
  {
    output.segment<3>(curIndex) = linVel();
    curIndex += 3;
  }
  if(angVel.isSet())
  {
    output.segment<3>(curIndex) = angVel();
    curIndex += 3;
  }
  if(linAcc.isSet())
  {
    output.segment<3>(curIndex) = linAcc();
    curIndex += 3;
  }
  if(angAcc.isSet())
  {
    output.segment<3>(curIndex) = angAcc();
  }

  return output;
}

inline void Kinematics::reset()
{
  position.reset();
  orientation.reset();
  linVel.reset();
  angVel.reset();
  linAcc.reset();
  angAcc.reset();
}

inline const Kinematics & Kinematics::update_deprecated(const Kinematics & newValue, double dt, Flags::Byte flags)
{
  bool flagPos = flags & Flags::position;
  bool flagLinVel = flags & Flags::linVel;
  bool flagLinAcc = flags & Flags::linAcc;
  bool flagOri = flags & Flags::orientation;
  bool flagAngVel = flags & Flags::angVel;
  bool flagAngAcc = flags & Flags::angAcc;

  {
    CheckedVector3 curPos, curLinVel;

    if(flagPos)
    {
      if(flagLinVel && !newValue.linVel.isSet() && position.isSet() && newValue.position.isSet())
      {
        curPos = position;
      }
      if(newValue.position.isSet())
      {
        position = newValue.position;
      }
      else
      {
        if(position.isSet())
        {
          if(linVel.isSet())
          {
            position() += linVel() * dt;
            if(linAcc.isSet())
            {
              position() += linAcc() * dt * dt / 2;
            }
          }
        }
        else
        {
          position.set().setZero();
        }
      }
    }

    if(flagLinVel)
    {
      if(flagLinAcc && !newValue.linAcc.isSet() && newValue.linVel.isSet())
      {
        curLinVel = linVel;
      }

      if(newValue.linVel.isSet())
      {
        linVel = newValue.linVel;
      }
      else
      {
        if(newValue.position.isSet() && (curPos.isSet() || (position.isSet() && !flagPos)))
        {
          if(curPos.isSet())
          {
            linVel = (newValue.position() - curPos()) / dt;
          }
          else
          {
            linVel = (newValue.position() - position()) / dt;
          }
        }
        else
        {
          if(linVel.isSet())
          {
            if(linAcc.isSet())
            {
              linVel() += linAcc() * dt;
            }
          }
          else
          {
            linVel.set().setZero();
          }
        }
      }
    }

    if(flagLinAcc)
    {
      if(newValue.linAcc.isSet())
      {
        linAcc = newValue.linAcc;
      }
      else
      {
        if(newValue.linVel.isSet() && (curLinVel.isSet() || (linVel.isSet() && !flagLinVel)))
        {
          if(curLinVel.isSet())
          {
            linAcc = (newValue.linVel() - curLinVel()) / dt;
          }
          else
          {
            linAcc = (newValue.linVel() - linVel()) / dt;
          }
        }
        else
        {
          if(!linAcc.isSet())
          {
            linAcc.set().setZero();
          }
        }
      }
    }
  }

  {
    Orientation curOri;
    CheckedVector3 curAngVel;

    if(flagOri)
    {
      if(flagAngVel && !newValue.angVel.isSet() && orientation.isSet() && newValue.orientation.isSet())
      {
        curOri = orientation;
      }
      if(newValue.orientation.isSet())
      {
        orientation = newValue.orientation;
      }
      else
      {
        if(orientation.isSet())
        {
          if(angVel.isSet())
          {
            Vector3 increment = Vector3::Zero();
            increment += angVel() * dt;
            if(angAcc.isSet())
            {
              increment += angAcc() * dt * dt / 2;
            }
            orientation.integrate(increment);
          }
        }
        else
        {
          orientation = Quaternion::Identity();
        }
      }
    }

    if(flagAngVel)
    {
      if(flagAngAcc && !newValue.angAcc.isSet() && newValue.angVel.isSet())
      {
        curAngVel = angVel;
      }

      if(newValue.angVel.isSet())
      {
        angVel = newValue.angVel;
      }
      else
      {
        if(newValue.orientation.isSet() && (curOri.isSet() || (orientation.isSet() && !flagOri)))
        {
          if(curOri.isSet())
          {
            angVel = curOri.differentiate(newValue.orientation) / dt;
          }
          else
          {
            angVel = orientation.differentiate(newValue.orientation) / dt;
          }
        }
        else
        {
          if(angVel.isSet())
          {
            if(angAcc.isSet())
            {
              angVel() += angAcc() * dt;
            }
          }
          else
          {
            angVel.set().setZero();
          }
        }
      }
    }

    if(flagAngAcc)
    {
      if(newValue.angAcc.isSet())
      {
        angAcc = newValue.angAcc;
      }
      else
      {
        if(newValue.angVel.isSet() && (curAngVel.isSet() || (angVel.isSet() && !flagAngVel)))
        {
          if(curAngVel.isSet())
          {
            angAcc = (newValue.angVel() - curAngVel()) / dt;
          }
          else
          {
            angAcc = (newValue.angVel() - angVel()) / dt;
          }
        }
        else
        {
          if(!angAcc.isSet())
          {
            angAcc.set().setZero();
          }
        }
      }
    }
  }

  return *this;
}

} // namespace kine

} // namespace stateObservation

inline std::ostream & operator<<(std::ostream & os, const stateObservation::kine::Kinematics & k)
{
  if(!k.position.isSet() && !k.orientation.isSet() && !k.linVel.isSet() && !k.angVel.isSet() && !k.linAcc.isSet()
     && !k.angAcc.isSet())
  {
    os << "empty kinematics" << std::endl;
  }

  if(k.position.isSet() || k.orientation.isSet())
  {
    if(k.position.isSet())
    {
      os << "pos   : " << k.position().transpose() << "   ";
    }
    else
    {
      os << "pos   :                                  ";
    }

    if(k.orientation.isSet())
    {
      os << "ori   : " << k.orientation.toRotationVector().transpose() << std::endl;
    }
    else
    {
      os << "ori   :" << std::endl;
    }
  }

  if(k.linVel.isSet() || k.angVel.isSet())
  {
    if(k.linVel.isSet())
    {
      os << "linVel: " << k.linVel().transpose() << "   ";
    }
    else
    {
      os << "linVel:                                  ";
    }

    if(k.angVel.isSet())
    {
      os << "angVel: " << k.angVel().transpose() << std::endl;
    }
    else
    {
      os << "angVel:" << std::endl;
    }
  }

  if(k.linAcc.isSet() || k.angAcc.isSet())
  {
    if(k.linAcc.isSet())
    {
      os << "linAcc: " << k.linAcc().transpose() << "   ";
    }
    else
    {
      os << "linAcc:                                  ";
    }

    if(k.angAcc.isSet())
    {
      os << "angAcc: " << k.angAcc().transpose() << std::endl;
    }
    else
    {
      os << "angAcc:" << std::endl;
    }
  }

  return os;
}
