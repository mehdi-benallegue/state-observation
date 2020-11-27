/**
 * \file      miscellaneous-algorithms.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief      Gathers many kinds of algorithms
 *
 *
 *
 */

#ifndef STATEOBSERVATIONTOOLSMISCELANEOUSALGORITHMS
#define STATEOBSERVATIONTOOLSMISCELANEOUSALGORITHMS
#include <boost/utility.hpp>
#include <cmath>

#include <state-observation/api.h>
#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{
namespace tools
{

/// computes the square of a value of any type
template<class T>
inline T square(const T & x)
{
  return T(x * x);
}

/// derivates any type with finite differences
template<class T>
inline T derivate(const T & o1, const T & o2, double dt)
{
  T o(o2 - o1);
  return o * (1 / dt);
}

/// gives the sign of a variable (1, 0 or -1)
template<typename T>
inline int signum(T x)
{
  return (T(0) < x) - (x < T(0));
}

template<typename T>
inline std::string toString(T val)
{
  std::stringstream ss("");
  ss << val;
  return ss.str();
}

/// @brief checks if the vector is already normalized or not
///
/// @param v                  the vector to normalize
/// @return true              the vector is normalized
/// @return false             The vector is not normalized
inline bool checkIfNormalized(const Vector3 & v)
{
  if(fabs(v.squaredNorm() - 1) > cst::epsilon1)
  {
    return false;
  }
  else
  {
    return true;
  }
}

/// @brief checks if the vector is already normalized or not
///
/// @param v                  the vector to normalize
/// @param outputSquaredNorm  the squared norm as an output
/// @return true              the vector is normalized
/// @return false             The vector is not normalized
inline bool checkIfNormalized(const Vector3 & v, double & outputSquaredNorm)
{
  outputSquaredNorm = v.squaredNorm();
  if(fabs(outputSquaredNorm - 1) > cst::epsilon1)
  {
    return false;
  }
  else
  {
    return true;
  }
}

/// @brief returns the value clamped between min and max
///
/// @tparam T The type of the scalar
/// @param x the input value
/// @param max the max value
/// @param min the min value
/// @return T the calmped value
template<typename T>
inline T clampScalar(const T & x, const T & max, const T & min)
{
  if(x > max)
  {
    return max;
  }
  else
  {
    if(x < min)
    {
      return min;
    }
    else
    {
      return x;
    }
  }
}

/// @brief returns the value clamped between limit and -limit
///
/// @tparam T The type of the scalar
/// @param x the input value
/// @param limit is the maximum absolute value allowed, has to be positive
/// @return T the clamped value
template<typename T>
inline T clampScalar(const T & x, const T & limit)
{
  return clampScalar(x, limit, -limit);
}

/// @brief normalize the vector only if it is not normalized already. Useful if the vector is likely to be normalized
///
/// @param v the input vector
/// @return Vector3
inline Vector3 normalizedLazy(const Vector3 & v)
{
  double squaredNorm;
  if(checkIfNormalized(v, squaredNorm))
  {
    return v;
  }
  else
  {
    return v / sqrt(squaredNorm);
  }
}

/// provides an acceleration giving a finite time convergence to zero
/// the state is the position x and the derivative xd and the output is the
/// acceleration. The gains kp, kv must be negative
inline double STATE_OBSERVATION_DLLAPI finiteTimeAccControl(double x, double xd, double kp = -1, double kv = -1)
{
  double sax = sqrt(fabs(x));
  double xdr = kp * signum(x) * sax;
  double y = xd - xdr;
  double ydr = -kv * signum(y) * sqrt(fabs(y));
  return ydr - kp * xd / (2 * sax);
}

/// sqme as the scalar version but for every member of the vector
inline Vector STATE_OBSERVATION_DLLAPI finiteTimeAccControl(const Vector & x,
                                                            const Vector & xd,
                                                            double kp = -1,
                                                            double kv = -1)
{
  Vector xdd(x.size());
  for(Index i = 1; i < x.size(); ++i)
  {
    xdd(i) = finiteTimeAccControl(x(i), xd(i), kp, kv);
  }
  return xdd;
}

namespace Detail
{
double constexpr sqrtNewtonRaphson(double x, double curr, double prev)
{
  return curr == prev ? curr : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
}
} // namespace Detail

/// @brief Constexpr version of the square root
/// @details For a finite and non-negative value of "x", returns an approximation for the square root of "x"
///-Otherwise, returns NaN
///
/// @param x
/// @return double constexpr
double constexpr sqrt(double x)
{
  return x >= 0 && x < std::numeric_limits<double>::infinity() ? Detail::sqrtNewtonRaphson(x, x, 0)
                                                               : std::numeric_limits<double>::quiet_NaN();
}

} // namespace tools

} // namespace stateObservation

#endif // STATEOBSERVATIONTOOLSMISCELANEOUSALGORITHMS
