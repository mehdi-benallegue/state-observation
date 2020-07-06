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

#include <state-observation/api.h>
#include <state-observation/tools/definitions.hpp>


namespace stateObservation
{
    namespace tools
    {

        ///computes the square of a value of any type
        template <class T>
        inline T square (const T & x)
        {
            return T(x*x);
        }

        ///derivates any type with finite differences
        template <class T>
        inline T derivate(const T & o1 , const T & o2 , double dt)
        {
            T o(o2-o1);
            return o*(1/dt);
        }

        ///gives the sign of a variable (1, 0 or -1)
        template <typename T> inline
        int signum(T x)
        {
          return (T(0) < x) - (x < T(0));
        }

        template<typename T>
        std::string toString(T val)
        {
          std::stringstream ss("");
          ss << val;
          return ss.str();
        }



        ///provides an acceleration giving a finite time convergence to zero
        ///the state is the position x and the derivative xd and the output is the
        ///acceleration. The gains kp, kv must be negative
        inline double STATE_OBSERVATION_DLLAPI finiteTimeAccControl(double x, double xd, double kp=-1, double kv=-1)
        {
          double sax = sqrt(fabs(x));
          double xdr = kp*signum(x)*sax;
          double y = xd - xdr;
          double ydr = -kv*signum(y)*sqrt(fabs(y));
          return ydr - kp*xd/(2*sax);
        }


        ///sqme as the scalar version but for every member of the vector
        inline Vector STATE_OBSERVATION_DLLAPI finiteTimeAccControl(const Vector &x, const Vector &xd, double kp=-1, double kv=-1)
        {
          Vector xdd(x.size());
          for (size_t i=1;i<size_t(x.size());++i)
          {
            xdd(i)= finiteTimeAccControl(x(i),xd(i),kp,kv);
          }
          return xdd;

        }
    }


}



#endif //STATEOBSERVATIONTOOLSMISCELANEOUSALGORITHMS
