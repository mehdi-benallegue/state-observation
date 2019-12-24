/**
 * \file      dynamics-functor-base.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief      Gives the base class to be derived in order to define any
 *             dynamics system.
 *
 * \details
 *
 *
 */

#ifndef STATEOBSERVERDYNAMICALSYSTEMFUNCTORBASE_H
#define STATEOBSERVERDYNAMICALSYSTEMFUNCTORBASE_H

#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{

    /**
    * \class  DynamicsFunctorBase
    * \brief
    *        This is the base class of any functor that describes the dynamics
    *        of the state and the measurement.
    *        This class is to be derived in order to be given
    *        to the Extended Kalman Filter.
    *
    */

    class DynamicalSystemFunctorBase
    {
    public:
        DynamicalSystemFunctorBase();
        virtual ~DynamicalSystemFunctorBase();

        ///The function to oberload to describe the dynamics of the state
        virtual Vector stateDynamics
        (const Vector& x, const Vector& u, TimeIndex k)=0;

        ///The function to overload to describe the dynamics of the sensor (measurements)
        virtual Vector measureDynamics
        (const Vector& x, const Vector& u, TimeIndex k)=0;

        ///The method to overload if the functor needs to be reset when the
        ///Exteded Kalman filter is reset itself
        virtual void reset(){}

        ///gets the state size
        virtual unsigned getStateSize()const =0 ;
        ///gets the input size
        virtual unsigned getInputSize()const =0 ;
        ///gets the measurements size
        virtual unsigned getMeasurementSize() const=0;

        ///Gives a boolean answer on whether or not the vector is correctly sized to be a state vector
        virtual bool checkStateVector(const Vector &);
        ///Gives a boolean answer on whether or not the vector is correctly sized to be an input vector
        virtual bool checkInputvector(const Vector &);

    protected:

        inline void assertStateVector_(const Vector & v)
        {
            (void)v;//avoid warning
            BOOST_ASSERT(checkStateVector(v) && "ERROR: The state vector has the wrong size");
        }
        inline void assertInputVector_(const Vector & v)
        {
            (void)v;//avoid warning
            BOOST_ASSERT(checkInputvector(v) && "ERROR: The input vector has the wrong size");
        }

    private:
    };

    /**
     * \class  StateArithmetics
     * \brief
     *        This class is used to customize the way the difference between measurements,
     *        the state update function and the differentiation are performed. 
     *         default is the usual natual arithmetics. overload any ohter one
     *
     */
    class StateArithmetics
    {
    public:
        virtual void stateSum(const  Vector& stateVector, const Vector& tangentVector, Vector& sum);

        virtual void stateDifference(const Vector& stateVector1, const Vector& stateVector2, Vector& difference);

        virtual void measurementDifference(const Vector& measureVector1, const Vector& measureVector2, Vector& difference);
    };
}

#endif // STATEOBSERVERDYNAMICSYSTEMFUNCTORBASE_H
