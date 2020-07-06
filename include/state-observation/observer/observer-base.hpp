/**
 * \file      observer-base.hpp
 * \author    Mehdi Benallegue
 * \date       2012
 * \brief      Defines the base class of a state
 *             observer.
 *             The observer is destinated to any dynamical system with a vector
 *             state representation
 *             The file describes also
 *             the used data structures and exceptions of the derivated classes
 *
 * \details
 *
 *
 */



#ifndef OBSERVERBASEHPP
#define OBSERVERBASEHPP



#include <boost/static_assert.hpp>
#include <boost/assert.hpp>

#include <state-observation/api.h>
#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{
     /**
     * \class  ObserverBase
     * \brief  The base class for observers.
     *         The observer is destinated to any dynamical system with a vector
     *         state representation. This class mostly defined an abstract
     *         interface, static constants and types. It is templated by:
     *
     *
     */


    class STATE_OBSERVATION_DLLAPI ObserverBase
    {
    public:

        ///StateVector is the type of state vector
        typedef Vector StateVector;

        ///MeasureVector is the type of measurements vector
        typedef Vector MeasureVector;

        ///InputVector is the type of the input vector
        typedef Vector InputVector;


        ///constructor
        ///      \li n : size of the state vector
        ///      \li m : size of the measurements vector
        ///      \li p : size of the input vector
        ObserverBase(unsigned n, unsigned m, unsigned p=0);

        ///default constructor (default values for n,m,p are zero)
        ObserverBase();

        ///Destructor
        virtual ~ObserverBase(){};

        /// Changes the size of the state vector
        virtual void setStateSize(unsigned n);

        /// gets the size of the state vector
        virtual unsigned getStateSize() const;

        /// Changes the size of the measurement vector
        virtual void setMeasureSize(unsigned m);

        /// gets the size of the measurement vector
        virtual unsigned getMeasureSize() const;

        /// Changes the size of the input vector
        virtual void setInputSize(unsigned p);

        /// gets the size of the input vector
        virtual unsigned getInputSize() const;

        ///Set the value of the state vector at time index k
        virtual void setState(const StateVector& x_k,TimeIndex k)=0;

        ///Remove all the given past values of the state
        virtual void clearStates()=0;

        ///Set the value of the measurements vector at time index k
        virtual void setMeasurement(const MeasureVector& y_k,TimeIndex k)=0;

        ///Remove all the given past values of the measurements
        virtual void clearMeasurements()=0;

        ///Set the value of the input vector at time index k
        virtual void setInput(const InputVector& x_k,TimeIndex k)=0;

        ///Remove all the given past values of the inputs
        virtual void clearInputs()=0;

        ///Run the observer loop and gets the state estimation of the state at
        ///instant k
        virtual StateVector getEstimatedState(TimeIndex k)=0;

        ///Reinitializes the whole observer
        ///default behavior is to call the three "ObserverBase::clear*" methods
        virtual void reset();

        ///These are initializer  for vectors
        /// We do not use the get prefix to be consistent with eigen initializers
        ///Gives a vector of state vector size having duplicated "c" value
        virtual StateVector stateVectorConstant( double c ) const;

        ///Gives a vector of state vector size having random values
        virtual StateVector stateVectorRandom() const;

        ///Gives a vector of state vector size having zero values
        virtual StateVector stateVectorZero() const;

        ///Tells whether or not the vector has the dimensions of a state vector
        virtual bool checkStateVector(const StateVector & v ) const;

        ///Gives a vector of measurement vector size having duplicated "c" value
        virtual MeasureVector measureVectorConstant( double c ) const;

        ///Gives a vector of measurement vector size having random values
        virtual MeasureVector measureVectorRandom() const;

        ///Gives a vector of measurement vector size having zero values
        virtual MeasureVector measureVectorZero() const;

        ///Tells whether or not the vector has the dimensions of a measurement vector
        virtual bool checkMeasureVector(const MeasureVector &) const;

        ///Gives a vector of input vector size having duplicated "c" value
        virtual InputVector inputVectorConstant( double c ) const;

        ///Gives a vector of input vector size having random values
        virtual InputVector inputVectorRandom() const;

        ///Gives a vector of input vector size having zero values
        virtual InputVector inputVectorZero() const;

        ///Tells whether or not the vector has the dimensions of a input vector
        virtual bool checkInputVector(const InputVector &) const;

    protected:

        ///stateSize is the size of the state vector
        unsigned n_;

        ///measureSize is the size of measurements vector
        unsigned m_;

        ///inputSize is the size of the input vector
        unsigned p_;

    };

}

#endif//OBSERVERBASEHPP
