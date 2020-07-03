/**
 * \file      sensor-base.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief      Implements the base class of all sensors.
 *
 *
 */



#ifndef SIMULATIONSENSORBASEHPP
#define SIMULATIONSENSORBASEHPP

#include <Eigen/Core>
#include <boost/assert.hpp>

#include <state-observation/api.h>
#include <state-observation/noise/noise-base.hpp>

namespace stateObservation
{
    /**
     * \class  SensorBase
     * \brief  The base class for sensors. This must be derived to implement a
     *         sensor
     *
     * \details
     *
     */

    class STATE_OBSERVATION_DLLAPI SensorBase
    {
    public:
        ///default constructor
        SensorBase();

        ///virtual destructor
        virtual ~SensorBase(){}

        ///gets the measurement of the current time. We can choose to consider
        ///noise or not (default is noisy)
        virtual Vector getMeasurements(bool noisy=true)=0;

        ///Sets the value of the state at instant k
        virtual void setState(const Vector & state, TimeIndex k)=0;

        ///Sets a pointer on the noise on the measurements. The class does NOT destroy the noise
        ///when the destructor is called.
        virtual void setNoise(NoiseBase *);

        ///gets the pointer on the measurements noise
        virtual NoiseBase* getNoise() const;

        ///removes the noise
        virtual void resetNoise();

        ///gets the current time, pure virtual method
        virtual TimeIndex getTime() const=0;

        ///gets the state vector size. Pure virtual method.
        virtual unsigned getStateSize() const=0;

        ///get the size of the measurements. Pure virtual method.
        virtual unsigned getMeasurementSize() const=0;

        ///gets a zero vector of the size of a state vector
        virtual Vector stateVectorZero() const;

        ///checks whether a vector is correctly sized or not
        virtual bool checkStateVector(const Vector &) const;


    protected:

        NoiseBase *  noise_;

    };



}

#endif //SIMULATIONSENSORBASEHPP
