/**
 * \file      algebraic-sensor.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief     Gives a base class for algebraic sensors
 *
 *
 */



#ifndef SIMULATIONALGEBRAICSENSORHPP
#define SIMULATIONALGEBRAICSENSORHPP

#include <Eigen/Core>
#include <boost/assert.hpp>

#include <state-observation/api.h>
#include <state-observation/sensors-simulation/sensor-base.hpp>

namespace stateObservation
{
    /**
     * \class  AlgebraicSensor
     * \brief  The base class for algebraic sensors. Algebraic sensors are sensors
     *         which depend only on the state value and the current time
     *          and do not have internal dynamics
     *         (or a dynamics which converges fast enough to be ignored). This class
     *         implements mostly the containers and the interface to algebraic sensors.
     *         Algebraic sensors must be derived from this class.
     *
     * \details
     *
     */

    class STATE_OBSERVATION_DLLAPI AlgebraicSensor: public SensorBase
    {
    public:
        /// Default constructor
        AlgebraicSensor();

        ///virtual destructor
        virtual ~AlgebraicSensor(){}

        ///gets the measurement of the current time. We can choose to consider
        ///noise or not (default is noisy)
        virtual Vector getMeasurements(bool noisy=true);

        ///Sets the value of the state at instant k
        virtual void setState(const Vector & state, TimeIndex k);

        ///gets the current time
        virtual TimeIndex getTime() const;

        ///gets the state vector size. Pure virtual method.
        virtual Index getStateSize() const;

        ///get the size of the measurements. Pure virtual method.
        virtual Index getMeasurementSize() const;

        ///concatenates the n last components of the state in the measurement
        ///(useful when the measurements are already computed or
        ///when they come from external source)
        virtual Index concatenateWithInput( Index n);

    protected:
        ///the actual algorithm for the computation of the measurements, must
        ///be overloaded to implement any sensor
        virtual Vector computeNoiselessMeasurement_()=0;

        virtual Index getStateSize_() const=0;

        virtual Index getMeasurementSize_() const=0;


        Vector computeNoisyMeasurement_();

        virtual void checkState_(const Vector &);

        TimeIndex time_;

        Index concat_;

        Vector state_;

        Vector directInputToOutput_;

        bool storedNoisyMeasurement_;

        Vector noisyMeasurement_;

        bool storedNoiselessMeasurement_;

        Vector noiselessMeasurement_;

    };

}

#endif //SIMULATIONALGEBRAICSENSORHPP
