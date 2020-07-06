/**
 * \file      accelerometer-gyrometer-magnetometer.hpp
 * \author    Mehdi Benallegue, Joseph Mirabel
 * \date      2015
 * \brief     Implements the accelerometer-gyrometer-magnetometer inertial
 *            measurement unit
 *
 *
 */



#ifndef SIMULATIONACCELEROMETERGYROMETERMAGNETOMETERSENSORHPP
#define SIMULATIONACCELEROMETERGYROMETERMAGNETOMETERSENSORHPP

#include <Eigen/Core>
#include <boost/assert.hpp>

#include <state-observation/api.h>
#include <state-observation/sensors-simulation/algebraic-sensor.hpp>
#include <state-observation/sensors-simulation/algorithm/linear-acceleration.hpp>
#include <state-observation/sensors-simulation/algorithm/rotation-velocity.hpp>
#include <state-observation/sensors-simulation/algorithm/magnetic-field.hpp>

namespace stateObservation
{
    /**
     * \class  AccelerometerGyrometerMagnetometer
     * \brief  Implements the accelerometer-gyrometer-magnetometer measurements
     *
     *
     *
     * \details
     *
     */

    class STATE_OBSERVATION_DLLAPI AccelerometerGyrometerMagnetometer : public AlgebraicSensor,
        protected algorithm::LinearAcceleration,
        protected algorithm::RotationVelocity,
        protected algorithm::MagneticField
    {
    public:
        AccelerometerGyrometerMagnetometer();

        ///Virtual destructor
        virtual ~AccelerometerGyrometerMagnetometer(){}



        void setMatrixMode(bool matrixMode);


    protected:
        ///Gets the state vector Size
        virtual size_t getStateSize_() const;

        ///Gets the measurements vector size
        virtual size_t getMeasurementSize_() const;


        virtual Vector computeNoiselessMeasurement_();

        Matrix3 r_;
        Vector3 acc_;
        Vector3 omega_;
        Vector3 magne_;
        Vector output_;

        bool matrixMode_;

        static const size_t stateSize_= 10;
        static const size_t stateSizeMatrix_= 15;

        static const size_t measurementSize_=9;

        size_t currentStateSize_;

    };

}

#endif // SIMULATIONACCELEROMETERGYROMETERMAGNETOMETERSENSORHPP
