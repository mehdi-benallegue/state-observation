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
        virtual Index getStateSize_() const;

        ///Gets the measurements vector size
        virtual Index getMeasurementSize_() const;


        virtual Vector computeNoiselessMeasurement_();

        Matrix3 r_;
        Vector3 acc_;
        Vector3 omega_;
        Vector3 magne_;
        Vector output_;

        bool matrixMode_;

        static const Index stateSize_= 10;
        static const Index stateSizeMatrix_= 15;

        static const Index measurementSize_=9;

        Index currentStateSize_;

    };

}

#endif // SIMULATIONACCELEROMETERGYROMETERMAGNETOMETERSENSORHPP
