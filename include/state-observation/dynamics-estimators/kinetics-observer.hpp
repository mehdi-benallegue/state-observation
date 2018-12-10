/**
 * \file      kinetics-observer-base.hpp
 * \author    Mehdi Benallegue
 * \date      2018
 * \brief     Unified Kinetics estimator
 *
 * \details
 *
 *
 */

#ifndef FLEXIBILITYESTIMATION_EKFFLEXIBILITYESTIMATORBASE_H
#define FLEXIBILITYESTIMATION_EKFFLEXIBILITYESTIMATORBASE_H

#include <boost/utility.hpp>

#include <state-observation/observer/extended-kalman-filter.hpp>

#include <state-observation/flexibility-estimation/flexibility-estimator-base.hpp>


namespace stateObservation
{
namespace  kineticsObserver
{
   /**
    * \class  EKFFlexibilityEstimatorBase
    * \brief  This class is the base class of the flexibility estimators that
    *         use an extended Kalman Filter. Several methods require to be overloaded
    *         to derive an implementation from this base class.
    *
    */

    class KineticsObserver
    {
    public:
        ///virtual destructor
        virtual ~KineticsObserver();

        /// The constructor.
        ///  \li maxContacts : maximum number of contacts,
        ///  \li dx gives the derivation step for a finite differences derivation method
        KineticsObserver(int maxContacts=4, const Vector & dx=Vector::Zero(0));




        /// ///////////////////////////////////////////////////////////
        /// Getting, setting the current time and running the estimation
        /// //////////////////////////////////////////////////////////

        /// gets the time
        double getTime() const;

        /// sets the time
        void setTime(double) ;

        /// this function triggers the estimation itself
        void update();




        /// ////////////////////////////////
        /// Getting and setting the state
        /// ////////////////////////////////

        /// Gets an estimation of the state in the form of a state vector $\hat{x_{k+1}}$
        virtual Vector getState() const;

        /// Get the kinematics of a given frame
        /// the input is the state of the desired frame  expressed in local cordinates
        /// the output is the state of the frame expressed in global cordinates
        /// the type of the vector and the output depends on its size
        /// if size == 0 then the kinematic state is provided up to 2nd order
        ///              (position, velocity and acceleration)
        /// if size == 3 it is a position and the output will be the estimated global pos
        /// if size == 6 then it is a position/orientation and
        ///              the output will be the global position/orientation
        /// if size == 9 then it is a position/orientation/linVelocity
        ///               the output will be the global value accordingly
        /// if size == 12 then it is a position/orientation/linVelocity/angVelocity
        ///               the output will be the global value accordingly
        /// if size == 15 if the acceleration is estimated then it is a
        ///                 position/orientation/linVelocity/angVelocity/linacc
        ///                 the output will be the global value accordingly
        ///               otherwise
        ///                 it acts as size == 12
        /// if size == 15 if the acceleration is estimated then it is a
        ///                 position/orientation/linVelocity/angVelocity/linacc/angacc
        ///                 the output will be the global value accordingly
        ///               otherwise
        ///                 it acts as size == 12
        virtual Vector getKinematics(const Vector &) const;


        virtual Vector getExternalForces() const;

        ///This function allows to estimate the acceleration even if
        /// withAccelerationEstimation is not set
        virtual void estimateAcceleration();

        /// Gets the state size
        virtual unsigned getStateSize() const;


        ///Sets a value of the state x_k provided from another source
        /// can be used for initialization of the estimator
        ///can be a 4x4 homogenous matrix or a vector of size stateSize
        virtual void setState(const Matrix &, bool resetCovariance=true);


        ///Sets a value for the kinematics part of the state
        virtual void setKinematics(const Matrix &, bool resetForces=true,
                                   bool resetCovariance=true);



        /// /////////////////////////////////////////////////////
        /// Setting and getting the state of the estimation
        /// ////////////////////////////////////////////////////
        void setWithExternalForces(bool b = true);

        void setWithAccelerationEstimation(bool b = true);

        void setWithFilteringMeasurements(bool b = true);



        /// ///////////////////////////////////////////////
        /// Getting and setting input data and measurements
        /// /////////////////////////////////////////////

        /// sets the measurement of the IMU (gyrometer, accelerometer,)
        virtual void setIMU(const Vector3 & accelero, const  Vector3 & gyrometer, const Vector3 & position,
                            const Matrix & orientation, const  Vector3 & linearVel= Vector::Zero(0),
                                      Vector3 angularVel= Vector::Zero(0),int num=0);

        /// Gets an estimation of the state in the form of a state vector $\hat{x_{k+1}}$
        /// if the IMU contains both
        virtual void setGyrometer(Vector3 measurement, Vector3 position,
                                      Matrix orientation, Vector3 linearVel= Vector::Zero(0),
                                      Vector3 angularVel= Vector::Zero(0),int num=0);

        /// Gets an estimation of the state in the form of a state vector $\hat{x_{k+1}}$
        virtual void setAccelerometer(Vector3 measurement, Vector3 position,
                                      Matrix orientation, Vector3 linearVel= Vector::Zero(0),
                                      Vector3 angularVel= Vector::Zero(0),
                                      Vector3 linearAcc=Vector::Zero(0),int num=0);

        ///this function allows to compute the filtered sensors even if
        /// WithFilteringMeasurements is not set
        virtual void filterMeasurements();



        ///Sets the covariance matrix of the flexibility Guess
        virtual void setStateCovariance(const Matrix & P);

        ///Gets the covariance matrix of the flexibility
        virtual Matrix getStateCovariance() const;

        ///Sets the covariance matrices for the process noises
        /// \li Q process noise
        virtual void setProcessNoiseCovariance(const Matrix & Q);

        ///Sets the covariance matrices for the sensor noises
        /// \li R sensor noise
        virtual void setMeasurementNoiseCovariance(const Matrix & R);

        ///gets the covariance matrices for the process noises
        virtual Matrix getProcessNoiseCovariance() const ;

        ///gets the covariance matrices for the sensor noises
        virtual Matrix getMeasurementNoiseCovariance() const ;

        /// Sets the value of the next sensor measurement y_{k+1}
        virtual void setMeasurement(const Vector & y);

        virtual Vector getMeasurement();

        /// Sets the value of the next input for the state process dynamics
        /// i.e. : gives u_k such that x_{k+1} = f(x_k,u_k)
        virtual void setInput(const Vector & u);

        /// Sets the value of the next  measurement
        /// i.e. : gives u_{k+1} such that y_{k+1}=h(x_{k+1},u_{k+1})
        virtual void setMeasurementInput(const Vector & u);

        virtual Vector getInput();

        virtual Vector getMeasurementInput();



        /// Gets a const reference on the extended Kalman filter
        virtual const stateObservation::ExtendedKalmanFilter & getEKF() const;

        /// Gets a reference on the extended Kalman filter
        virtual stateObservation::ExtendedKalmanFilter & getEKF();



        /// Gets the measurements size
        /// this method is pure virtual and reauires to be overloaded in implementation
        virtual unsigned getMeasurementSize() const =0;

        /// Gets the input size
        /// this method is pure virtual and reauires to be overloaded in implementation
        virtual unsigned getInputSize() const =0;

        /// Gets a simulation of the
        virtual Vector getSimulatedMeasurement();

        ///Resets the covariance matrices to their original values
        virtual void resetCovarianceMatrices()=0;

        ///Get the last vector of inovation of the Kalman filter
        virtual Vector getInovation();

        ///Get the simulated measurement of the predicted state
        virtual Vector getPredictedMeasurement();

        ///Get the predicted state
        virtual Vector getPrediction();

        ///Get the last simulated measurement
        virtual Vector getLastPredictedMeasurement();

        ///Get the last predicted state
        virtual Vector getLastPrediction();

    protected:
        virtual void setJacobians(const Matrix & A, const Matrix & C);

        virtual void useFiniteDifferencesJacobians(Vector dx);

        stateObservation::ExtendedKalmanFilter ekf_;

        bool finiteDifferencesJacobians_;

        Vector dx_;

        Vector lastX_;

        TimeIndex k_;

    private:

    };
}
}
#endif // FLEXIBILITYESTIMATION_EKFFLEXIBILITYESTIMATORBASE_H
