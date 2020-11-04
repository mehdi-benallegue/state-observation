/**
 * \file      fixed-contact-ekf-flex-estimator-imu.hpp
 * \author    Mehdi Benallegue
 * \date      2013
 * \brief     Declares the class of the estimation of the flexibility using an
 *            extended Kalman filter and a fixed contact hypothesis
 *
 * \details
 *
 *
 */


#ifndef FLEXBILITYESTMATOR_FIXEDCONTACTEKFFLEXIBILITYESTIMATOR_IMU_H
#define FLEXBILITYESTMATOR_FIXEDCONTACTEKFFLEXIBILITYESTIMATOR_IMU_H

#include <state-observation/api.h>
#include <state-observation/flexibility-estimation/ekf-flexibility-estimator-base.hpp>
#include <state-observation/flexibility-estimation/imu-fixed-contact-dynamical-system.hpp>

namespace stateObservation
{
namespace flexibilityEstimation
{

    /**
    * \class  FixedContactEKFFlexEstimatorIMU
    * \brief  This class implements the flexibility estimation of a robot with
    *         the hypothesis that the contact positions do not move. This constraint
    *         is expressed using fictious measurements but the interface is transparent
    *         to this assumption, the state is expressed using classical representation
    *         of position, velocity, acceleration, orientation (using (theta x mu) representation)
    *         angular velocity (omega) and acceleration (omega dot)
    *
    */

    class STATE_OBSERVATION_DLLAPI FixedContactEKFFlexEstimatorIMU :
                                            public EKFFlexibilityEstimatorBase,
                                            private boost::noncopyable
    {
    public:

        ///The constructor, it requires the value of the time discretization period
        explicit FixedContactEKFFlexEstimatorIMU( double dt=0.005 );

        ///Virtual destructor
        virtual ~FixedContactEKFFlexEstimatorIMU();

        ///Sets the number of contacts can be changed online
        void setContactsNumber(unsigned i);

        ///Sets the position of the i-th contact
        void setContactPosition(unsigned i, Vector3 position);

        /// Sets the value of the next sensor measurement y_{k+1}
        virtual void setMeasurement(const Vector & y);

        ///Sets the covariance of the fictious measurements (not mandatory)
        virtual void setVirtualMeasurementsCovariance(double c_);

        ///Sets the covariance of the fictious measurements (not mandatory)
        virtual double getVirtualMeasurementsCovariance() const;

        ///Sets the process covariance matrice
        virtual void setProcessNoiseCovariance(const Matrix & Q);

        ///Sets the measurements covariance matrice
        virtual void setMeasurementNoiseCovariance(const Matrix & R);

        ///gets the covariance matrices for the process noises
        virtual Matrix getProcessNoiseCovariance() const ;

        ///gets the covariance matrices for the sensor noises
        virtual Matrix getMeasurementNoiseCovariance() const ;

        ///Sets a value of the flexibility x_k provided from another source
        /// can be used for initialization of the estimator
        virtual void setFlexibilityGuess(const Matrix & x);

        /// Gets an estimation of the flexibility in the form of a homogeneous matrix
        virtual Matrix4 getFlexibility();

        /// Gets an estimation of the flexibility in the form of a state vector \hat{x_{k+1}}
        virtual const Vector& getFlexibilityVector();


        virtual Index getMeasurementSize() const ;

        virtual Index getStateSize() const ;

        virtual Index getInputSize() const ;

        /// sets the sampling period
        virtual void setSamplingPeriod(double);

        ///Resets the covariance matrices to their original values
        virtual void resetCovarianceMatrices();


    protected:

        typedef kine::indexes<kine::rotationVector> indexes;

        virtual void updateCovarianceMatrix_();

        IMUFixedContactDynamicalSystem functor_;

        double virtualMeasurementCovariance_;

        Matrix R_,Q_;

        static const Index stateSizeConst_=18;
        static const Index measurementSizeConst_=6;
        static const Index inputSizeConst_=15;

        double dt_;//sampling period

    private:
    };

}
}
#endif // FLEXBILITYESTMATOR_FIXEDCONTACTEKFFLEXIBILITYESTIMATOR_H
