#include <state-observation/flexibility-estimation/fixed-contact-ekf-flex-estimator-imu.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>

const double initialVirtualMeasurementCovariance=1.e-10;

const double dxFactor = 1.0e-8;

namespace stateObservation
{
namespace flexibilityEstimation
{
    FixedContactEKFFlexEstimatorIMU::FixedContactEKFFlexEstimatorIMU(double dt):
        EKFFlexibilityEstimatorBase
            (stateSizeConst_,measurementSizeConst_,inputSizeConst_,
                            Matrix::Constant(getStateSize(),1,dxFactor)),
        functor_(dt),
        virtualMeasurementCovariance_(initialVirtualMeasurementCovariance)
    {
        ekf_.setDirectInputStateFeedthrough(false);

        ekf_.setMeasureSize(functor_.getMeasurementSize());

        FixedContactEKFFlexEstimatorIMU::resetCovarianceMatrices();

        Vector dx = Matrix::Constant(getStateSize(),1,dxFactor);

        dx.segment(indexes::ori,3) = Vector3::Constant(1e-4) ;
        dx.segment(indexes::angVel,3) = Vector3::Constant(1e-4) ;

        FixedContactEKFFlexEstimatorIMU::useFiniteDifferencesJacobians(dx);

        Vector x0=ekf_.stateVectorZero();

        lastX_=x0;

        ekf_.setState(x0,0);

        ekf_.setStateCovariance(Q_);

        ekf_.setFunctor(& functor_);
    }



    FixedContactEKFFlexEstimatorIMU::~FixedContactEKFFlexEstimatorIMU()
    {
        //dtor
    }

    void FixedContactEKFFlexEstimatorIMU::resetCovarianceMatrices()
    {
        R_=Matrix::Identity(getMeasurementSize(),getMeasurementSize());
        R_.block(0,0,3,3)=Matrix3::Identity()*1.e-6;//accelerometer
        R_.block(3,3,3,3)=Matrix3::Identity()*1.e-6;//gyrometer

        updateCovarianceMatrix_();

        Q_=ekf_.getQmatrixIdentity();
        Q_=Q_*1.e-8;
        Q_.block(indexes::linVel,indexes::linVel,3,3)=Matrix3::Identity()*1.e-4;
        Q_.block(indexes::angVel,indexes::angVel,3,3)=Matrix3::Identity()*1.e-4;
        Q_.block(indexes::linAcc,indexes::linAcc,3,3)=Matrix3::Identity()*1.e-2;
        Q_.block(indexes::angAcc,indexes::angAcc,3,3)=Matrix3::Identity()*1.e-2;

        ekf_.setQ(Q_);

        Matrix P0 (ekf_.getQmatrixIdentity());
        P0=P0*1e-2;
        P0.block(indexes::linVel,indexes::linVel,3,3)=Matrix3::Identity()*1.e-2;
        P0.block(indexes::angVel,indexes::angVel,3,3)=Matrix3::Identity()*1.e-2;
        P0.block(indexes::linAcc,indexes::linAcc,3,3)=Matrix3::Identity()*1.e-2;
        P0.block(indexes::angAcc,indexes::angAcc,3,3)=Matrix3::Identity()*1.e-2;

        ekf_.setStateCovariance(P0);

    }

    void FixedContactEKFFlexEstimatorIMU::setContactsNumber(unsigned i)
    {
        functor_.setContactsNumber(i);
        ekf_.setMeasureSize(functor_.getMeasurementSize());
        updateCovarianceMatrix_();
    }

    void FixedContactEKFFlexEstimatorIMU::setContactPosition
                                            (unsigned i, Vector3 position)
    {
        functor_.setContactPosition(i,position);
    }

    void FixedContactEKFFlexEstimatorIMU::setMeasurement(const Vector & y)
    {
        BOOST_ASSERT((getMeasurementSize()==y.size()) &&
                "ERROR: The measurement vector has incorrect size");


        Vector y2 = ekf_.measureVectorZero();
        y2.head(getMeasurementSize()) = y;
        ekf_.setMeasurement(y2,k_+1);
    }

    void FixedContactEKFFlexEstimatorIMU::setVirtualMeasurementsCovariance
                                                                    (double c)
    {
        virtualMeasurementCovariance_=c;
        updateCovarianceMatrix_();
    }

    double FixedContactEKFFlexEstimatorIMU::getVirtualMeasurementsCovariance() const
    {
        return virtualMeasurementCovariance_;
    }

    void FixedContactEKFFlexEstimatorIMU::setFlexibilityGuess(const Matrix & x)
    {
        bool bstate =ekf_.checkStateVector(x);
        bool b6= (x.rows()==6 && x.cols()==1);
        bool bhomogeneous = (x.rows()==4 && x.cols()==4);

        (void)b6;//avoid warning

        BOOST_ASSERT((bstate||b6||bhomogeneous) &&
                "ERROR: The flexibility state has incorrect size \
                    must be 18x1 vector, 6x1 vector or 4x4 matrix");

        Vector x0=x;

        if (bstate)
        {
            ekf_.setState(x0,k_);
        }
        else
        {
            if (bhomogeneous)
                x0=kine::homogeneousMatrixToVector6(x);

            Vector x_s = ekf_.stateVectorZero();

            x_s.segment(indexes::pos,3)=x0.head(3);

            x_s.segment(indexes::ori,3)=x0.tail(3);

            ekf_.setState(x_s,k_);

            ekf_.setQ(Q_);
        }
    }

    void FixedContactEKFFlexEstimatorIMU::setMeasurementNoiseCovariance
                                            (const Matrix & R)
    {
        BOOST_ASSERT(R.rows()==getMeasurementSize() &&
                     R.cols()==getMeasurementSize() &&
                    "ERROR: The measurement noise covariance matrix R has \
                        incorrect size");

        R_=R;
        updateCovarianceMatrix_();
    }


    void FixedContactEKFFlexEstimatorIMU::setProcessNoiseCovariance
                                            (const Matrix & Q)
    {
        Q_=Q;
        ekf_.setQ(Q_);
    }

    Matrix FixedContactEKFFlexEstimatorIMU::getProcessNoiseCovariance() const
    {
        return Q_;
    }


    Matrix FixedContactEKFFlexEstimatorIMU::getMeasurementNoiseCovariance() const
    {
        return R_;
    }

    void FixedContactEKFFlexEstimatorIMU::updateCovarianceMatrix_()
    {
        Matrix R = ekf_.getRmatrixIdentity();
        R.block(0,0, getMeasurementSize() , getMeasurementSize())=R_;

        for (Index i= getMeasurementSize() ; i<ekf_.getMeasureSize() ; ++i)
        {
            R(i,i)=virtualMeasurementCovariance_;
        }

        ekf_.setR(R);
    }

    Index FixedContactEKFFlexEstimatorIMU::getStateSize() const
    {
        return stateSizeConst_;
    }

    Index FixedContactEKFFlexEstimatorIMU::getInputSize() const
    {
        return inputSizeConst_;
    }

    Index FixedContactEKFFlexEstimatorIMU::getMeasurementSize() const
    {
        return measurementSizeConst_;
    }

    Matrix4 FixedContactEKFFlexEstimatorIMU::getFlexibility()
    {
        Vector v (getFlexibilityVector());
        Vector6 v2;

        v2.head(3) = v.segment(indexes::pos,3);
        v2.tail(3) = v.segment(indexes::ori,3);

        return kine::vector6ToHomogeneousMatrix(v2);
    }

    const Vector& FixedContactEKFFlexEstimatorIMU::getFlexibilityVector()
    {
        if (ekf_.getMeasurementsNumber()>0)
        {
            lastX_ =EKFFlexibilityEstimatorBase::getFlexibilityVector();

            ///regulate the part of orientation vector in the state vector
            lastX_.segment(indexes::ori,3)=
                kine::regulateOrientationVector(lastX_.segment(indexes::ori,3));

            ekf_.setState(lastX_,ekf_.getCurrentTime());
        }
        return lastX_;
    }

    void FixedContactEKFFlexEstimatorIMU::setSamplingPeriod(double dt)
    {
        dt_=dt;
    }


}
}
