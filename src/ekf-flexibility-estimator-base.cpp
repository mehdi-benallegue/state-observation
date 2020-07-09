#include <state-observation/flexibility-estimation/ekf-flexibility-estimator-base.hpp>
namespace stateObservation
{
namespace flexibilityEstimation
{
    EKFFlexibilityEstimatorBase::EKFFlexibilityEstimatorBase(Index stateSize,
                                    Index measurementSize,
                                    Index inputSize,
                                    const Vector & dx):
        FlexibilityEstimatorBase(),
        ekf_(stateSize,measurementSize,inputSize),
        finiteDifferencesJacobians_(true),
        dx_(dx),
        k_(0)
    {
        //ctor
    }

    EKFFlexibilityEstimatorBase::~EKFFlexibilityEstimatorBase()
    {
        //dtor
    }

    void EKFFlexibilityEstimatorBase::setFlexibilityCovariance
                                                            (const Matrix & P)
    {
        ekf_.setStateCovariance(P);
    }

    Matrix EKFFlexibilityEstimatorBase::getFlexibilityCovariance
                                                            () const
    {
        return ekf_.getStateCovariance();
    }

    void EKFFlexibilityEstimatorBase::setMeasurement(const Vector & y)
    {
        ekf_.setMeasurement(y,k_+1);
    }

    Vector EKFFlexibilityEstimatorBase::getMeasurement()
    {
        return ekf_.getMeasurement(k_+1);
    }

    void EKFFlexibilityEstimatorBase::setProcessNoiseCovariance(const Matrix & Q)
    {
        ekf_.setQ(Q);
    }

    void EKFFlexibilityEstimatorBase::setMeasurementNoiseCovariance(const Matrix & R)
    {
        ekf_.setR(R);
    }

    ///gets the covariance matrices for the process noises
    Matrix EKFFlexibilityEstimatorBase::getProcessNoiseCovariance() const
    {
        return ekf_.getQ();
    }

    ///gets the covariance matrices for the sensor noises
    Matrix EKFFlexibilityEstimatorBase::getMeasurementNoiseCovariance() const
    {
        return ekf_.getR();
    }

    void EKFFlexibilityEstimatorBase::useFiniteDifferencesJacobians(Vector x)
    {
        finiteDifferencesJacobians_= true ;
        dx_ = x;
    }

    void EKFFlexibilityEstimatorBase::setJacobians
                                          (const Matrix & A, const Matrix & C)
    {
        finiteDifferencesJacobians_= false ;
        ekf_.setA(A);
        ekf_.setC(C);
    }

    void EKFFlexibilityEstimatorBase::setInput(const Vector & u)
    {
        ekf_.setInput(u,k_);
    }


    void EKFFlexibilityEstimatorBase::setMeasurementInput(const Vector & u)
    {
        if (ekf_.getInputsNumber()==0)
          ekf_.setInput(u,k_);
        ekf_.setInput(u,k_+1);
    }

    Vector EKFFlexibilityEstimatorBase::getInput()
    {
        return ekf_.getInput(k_);
    }


    Vector EKFFlexibilityEstimatorBase::getMeasurementInput()
    {
        return ekf_.getInput(k_+1);
    }

    const Vector & EKFFlexibilityEstimatorBase::getFlexibilityVector()
    {
        if (ekf_.getMeasurementsNumber()>0)
        {
            k_=ekf_.getMeasurementTime();
           // std::cout << "\n\n\n\n\n k " << k_ << std::endl;

            for (TimeIndex i=ekf_.getCurrentTime()+1; i<=k_; ++i)
            {
                if (finiteDifferencesJacobians_)
                {
                  ekf_.setA(ekf_.getAMatrixFD(dx_));
                  ekf_.setC(ekf_.getCMatrixFD(dx_));
                }
                ekf_.getEstimatedState(i);
            }

            Vector x(ekf_.getEstimatedState(k_));

            if (x==x)//detect NaN values
            {
                lastX_=x;
            }
            else //delete NaN values
            {
                ekf_.setState(lastX_,k_);
               // std::cout << "\n\n\n\n\n Need to reset covariance matrix \n\n\n\n\n" << std::endl;
//                std::cout << "k_: " << k_ << std::endl;
                if(k_>1) //the first iteration give always nan without initialization
                {
                    resetCovarianceMatrices();
                }
            }
        }
        return lastX_;
    }

    const stateObservation::ExtendedKalmanFilter &
        EKFFlexibilityEstimatorBase::getEKF() const
    {
        return ekf_;
    }

    stateObservation::ExtendedKalmanFilter & EKFFlexibilityEstimatorBase::getEKF()
    {
        return ekf_;
    }

    Vector EKFFlexibilityEstimatorBase::getSimulatedMeasurement()
    {

        getFlexibilityVector();

        return ekf_.getSimulatedMeasurement(ekf_.getCurrentTime());
    }

    Vector EKFFlexibilityEstimatorBase::getInnovation()
    {
        return ekf_.getInnovation();
    }

    Vector EKFFlexibilityEstimatorBase::getPredictedMeasurement()
    {
        return ekf_.updateStateAndMeasurementPrediction();
    }

    Vector EKFFlexibilityEstimatorBase::getPrediction()
    {
        return ekf_.updateStatePrediction();
    }


    Vector EKFFlexibilityEstimatorBase::getLastPredictedMeasurement()
    {
        return ekf_.getLastPredictedMeasurement();
    }

    Vector EKFFlexibilityEstimatorBase::getLastPrediction()
    {
        return ekf_.getLastPrediction();
    }

}
}
