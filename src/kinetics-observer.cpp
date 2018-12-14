
#include <state-observation/dynamics-estimators/kinetics-observer.hpp>

static double defaultdx=1e-6;

namespace stateObservation
{
  KineticsObserver::KineticsObserver(int maxContacts):
  acceleroDefaultCovMat_(Matrix3::Identity()*acceleroVariance),
  gyroCovDefaultMat_(Matrix3::Identity()*gyroVariance),
  contactFTSensorCovMat_(Matrix6::Identity()),
  poseSensorCovMat_(Matrix6::Identity()),
  ekf_(stateSizeBase + maxContacts*stateSizePerContact, 
        stateTangentSizeBase + stateTangentSizePerContact* maxContacts,
        signalSizeIMU,0,false,false),
  k_est(0),
  k_input(0)
  {
    contactFTSensorCovMat_.block<3,3>(0,0) *= forceSensorVariance;
    contactFTSensorCovMat_.block<3,3>(3,3) *= torqueSensorVariance;
    poseSensorCovMat_.block<3,3>(3,3) *= positionSensorVariance;
    poseSensorCovMat_.block<3,3>(3,3) *= orientationSensorVariance;

    dx_ = Vector::Constant(stateSizeBase + maxContacts*stateSizePerContact,defaultdx)  ;
  }

  void KineticsObserver::update()
  {



    ///last instruction
    ++k_est;
  }

  int KineticsObserver::setIMU(const Vector3 & accelero, const  Vector3 & gyrometer, const Kinematics &localKine, int num)
  {
    ///ensure the measuements are labeled with the good time stamp
    onSetMeasurement_();
    
  }

  void KineticsObserver::onSetMeasurement_()
  {
    if (k_input!=k_est) ///This is a new iteration
    {
      resetIteration();
      k_input = k_est;
    }
  }

  void KineticsObserver::resetIteration()
  {
    mapIMU_.clear();
    mapFT_.clear();
    
  }

  
}




