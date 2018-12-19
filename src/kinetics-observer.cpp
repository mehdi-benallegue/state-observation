#include <state-observation/dynamics-estimators/kinetics-observer.hpp>





namespace stateObservation
{

  KineticsObserver::KineticsObserver(int maxContacts):
  stateSize_(stateSizeBase + maxContacts*stateSizePerContact),
  stateTangentSize_(stateTangentSizeBase + stateTangentSizePerContact* maxContacts),
  measurementSize_(0),
  stateVector_(stateSize_),
  stateVectorDx_(Vector::Constant(stateSizeBase + maxContacts*stateSizePerContact,defaultdx)),
  oldStateVector_(stateSize_),
  acceleroDefaultCovMat_(Matrix3::Identity()*acceleroVariance),
  gyroCovDefaultMat_(Matrix3::Identity()*gyroVariance),
  contactFTSensorDefaultCovMat_(Matrix6::Identity()),
  poseSensorCovMat_(Matrix6::Identity()),
  ekf_(stateSize_, stateTangentSize_, sizeIMUSignal, sizeIMUSignal,  0,false,false),
  finiteDifferencesJacobians_(true),
  k_est(0),k_data(0)
  {
    contactFTSensorDefaultCovMat_.block<3,3>(0,0) *= forceSensorVariance;
    contactFTSensorDefaultCovMat_.block<3,3>(3,3) *= torqueSensorVariance;
    poseSensorCovMat_.block<3,3>(3,3) *= positionSensorVariance;
    poseSensorCovMat_.block<3,3>(3,3) *= orientationSensorVariance;

    ekf_.setFunctor(this);
    
    Contact::numberOfRealSensors = 0;
  }

  KineticsObserver::~KineticsObserver()
  {}

  void KineticsObserver::update()
  {

    if (k_est!=k_data)
    {
      /// Synchronizing the sensors
      MapIMUIterator i = imuSensors_.begin();
      while (i != imuSensors_.end()) 
      {
        if (i->second.time!=k_data) 
        {
          imuSensors_.erase(i++); /// remove the i-th element and move to the next
        } 
        else 
        {
          ++i;
        }
      }
      
      for (MapContactIterator i= contacts_.begin(), ie = contacts_.end(); i!=ie ; ++i) 
      {
        BOOST_ASSERT((i->second.index == k_data) && "Exception in Kinetics Observer: the contacts have not all been updated. \
              Either remove lost contacts using removeContact \
              or Run setContactFTSensor or setContactWithNoSensor on every added contact");

        /// the following code is only an attempt to maintain a coherent state of the state observer
        if (i->second.index != k_data && i->second.withRealSensor)
        {
          i->second.withRealSensor=false;
          Contact::numberOfRealSensors--;
        }
      }

      if (absPoseSensor_.time != k_data)
      {
        absPoseSensor_.isSet=false;
      }

      ///////////// initialize the measurement Vector and matrix //////////////

      measurementSize_ = sizeIMUSignal*int(imuSensors_.size()) + sizeFTSignal*Contact::numberOfRealSensors;
      int measurementTangentSize = measurementSize_;
      if (absPoseSensor_.isSet)
      {
        measurementSize_ +=sizePoseSignal;
        measurementTangentSize+=sizePoseSignalTangent;
      }
      
      measurementVector_.resize(measurementSize_);
      measurementCovMatrix_.resize(measurementTangentSize,measurementTangentSize);
      measurementCovMatrix_.setZero();
    
      int localIndex = 0;
      
      for (MapIMUIterator i= imuSensors_.begin(), ie = imuSensors_.end(); i!=ie ; ++i)
      {
        i->second.index = localIndex;
        measurementVector_.segment<sizeIMUSignal>(localIndex) = i->second.acceleroGyro;
        measurementCovMatrix_.block<sizeAcceleroSignal,sizeAcceleroSignal>(localIndex,localIndex)=i->second.covMatrixAccelero;
        localIndex+=sizeAcceleroSignal;
        measurementCovMatrix_.block<sizeGyroSignal,sizeGyroSignal>(localIndex,localIndex)=i->second.covMatrixGyro;
        localIndex+=sizeGyroSignal;
      }

      for (MapContactIterator i=contacts_.begin(), ie = contacts_.end();i!=ie;++i) 
      {
        if (i->second.withRealSensor)
        {
          i->second.index = localIndex;
          measurementVector_.segment<sizeFTSignal>(localIndex) = i->second.forceTorque;
          measurementCovMatrix_.block<sizeFTSignal,sizeFTSignal>(localIndex,localIndex)=i->second.covMatrix;
          localIndex+=sizeFTSignal;
        }
      }

      if (absPoseSensor_.isSet)
      {
        absPoseSensor_.index= localIndex;
        BOOST_ASSERT(absPoseSensor_.pose.position.isSet() && absPoseSensor_.pose.orientation.isSet() \
                    && "The absolute pose needs to contain the position and the orientation")
        measurementVector_.segment<sizePoseSignal>(localIndex) = absPoseSensor_.pose.toVector(Kinematics::Flags::position | Kinematics::Flags::orientation);
        measurementCovMatrix_.block<sizePoseSignalTangent,sizePoseSignalTangent>(localIndex,localIndex)=absPoseSensor_.covMatrix;
      }

      ekf_.setMeasureSize(measurementSize_,measurementTangentSize);
      ekf_.setMeasurement(measurementVector_);
      ekf_.setR(measurementCovMatrix_);
      if (useFiniteDifferencesJacobians)
      {
        ekf_.setA(ekf_.getAMatrixFD(dx_));
        ekf_.setC(ekf_.getCMatrixFD(dx_));
      }
      
      stateVector_ = ekf_.getEstimatedState(k_data);


      if (stateVector_.hasNaN())
      {
  #ifndef NDEBUG
        std::cout << "Kinetics observer: NaN value detected" << std::endl;
  #endif 
      }
      else
      {
        oldStateVector_ = stateVector_;
      }

      ++k_est; //the timestamp of the state we estimated

    }
  }

  int KineticsObserver::setIMU(const Vector3 & accelero, const  Vector3 & gyrometer, const Kinematics &localKine, int num)
  {
    ///ensure the measuements are labeled with the good time stamp
    startNewIteration_();    
  }

  void KineticsObserver::setContactFTSensor(const Vector3& force, const Vector3& torque, const Matrix6 ForcetorqueCovMatrix, 
                                                                                    const Kinematics &localKine, int contactNumber)
  {
    startNewIteration_();
  }
   
  void KineticsObserver::setContactWithNoSensor(const Kinematics &localKine, int contactNumber)
  {
    startNewIteration_();
  }

  void KineticsObserver::setPoseSensor(const Kinematics &, const Matrix6 CovarianceMatrix)
  {
    startNewIteration_();
  }

  void KineticsObserver::resetIteration()
  {
    imuSensors_.clear();
    contacts_.clear();
    absPoseSensor_.isSet=false;
  }

  void KineticsObserver::startNewIteration_()
  {
    if (k_est==k_data)
    {
      ++k_data;
    }
    
  }

  
}
