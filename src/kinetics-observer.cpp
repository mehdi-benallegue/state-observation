#include <state-observation/dynamics-estimators/kinetics-observer.hpp>





namespace stateObservation
{
  inline Matrix6 blockMat6(const Matrix3 & m1, const Matrix3 & m2, const Matrix3 & m3, const Matrix3 & m4)
  {
    Matrix6 m;
    m<< m1,m2,
        m3,m4;
    
    return m;
  }

  const double KineticsObserver::defaultMass = 50;

  const double KineticsObserver::statePoseInitVarianceDefault = 1e-4;
  const double KineticsObserver::stateOriInitVarianceDefault = 1e-4;
  const double KineticsObserver::stateLinVelInitVarianceDefault = 1e-6;
  const double KineticsObserver::stateAngVelInitVarianceDefault = 1e-6;
  const double KineticsObserver::gyroBiasInitVarianceDefault = 1e-10;
  const double KineticsObserver::unmodeledWrenchInitVarianceDefault = 1e100;
  const double KineticsObserver::contactForceInitVarianceDefault = 1e100;
  const double KineticsObserver::contactTorqueInitVarianceDefault = 1e100;

  const double KineticsObserver::statePoseProcessVarianceDefault = 1e-8;
  const double KineticsObserver::stateOriProcessVarianceDefault = 1e-8;
  const double KineticsObserver::stateLinVelProcessVarianceDefault = 1e-8;
  const double KineticsObserver::stateAngVelProcessVarianceDefault = 1e-8;
  const double KineticsObserver::gyroBiasProcessVarianceDefault = 1e-12;
  const double KineticsObserver::unmodeledWrenchProcessVarianceDefault = 1e-8;
  const double KineticsObserver::contactForceProcessVarianceDefault = 1e-8;
  const double KineticsObserver::contactTorqueProcessVarianceDefault = 1e-8;

  const double KineticsObserver::acceleroVarianceDefault = 1e-4;
  const double KineticsObserver::gyroVarianceDefault = 1e-8;
  const double KineticsObserver::forceSensorVarianceDefault = 1e-8;
  const double KineticsObserver::torqueSensorVarianceDefault = 1e-10;
  const double KineticsObserver::positionSensorVarianceDefault = 1e-4;
  const double KineticsObserver::orientationSensorVarianceDefault = 1e-3;

  const double linearStiffnessDefault = 40000;
  const double angularStiffnessDefault = 400;
  const double linearDampingDefault = 120;
  const double angularDampingDefault = 12;

  const double KineticsObserver::defaultdx = 1e-6;

  KineticsObserver::KineticsObserver(int maxContacts):
    stateSize_(sizeStateBase + maxContacts*sizeStatePerContact),
    stateTangentSize_(sizeStateTangentBase + sizeStateTangentPerContact* maxContacts),
    measurementSize_(0),
    stateVector_(stateSize_),
    stateVectorDx_(stateTangentSize_),
    oldStateVector_(stateSize_),
    additionalWrench_(Vector6::Zero()),    
    ekf_(stateSize_, stateTangentSize_, sizeIMUSignal, sizeIMUSignal,  0,false,false),
    finiteDifferencesJacobians_(true),
    withGyroBias_(true),
    withUnmodeledWrench_(false),
    k_est(0),k_data(0), mass_(defaultMass), dt_(defaultdx),
    linearStiffnessMatDefault_(Matrix3::Identity()*linearStiffnessDefault),
    angularStiffnessMatDefault_(Matrix3::Identity()*angularStiffnessDefault),
    linearDampingMatDefault_(Matrix3::Identity()*linearDampingDefault),
    angularDampingMatDefault_(Matrix3::Identity()*angularDampingDefault),
    acceleroCovMatDefault_(Matrix3::Identity()*acceleroVarianceDefault),
    gyroCovMatDefault_( Matrix3::Identity()*gyroVarianceDefault),
    contactWrenchSensorCovMatDefault_(blockMat6( Matrix3::Identity()*forceSensorVarianceDefault, Matrix3::Zero(),
                                Matrix3::Zero(), Matrix3::Identity()*torqueSensorVarianceDefault )),
    poseSensorCovMatDefault_(blockMat6( Matrix3::Identity()*positionSensorVarianceDefault, Matrix3::Zero(),
                                Matrix3::Zero(), Matrix3::Identity()*orientationSensorVarianceDefault )),
    statePosInitCovMat_(Matrix3::Identity()*statePoseInitVarianceDefault),
    stateOriInitCovMat_(Matrix3::Identity()*stateOriInitVarianceDefault),
    stateLinVelInitCovMat_(Matrix3::Identity()*stateLinVelInitVarianceDefault),
    stateAngVelInitCovMat_(Matrix3::Identity()*stateAngVelInitVarianceDefault),
    gyroBiasInitCovMat_(Matrix3::Identity()*gyroBiasInitVarianceDefault),
    unmodeledWrenchInitCovMat_(Matrix6::Identity()*unmodeledWrenchInitVarianceDefault),
    contactWrenchInitCovMat_(blockMat6(Matrix3::Identity()*contactForceInitVarianceDefault, Matrix3::Zero(),
                                Matrix3::Zero(), Matrix3::Identity()*contactTorqueInitVarianceDefault)),
    statePosProcessCovMat_(Matrix3::Identity()*statePoseProcessVarianceDefault),
    stateOriProcessCovMat_(Matrix3::Identity()*stateOriProcessVarianceDefault),
    stateLinVelProcessCovMat_(Matrix3::Identity()*stateLinVelProcessVarianceDefault),
    stateAngVelProcessCovMat_(Matrix3::Identity()*stateAngVelProcessVarianceDefault),
    gyroBiasProcessCovMat_(Matrix3::Identity()*gyroBiasProcessVarianceDefault),
    unmodeledWrenchProcessCovMat_(Matrix6::Identity()*unmodeledWrenchProcessVarianceDefault),
    contactWrenchProcessCovMat_(blockMat6(Matrix3::Identity()*contactForceProcessVarianceDefault, Matrix3::Zero(),
                                Matrix3::Zero(), Matrix3::Identity()*contactTorqueProcessVarianceDefault))
  {
    ekf_.setFunctor(this);

    stateVector_.setZero();

    ekf_.setState(stateVector_,k_est);

    const Matrix3 & Zero3 = Matrix3::Zero();

    stateKineMatricsInitCovMat_ << 
        statePosInitCovMat_, Zero3	            , Zero3	                , Zero3,
        Zero3              , stateOriInitCovMat_, Zero3	                , Zero3,
        Zero3              , Zero3              , stateLinVelInitCovMat_, Zero3,
        Zero3              , Zero3              , Zero3	                , stateAngVelInitCovMat_;


    stateKineMatricsProcessCovMat_ << 
        statePosProcessCovMat_, Zero3	                , Zero3	                   , Zero3,
        Zero3                 , stateOriProcessCovMat_, Zero3	                   , Zero3,
        Zero3                 , Zero3                 , stateLinVelProcessCovMat_, Zero3,
        Zero3                 , Zero3                 , Zero3	                   , stateAngVelProcessCovMat_;

    resetStateCovarianceMat();
    resetProcessCovarianceMat();

    updateKine_();

    Contact::numberOfRealSensors = 0;
  }

  KineticsObserver::~KineticsObserver()
  {
  }

  
  unsigned KineticsObserver::getStateSize() const
  {
    return stateSize_;
  }

  unsigned KineticsObserver::getMeasurementSize() const
  {
    int size = 0;
  /// Synchronizing the sensors
    MapIMUConstIterator i = imuSensors_.begin();
    while (i != imuSensors_.end()) 
    {
      if (i->second.time==k_data) 
      {
        size += sizeIMUSignal;
      } 
    }
    
    for (MapContactConstIterator i= contacts_.begin(), ie = contacts_.end(); i!=ie ; ++i) 
    {
      if (i->second.index == k_data && i->second.withRealSensor)
      {
        size += sizeWrench;
      }
    }

    if (absPoseSensor_.isSet && absPoseSensor_.time == k_data)
    {
      size += sizePose;
    }
    return size;
  }

  double KineticsObserver::getSamplingTime() const
  {
    return dt_;
  }

  void KineticsObserver::setSamplingTime( double dt) 
  {
    dt_ = dt;
  }

  void KineticsObserver::setMass(double m)
  {
    mass_=m;
  }
  
  Vector KineticsObserver::update()
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

      if (absPoseSensor_.isSet && absPoseSensor_.time != k_data)
      {
        absPoseSensor_.isSet=false;
      }

      ///////////// initialize the measurement Vector and matrix //////////////

      measurementSize_ = sizeIMUSignal*int(imuSensors_.size()) + sizeWrench*Contact::numberOfRealSensors;
      int measurementTangentSize = measurementSize_;
      if (absPoseSensor_.isSet)
      {
        measurementSize_ +=sizePose;
        measurementTangentSize+=sizePoseTangent;
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
          measurementVector_.segment<sizeWrench>(localIndex) = i->second.forceTorque;
          measurementCovMatrix_.block<sizeWrench,sizeWrench>(localIndex,localIndex)=i->second.covMatrix;
          localIndex+=sizeWrench;
        }
      }

      if (absPoseSensor_.isSet)
      {
        absPoseSensor_.index= localIndex;
        BOOST_ASSERT(absPoseSensor_.pose.position.isSet() && absPoseSensor_.pose.orientation.isSet() \
                    && "The absolute pose needs to contain the position and the orientation");
        measurementVector_.segment<sizePose>(localIndex) = absPoseSensor_.pose.toVector(flagsKineSensor);
        measurementCovMatrix_.block<sizePoseTangent,sizePoseTangent>(localIndex,localIndex)=absPoseSensor_.covMatrix;
      }

      ekf_.setMeasureSize(measurementSize_,measurementTangentSize);
      ekf_.setMeasurement(measurementVector_,k_data);
      ekf_.setR(measurementCovMatrix_);
      if (finiteDifferencesJacobians_)
      {
        ekf_.setA(ekf_.getAMatrixFD(stateVectorDx_));
        ekf_.setC(ekf_.getCMatrixFD(stateVectorDx_));
      }

      stateVector_ = ekf_.getEstimatedState(k_data);

      if (stateVector_.hasNaN())
      {
  #ifndef NDEBUG
        std::cout << "Kinetics observer: NaN value detected" << std::endl;
  #endif 
        stateVector_ = stateNaNCorrection_();
      }
      else
      {
        oldStateVector_ = stateVector_;
      }

      ++k_est; //the timestamp of the state we estimated

      return stateVector_;
    }
  }

  Vector KineticsObserver::getStateVector() const
  {
    return stateVector_;
  }

  kine::Kinematics KineticsObserver::getKinematics() const
  {
    return stateKinematics_; 
  }  

  kine::Kinematics KineticsObserver::getKinematics( const Kinematics & local) const
  {
    return getKinematics()*local;
  }

  Vector6 KineticsObserver::getContactWrench(int contactNbr) const
  {
    return stateVector_.segment<sizeWrench>(contactWrenchIndex(contactNbr));
  }

  kine::Kinematics KineticsObserver::getContactPosition(int contactNbr) const
  {
    return Kinematics(stateVector_.segment<sizeStateKine>(contactKineIndex(contactNbr)),
                                                          flagsContactKine);
  }

  Vector6 KineticsObserver::getUnmodeledWrench() const
  {
    return stateVector_.segment<sizeWrench>(unmodeledWrenchIndex());
  }

  kine::Kinematics KineticsObserver::estimateAccelerations()
  {
    if (!stateKinematics_.linAcc.isSet())
    {
      stateKinematics_.fromVector(computeAccelerations_(),
                                  Kinematics::Flags::linAcc|
                                  Kinematics::Flags::angAcc);
    }

    return stateKinematics_;
  }

  void KineticsObserver::setStateKinematics(const Kinematics & kine, bool resetForces,
                                            bool resetCovariance)
  {
    stateKinematics_ = kine;
    stateVector_.segment<sizeStateKine>(kineIndex()) = stateKinematics_.toVector(flagsStateKine);

    if (resetForces)
    {
      for (MapContactIterator i = contacts_.begin(); i != contacts_.end(); ++i)
      {
        stateVector_.segment<sizeWrench>(contactWrenchIndex(i->first)).setZero();
      }
    }

    ekf_.setState(stateVector_,k_est);

    if (resetCovariance)
    {
      Matrix stateCovariance = ekf_.getStateCovariance();
      setBlockStateCovariance<sizeStateKine>(stateCovariance,stateKineMatricsInitCovMat_,kineIndex(),stateSize_);

      if (resetForces)
      {
        for (MapContactIterator i = contacts_.begin(); i != contacts_.end(); ++i)
        {
          setBlockStateCovariance<sizeStateKine>(stateCovariance,contactWrenchInitCovMat_,contactWrenchIndex(i->first),stateSize_);
        }
      }
      ekf_.setStateCovariance(stateCovariance);
    }

  }

  void KineticsObserver::setGyroBias(const Vector3 &,  bool resetCovariance)
  {
    stateVector_.segment<sizeGyroBias>(gyroBiasIndex());
    ekf_.setState(stateVector_,k_est);

    if (resetCovariance)
    {
      Matrix stateCovariance = ekf_.getStateCovariance();
      setBlockStateCovariance<sizeGyroBias>(stateCovariance,gyroBiasInitCovMat_,gyroBiasIndex(),stateSize_);
      
      ekf_.setStateCovariance(stateCovariance);
    }
  }

  void KineticsObserver::setStateUnmodeledWrench(const Vector6 &, bool resetCovariance)
  {
    stateVector_.segment<sizeWrench>(unmodeledWrenchIndex());
    ekf_.setState(stateVector_,k_est);

    if (resetCovariance)
    {
      Matrix stateCovariance = ekf_.getStateCovariance();
      setBlockStateCovariance<sizeWrench>(stateCovariance,unmodeledWrenchInitCovMat_,unmodeledWrenchIndex(),stateSize_);
      
      ekf_.setStateCovariance(stateCovariance);
    }
  }


  void KineticsObserver::setStateVector(const Vector & v, bool resetCovariance)
  {
    stateVector_ = v;
    ekf_.setState(v,k_est);
    updateKine_();

    if (resetCovariance)
    {
      resetStateCovarianceMat();
    }
  }

  void KineticsObserver::setAdditionalWrench(const Vector6& wrench)
  {
    
    additionalWrench_=wrench;
  }

///////////////////////
//////////////////////////
///////////////////////


  Vector KineticsObserver::stateNaNCorrection_()
  {
    ///TODO implement this function
    assert(false && "NaN Correction not yet implemented. Please Contact mehdi.benallegue@gmail.com");
    return oldStateVector_;
  }


  int KineticsObserver::setIMU(const Vector3 & accelero, const  Vector3 & gyrometer, const Kinematics &localKine, int num)
  {
    ///ensure the measuements are labeled with the good time stamp
    startNewIteration_();    
  }

  void KineticsObserver::setContactWrenchSensor(const Vector3& force, const Vector3& torque, const Matrix6 ForcetorqueCovMatrix, 
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

  void KineticsObserver::resetInputs()
  {
    imuSensors_.clear();
    absPoseSensor_.isSet=false;

    for (MapContactIterator i=contacts_.begin(), ie = contacts_.end();i!=ie;++i) 
    {
      i->second.withRealSensor=false;
    }
  }

  void KineticsObserver::startNewIteration_()
  {
    if (k_est==k_data)
    {
      ++k_data;
    }    
  } 

  void KineticsObserver::updateKine_()
  {
    stateKinematics_.fromVector(stateVector_.segment<sizeStateKine>(kineIndex()), 
                                flagsStateKine);
  } 
}
