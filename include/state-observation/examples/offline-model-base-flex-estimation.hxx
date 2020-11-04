//const double acc_cov_const=5e-3;
//const double gyr_cov_const=5e-6;
//const double force_sensor_const=1e-10;
//const double torque_sensor_const=1e-30;
//const double state_fc_const=5e-4;
//const double pos_state_cov_const=0;
//const double vel_state_cov_const=1e-10;
//const double ori_state_cov_const=0;
//const double angv_state_cov_const=1e-8;

typedef flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::state state;

stateObservation::IndexedVectorArray offlineModelBaseFlexEstimation(
  const stateObservation::IndexedVectorArray & y,
  const stateObservation::IndexedVectorArray & u,
  const Matrix & xh0,
  const stateObservation::IndexedVectorArray numberOfContacts,
  double dt,
  double mass,
  bool withForce,
  const stateObservation::IndexedMatrixArray & Q,
  const stateObservation::IndexedMatrixArray & R,
  const Matrix3 & kfe,
  const Matrix3 & kfv,
  const Matrix3 & kte,
  const Matrix3 & ktv,
  IndexedVectorArray * prediction,
  IndexedVectorArray * ino,
  IndexedVectorArray * premea,
  IndexedVectorArray * simumea,
  int verbose)
{


  flexibilityEstimation::ModelBaseEKFFlexEstimatorIMU estimator;

  estimator.setContactModel(flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::contactModel::elasticContact);
  estimator.setSamplingPeriod(dt);
  estimator.setRobotMass(mass);

  bool customQ = false,
       customR = false;

  if (Q.size()==1)
  {
    estimator.setProcessNoiseCovariance(Q[Q.getFirstIndex()]);
  }
  if (R.size()==1)
  {
    estimator.setMeasurementNoiseCovariance(R[R.getFirstIndex()]);
  }

  if (Q.size()>1)
  {
    customQ=true;
  }
  if (R.size()>1)
  {
    customR=true;
  }



//  Matrix Ri,Qi,P;
//
//  int measurementSize=6;
//  int stateSize=flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::state::size;
//
//  double acceleroCovariance=acc_cov_const;
//  double gyroCovariance=gyr_cov_const;
//  double posStateCov=pos_state_cov_const;
//  double oristateCov=ori_state_cov_const;
//  double velStateCov=vel_state_cov_const;
//  double angvStateCov=angv_state_cov_const;
//  double stateForceCov=state_fc_const;



//  Ri.noalias()=Matrix::Identity(measurementSize,measurementSize)*acceleroCovariance;
//  Ri(3,3)=R(4,4)=R(5,5)=gyroCovariance;
//  Qi.noalias()=Matrix::Identity(stateSize,stateSize)*posStateCov;
//
//  Qi.diagonal().segment<3>(state::ori).setConstant(oristateCov);
//
//  Qi.diagonal().segment<3>(state::linVel).setConstant(velStateCov);
//  Qi.diagonal().segment<3>(state::angVel).setConstant(angvStateCov);
//  Qi.diagonal().segment<12>(state::fc).setConstant(stateForceCov);
//
//  estimator.setProcessNoiseCovariance(Qi);
//  estimator.setMeasurementNoiseCovariance(Ri);
//  Matrix forcevariance= Matrix::Identity(6,6)*force_sensor_const;
//  forcevariance.block(3,3,3,3)=Matrix::Identity(3,3)*torque_sensor_const;
//  /// consider the vertical torque as reliable as a force
//  forcevariance(2,2)=torque_sensor_const;
//  /// consider the vertical force as reliable as torque
//  forcevariance(6,6)=force_sensor_const;
//  estimator.setForceVariance(forcevariance);



  estimator.setFlexibilityGuess(xh0);

  estimator.setWithUnmodeledForces(true);


  estimator.setWithForcesMeasurements(withForce);


  if (kfe!=Matrix3::Zero())
  {
    estimator.setKfe(kfe);
    estimator.setKfv(kfv);
    estimator.setKte(kte);
    estimator.setKtv(ktv);
  }

  ///the array of the state estimations over time
  stateObservation::IndexedVectorArray xh;
  if (y.getFirstIndex()>0)
    xh.setValue(xh0,y.getFirstIndex()-1);

  ///the reconstruction of the state
  for (TimeIndex i=y.getFirstIndex(); i<y.getNextIndex(); ++i)
  {

    estimator.setContactsNumber(unsigned(round(numberOfContacts[i](0))));

    if (verbose>0)
    {
      if (i%1000==0 || verbose >1)
      {
        std::cout << "iteration: " << i  << std::endl ;

        if (verbose >2)
        {
          std::cout << "number of contacts: "<< numberOfContacts[i]<< std::endl;
          std::cout << "size of the measurement: "<< y[i].size()
                    << ", supposed to be " << estimator.getMeasurementSize() << std::endl;
          std::cout << "size of the input: "<< u[i].size()
                    << ", supposed to be " << estimator.getInputSize() << std::endl;

          if (verbose > 3)
          {
            std::cout << "numberOfContacts: " << numberOfContacts[i].transpose() << std::endl;
            std::cout << "Measurement: " << y[i].transpose() << std::endl;
            std::cout << "Input: " << u[i].transpose() << std::endl;
          }
        }
      }
    }

    ///introduction of the measurement
    estimator.setMeasurement(Vector(y[i]).head(estimator.getMeasurementSize()));

    estimator.setMeasurementInput(u[i]);

    ///initialize flexibility
    if (customQ)
    {
      estimator.setProcessNoiseCovariance(Q[i]);
    }

    if (customR)
    {
      estimator.setMeasurementNoiseCovariance(R[i]);
    }

    ///get the estimation and give it to the array
    Vector xhk=estimator.getFlexibilityVector();

    xh.pushBack(xhk);

    if (prediction != 0x0)
    {
      prediction->setValue(estimator.getLastPrediction(),i);
    }

    if (ino != 0x0)
    {
      ino->setValue(estimator.getInnovation(),i);
    }

    if (premea != 0x0)
    {
      premea->setValue(estimator.getLastPredictedMeasurement(),i);
    }

    if (simumea != 0x0)
    {
      simumea->setValue(estimator.getSimulatedMeasurement(),i);
    }

    if (verbose >1)
    {
      std::cout << "Success";
      if (verbose >2)
      {
        std::cout << ", rebuilt state "<<std::endl;
        std::cout << xhk.transpose();
      }

      std::cout<<std::endl;

    }

  }

  return xh;
}
