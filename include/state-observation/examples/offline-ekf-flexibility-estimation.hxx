stateObservation::IndexedVectorArray offlineEKFFlexibilityEstimation(
            const stateObservation::IndexedVectorArray & y,
            const stateObservation::IndexedVectorArray & u,
            const Vector & xh0,
            unsigned numberOfContacts,
            const std::vector<Vector3,Eigen::aligned_allocator<Vector3> > & contactsPositions,
            double dt,
            IndexedVectorArray * ino,
            IndexedVectorArray * premea)
{

    (void)dt;
        ///Sizes of the states for the state, the measurement, and the input vector

    flexibilityEstimation::FixedContactEKFFlexEstimatorIMU estimator;


    Matrix R(estimator.getEKF().getRmatrixIdentity());
    R=R*1.e-30;

    Matrix Q(estimator.getEKF().getQmatrixIdentity());
    Q=Q*1.e-4;
    Q.block(6,6,3,3)=Matrix3::Identity()*1.e-8;
    Q.block(15,15,3,3)=Matrix3::Identity()*1.e-8;


    estimator.setFlexibilityCovariance(Q);
    estimator.setMeasurementNoiseCovariance(R);
    estimator.setProcessNoiseCovariance(Q);


    estimator.setContactsNumber(numberOfContacts);

    for (unsigned i = 0; i<numberOfContacts ; ++i)
    {
        estimator.setContactPosition(i,contactsPositions[i]);
    }

    ///
    estimator.setFlexibilityGuess(xh0);

    ///the array of the state estimations over time
    stateObservation::IndexedVectorArray xh;
    xh.setValue(xh0,y.getFirstIndex()-1);

    ///the reconstruction of the state
    for (TimeIndex i=y.getFirstIndex();i<y.getNextIndex();++i)
    {
        //std::cout << i << std::endl;

        ///introduction of the measurement
        estimator.setMeasurement(y[i]);

        estimator.setMeasurementInput(u[i]);

        ///get the estimation and give it to the array
        xh.pushBack(estimator.getFlexibilityVector());

        if (ino != 0)
        {
            ino->setValue(estimator.getInnovation(),i);
        }

        if (premea != 0)
        {
            premea->setValue(estimator.getLastPredictedMeasurement(),i);
        }

    }

    return xh;
}


stateObservation::IndexedVectorArray offlineEKFFlexibilityEstimation(
            const stateObservation::IndexedVectorArray & y,
            const Vector & xh0,
            unsigned numberOfContacts,
            const std::vector<Vector3, Eigen::aligned_allocator<Vector3> > & contactsPositions,
            double dt)
{
    const unsigned inputSize=15;

    ///initialization of a zero input
    stateObservation::IndexedVectorArray u;
    for (TimeIndex k=y.getFirstIndex()-1; k<y.getNextIndex(); ++k)
    {
        u.setValue(Vector::Zero(inputSize,1),k);
    }

    return offlineEKFFlexibilityEstimation
                        (y, u, xh0, numberOfContacts, contactsPositions, dt);
}

