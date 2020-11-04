IndexedVectorArray imuAttitudeTrajectoryReconstruction
    (
    const IndexedVectorArray & y,
    const IndexedVectorArray & u,
    const Vector & xh0,
    const Matrix & p,
    const Matrix & q,
    const Matrix & r,
    double dt)
{
    ///Sizes of the states for the state, the measurement, and the input vector
    const Index stateSize=18;
    const Index measurementSize=6;
    const Index inputSize=6;

    typedef kine::indexes<kine::rotationVector> indexes;

    ///initialization of the extended Kalman filter
    ExtendedKalmanFilter filter(stateSize, measurementSize, inputSize, false);

    ///initalization of the functor
    IMUDynamicalSystem imuFunctor;
    imuFunctor.setSamplingPeriod(dt);
    filter.setFunctor(& imuFunctor);

    ///the initalization of the estimation of the initial state
    filter.setState(xh0,y.getFirstIndex()-1);

    ///computation and initialization of the covariance matrix of the initial state
    filter.setStateCovariance(p);

    ///set initial input
    filter.setInput(u[y.getFirstIndex()-1],y.getFirstIndex()-1);

    ///The covariance matrix of the process noise and the measurement noise
    /// for the extended Kalman filter
    filter.setR(r);
    filter.setQ(q);

    ///set the derivation step for the finite difference method
    Vector dx=filter.stateVectorConstant(1)*1e-8;

    ///the array of the state estimations over time
    IndexedVectorArray xh;
    xh.setValue(xh0,y.getFirstIndex()-1);

    ///the reconstruction of the state
    for (TimeIndex i=y.getFirstIndex();i<y.getNextIndex();++i)
    {
        ///introduction of the measurement
        filter.setMeasurement(y[i],i);

        ///introduction of the input
        if (i<y.getLastIndex())
            filter.setInput(u[i],i);

        ///get the jacobians by finite differences and provide
        ///them to the Kalman filter
        Matrix a=filter.getAMatrixFD(dx);
        Matrix c= filter.getCMatrixFD(dx);
        filter.setA(a);
        filter.setC(c);

        ///get the estimation and give it to the array
        Vector xhk=filter.getEstimatedState(i);

        ///regulate the part of orientation vector in the state vector
        xhk.segment(indexes::ori,3)=kine::regulateOrientationVector
                                                    (xhk.segment(indexes::ori,3));

        ///give the new value of the state to the kalman filter.
        ///This step is usually unnecessary, unless we modify the
        ///value of the state esimation which is the case here.
        filter.setState(xhk,i);

        xh.setValue(xhk,i);
    }

    return xh;
}

IndexedVectorArray imuAttitudeTrajectoryReconstruction(
    const IndexedVectorArray & y,
    const Vector & xh0,
    const Matrix & p,
    const Matrix & q,
    const Matrix & r,
    double dt)
{
    const Index inputSize=6;

    ///initialization of a zero input
    IndexedVectorArray u;
    for (TimeIndex k=y.getFirstIndex()-1; k<y.getNextIndex(); ++k)
    {
        u.setValue(Vector::Zero(inputSize,1),k);
    }

    return imuAttitudeTrajectoryReconstruction (y, u, xh0, p, q, r, dt);
}
