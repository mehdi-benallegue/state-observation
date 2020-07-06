#include <iostream>
#include <fstream>

#include <state-observation/flexibility-estimation/model-base-ekf-flex-estimator-imu.hpp>
#include <state-observation/flexibility-estimation/imu-elastic-local-frame-dynamical-system.hpp>


using namespace stateObservation;

inline Matrix3 buildInertiaTensor(const Vector6 inputInertia, Matrix3& inertiaTensor)
{

  const double & Ixx=inputInertia[0];
  const double & Iyy=inputInertia[1];
  const double & Izz=inputInertia[2];
  const double & Ixy=inputInertia[3];
  const double & Ixz=inputInertia[4];
  const double & Iyz=inputInertia[5];

  inertiaTensor   <<    Ixx, Ixy, Ixz,
                  Ixy, Iyy, Iyz,
                  Ixz, Iyz, Izz;

  return inertiaTensor;
}

Vector computeZmp (unsigned footNumber, IndexedMatrixArray& forces, IndexedMatrixArray& sensorPositions)
{
    double fnormal = 0;
    double sumZmpx = 0;
    double sumZmpy = 0;
    double sumZmpz = 0;
    Vector zmp; zmp.setZero(); zmp.resize (3);

    for (unsigned int i=0; i<footNumber; ++i)
    {
        const Vector& f = forces [i];
        // Check that force is of dimension 6
        if (f.size () != 6)
        {
            zmp.fill (0.);
            return zmp;
        }

        const Matrix& M = sensorPositions[i];
        double fx = M (0,0) * f (0) + M(0,1) * f (1) + M (0,2) * f (2);
        double fy = M (1,0) * f (0) + M(1,1) * f (1) + M (1,2) * f (2);
        double fz = M (2,0) * f (0) + M(2,1) * f (1) + M (2,2) * f (2);

        if (fz > 0)
        {
            double Mx = M (0,0)*f(3) + M (0,1)*f(4) + M (0,2)*f(5) + M (1,3)*fz - M (2,3)*fy;
            double My = M (1,0)*f(3) + M (1,1)*f(4) + M (1,2)*f(5) + M (2,3)*fx - M (0,3)*fz;

            fnormal += fz;
            sumZmpx -= My;
            sumZmpy += Mx;
            sumZmpz += fz * M (2,3);
        }
    }

    if (fnormal != 0)
    {
        zmp (0) = sumZmpx / fnormal;
        zmp (1) = sumZmpy / fnormal;
        zmp (2) = sumZmpz / fnormal;
    }
    else
    {
        zmp.fill (0.);
    }

    return zmp;
}

Vector3 computeFc(unsigned nbContacts, stateObservation::Vector x, stateObservation::Vector u, stateObservation::flexibilityEstimation::ModelBaseEKFFlexEstimatorIMU&)
{
    stateObservation::Vector3 Fc;
    stateObservation::Vector3 oriV, angVel, angAcc, pos, vel, acc, fm, tm, AngMomentum, dotAngMomentum;
    stateObservation::Matrix3 rFlex, inertia, dotInertia;
    stateObservation::Vector3 cl, dcl, ddcl;
    IndexedMatrixArray contactPosV;
    IndexedMatrixArray contactOriV;
    IndexedMatrixArray efforts;
    Vector fc_;
    Vector tc_;
    stateObservation::Vector3 gmuz;

    gmuz << 0,
            0,
            9.81*56.8;
    fm.setZero(); tm.setZero();

    oriV = x.segment(stateObservation::flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::state::ori,3);
    angVel = x.segment(stateObservation::flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::state::angVel,3);
    rFlex = kine::rotationVectorToRotationMatrix(oriV);

    cl = u.segment(0,3);
    dcl = u.segment(3,3);
    ddcl = u.segment(6,3);
    AngMomentum=u.segment<3>(stateObservation::flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::input::angMoment);
    dotAngMomentum=u.segment<3>(stateObservation::flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::input::dotAngMoment);
    buildInertiaTensor(u.segment<6>(stateObservation::flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::input::inertia), inertia);
    buildInertiaTensor(u.segment<6>(stateObservation::flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::input::dotInertia), dotInertia);

    for (unsigned i = 0; i<nbContacts ; ++i)
    {
        contactPosV.setValue(u.segment<3>(stateObservation::flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::input::contacts + 12*i),i);
        contactOriV.setValue(u.segment<3>(stateObservation::flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::input::contacts +12*i+3),i);
    }

    for (unsigned i=0; i<hrp2::contact::nbModeledMax; ++i)
    {
        efforts[i]=x.segment(stateObservation::flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::state::fc+6*i,6);
        fc_.segment<3>(3*i) = efforts[i].block<3,1>(0,0);
        tc_.segment<3>(3*i) = efforts[i].block<3,1>(3,0);
    }

//    est1.getEKF().getFunctor()->computeAccelerations(cl, dcl, ddcl, AngMomentum, dotAngMomentum, inertia, dotInertia,  contactPosV,
//                                                    contactOriV, pos, vel, acc, oriV, rFlex, angVel, angAcc, fc_, tc_, fm, tm);

    Fc = gmuz + 56.8*(2*kine::skewSymmetric(angVel)*rFlex*dcl+rFlex*ddcl);

    return Fc;
}

int test()
{
    std::cout << "Starting" << std::endl;

    // For measurement vector
    bool withUnmodeledForces_ = false;
    bool withForceSensors_=true;
    bool withAbsolutePose_ = false;

    // For state vector
    bool withComBias_=false;

    // Time
    const double dt=5e-3;
    const unsigned kinit=10000; //3500; //3758; //
    const unsigned kmax=11500; //5700; //5100; //

    // Fix sizes
    const unsigned measurementSizeBase=6;
    const unsigned inputSizeBase=42;
    const unsigned stateSize=35;

    /// Definitions of input vectors
     // Measurement
     IndexedMatrixArray y;
     std::cout << "Loading measurements file" << std::endl;
     y.readFromFile("source_measurement.dat",1,30);

     // Input
     IndexedMatrixArray u;
     std::cout << "Loading input file" << std::endl;
     u.readFromFile("source_input.dat",1,66);

     // Number of support contacts
     IndexedMatrixArray nbSupport;
     std::cout << "Loading the number of supports file" << std::endl;
     nbSupport.readFromFile("source_nbSupport.dat",1,1);

     // CoM bias
     IndexedMatrixArray bias;
     std::cout << "Loading comBias file" << std::endl;
     bias.readFromFile("source_comBias.dat",1,35);

     // Zmp ref
     IndexedMatrixArray zmpRef;
     std::cout << "Loading zmpRef file" << std::endl;
     zmpRef.readFromFile("source_zmpRef.dat",1,3);

    /// Definition of ouptut vectors
     // State: what we want
     IndexedMatrixArray x_output;
     // State: what we want
     IndexedMatrixArray xPredicted_output;
     // Measurement
     IndexedMatrixArray y_output;
     // Input
     IndexedMatrixArray u_output;

     // Zmp from sensors
     IndexedMatrixArray sensorZmp_output;
     // Zmp from filter
     IndexedMatrixArray filteredZmp_output;
     // reference Zmp
     IndexedMatrixArray referenceZmp_output;

    std::cout << "Creating estimator" <<std::endl;

    stateObservation::flexibilityEstimation::ModelBaseEKFFlexEstimatorIMU est;
    est.setSamplingPeriod(dt);

    // Model
    est.setContactModel(1);
    est.setRobotMass(56.8);
    est.setKfe(40000*Matrix3::Identity());
    est.setKte(350*Matrix3::Identity());
    est.setKfv(600*Matrix3::Identity());
    est.setKtv(10*Matrix3::Identity());

    // Config
    est.setWithUnmodeledForces(withUnmodeledForces_);
    est.setWithForcesMeasurements(withForceSensors_);
    est.setWithAbsolutePos(withAbsolutePose_);
    est.setWithComBias(withComBias_);

    std::cout << "Setting covariances" << std::endl;

    // Measurement noise covariance
    stateObservation::Matrix R_; R_.resize(measurementSizeBase,measurementSizeBase); R_.setIdentity();
    R_.block(0,0,3,3)*=1.e-3;
    R_.block(3,3,3,3)*=1.e-6;
    est.setMeasurementNoiseCovariance(R_);
    est.setUnmodeledForceVariance(1e-13);
    est.setForceVariance(1.e-6);
    est.setAbsolutePosVariance(1e-4);

    // Process noise covariance
    stateObservation::Matrix Q_; Q_.resize(stateSize,stateSize); Q_.setIdentity();
    Q_.block(0,0,12,12)*=1.e-8;
    Q_.block(12,12,12,12)*=1.e-4;
    Q_.block(24,24,6,6)*=1.e-2;
    Q_.block(30,30,2,2)*=1.e-7;
    Q_.block(32,32,3,3)*=1.e-8;
    est.setProcessNoiseCovariance(Q_);

    // Temporary variables
    stateObservation::Vector input, measurement;
    Vector x; x.resize(stateSize);
    Vector xPredicted; xPredicted.resize(stateSize);
    unsigned contactNbr;
    size_t inputSize, measurementSize;

    stateObservation::IndexedMatrixArray measurementForces, filteredForces, inputFeetPositions;
    stateObservation::Vector sensorZmp, filteredZmp;
    stateObservation::Vector6 filteredForcesk, measurementForcesk;
    stateObservation::Matrix4 inputFeetPositionsk;

    std::cout << "Beginning reconstruction "<<std::endl;

    for (unsigned k=kinit;k<kmax;++k)
    {
        std::cout << "\n" << k << std::endl;

        if(nbSupport[k](0)!=est.getContactsNumber())
        {
            contactNbr = unsigned(nbSupport[k](0));
            est.setContactsNumber(contactNbr);
        }

        measurementSize = est.getMeasurementSize();
        measurement.resize(measurementSize);
        measurement = (y[k].block(0,0,1,measurementSize)).transpose();
        est.setMeasurement(measurement);

        inputSize = est.getInputSize();
        input.resize(inputSize);
        input = (u[k+1].block(0,0,1,inputSize)).transpose();
        est.setMeasurementInput(input);

        xPredicted = est.getEKF().getFunctor()->stateDynamics(x,(u[k-1].block(0,0,1,inputSize)).transpose(),0);
        xPredicted.segment(stateObservation::flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::state::fc,12)
               = x.segment(stateObservation::flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::state::fc,12);

        x = est.getFlexibilityVector();

        // Compute Zmp from sensor and filter
        for(unsigned i=0;i<contactNbr;++i)
        {
            measurementForcesk = y[k].block<1,6>(0,measurementSizeBase+withUnmodeledForces_*6+i*6).transpose();
            filteredForcesk = x.segment<6>(stateObservation::flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::state::fc+i*6);
            inputFeetPositionsk = kine::vector6ToHomogeneousMatrix(u[k].block<1,6>(0,inputSizeBase+i*12)).transpose();

//            measurementForcesk.segment<3>(0) = inputFeetPositionsk.block(0,0,3,3).transpose()*measurementForcesk.segment<3>(0);
//            measurementForcesk.segment<3>(3) = inputFeetPositionsk.block(0,0,3,3).transpose()*measurementForcesk.segment<3>(3)
//                                               +kine::skewSymmetric(inputFeetPositionsk.block(0,0,3,3).transpose()*inputFeetPositionsk.block(0,3,3,1))*measurementForcesk.segment<3>(0);
//
//            filteredForcesk.segment<3>(0) = inputFeetPositionsk.block(0,0,3,3).transpose()*filteredForcesk.segment<3>(0);
//            filteredForcesk.segment<3>(3) = inputFeetPositionsk.block(0,0,3,3).transpose()*filteredForcesk.segment<3>(3)
//                                               -kine::skewSymmetric(inputFeetPositionsk.block(0,0,3,3).transpose()*inputFeetPositionsk.block(0,3,3,1))*filteredForcesk.segment<3>(0);

            measurementForces.setValue(measurementForcesk,i);
            filteredForces.setValue(filteredForcesk,i);
            stateObservation::Matrix4 identity(stateObservation::Matrix4::Identity());
            inputFeetPositions.setValue(identity,TimeIndex(i));
        }
        sensorZmp=computeZmp(contactNbr, measurementForces, inputFeetPositions);
        filteredZmp=computeZmp(contactNbr, filteredForces, inputFeetPositions);

        x_output.setValue(x,k);
        xPredicted_output.setValue(xPredicted,k);
        y_output.setValue(y[k],k);
        u_output.setValue(u[k],k);
        sensorZmp_output.setValue(sensorZmp,k);
        filteredZmp_output.setValue(filteredZmp,k);
        referenceZmp_output.setValue(zmpRef[k],k);
    }

    std::cout << "Completed "<<std::endl;

    x_output.writeInFile("state.dat");
    xPredicted_output.writeInFile("statePredicted.dat");
    y_output.writeInFile("measurement.dat");
    u_output.writeInFile("input.dat");
    sensorZmp_output.writeInFile("sensorZmp.dat");
    filteredZmp_output.writeInFile("filteredZmp.dat");
    referenceZmp_output.writeInFile("referenceZmp.dat");

    return 1;
}

int main()
{
    return test();
}
