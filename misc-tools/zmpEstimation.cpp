#include <iostream>
#include <fstream>

//#include <state-observation/noise/gaussian-white-noise.hpp>
//#include <state-observation/examples/offline-ekf-flexibility-estimation.hpp>
//#include <state-observation/dynamical-system/dynamical-system-simulator.hpp>
//#include <state-observation/tools/miscellaneous-algorithms.hpp>

#include <state-observation/flexibility-estimation/model-base-ekf-flex-estimator-imu.hpp>

#include <time.h>

using namespace stateObservation;

timespec diff(const timespec & start, const timespec & end)
{
        timespec temp;
        if ((end.tv_nsec-start.tv_nsec)<0) {
                temp.tv_sec = end.tv_sec-start.tv_sec-1;
                temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
        } else {
                temp.tv_sec = end.tv_sec-start.tv_sec;
                temp.tv_nsec = end.tv_nsec-start.tv_nsec;
        }
        return temp;
}

Vector computeZmp (Vector& forces, Matrix& sensorPosition1, Matrix& sensorPosition2, unsigned contactNbr)
{
    double fnormal = 0;
    double sumZmpx = 0;
    double sumZmpy = 0;
    double sumZmpz = 0;
    Vector zmp;
    zmp.resize (3);

    for (unsigned int i=0; i<contactNbr; ++i) {
      const Vector& f = forces.block(i*6,0,6,1);
      // Check that force is of dimension 6
      if (f.size () != 6) {
        zmp.fill (0.);
        return zmp;
      }

      Matrix M;
      M.resize(4,4);
      if (i==0)
      {
        M = sensorPosition1;

      }
      if (i==1)
      {
        M = sensorPosition2;

      }

     double fx = M (0,0) * f (0) + M(0,1) * f (1) + M (0,2) * f (2);
     double fy = M (1,0) * f (0) + M(1,1) * f (1) + M (1,2) * f (2);
     double fz = M (2,0) * f (0) + M(2,1) * f (1) + M (2,2) * f (2);

      if (fz > 0) {
    double Mx = M (0,0)*f(3) + M (0,1)*f(4) + M (0,2)*f(5)
                + M (1,3)*fz - M (2,3)*fy;
        double My = M (1,0)*f(3) + M (1,1)*f(4) + M (1,2)*f(5)
                + M (2,3)*fx - M (0,3)*fz;
        fnormal += fz;
        sumZmpx -= My;
        sumZmpy += Mx;
        sumZmpz += fz * M (2,3);
      }
    }
    if (fnormal != 0) {
      zmp (0) = sumZmpx / fnormal;
      zmp (1) = sumZmpy / fnormal;
      zmp (2) = sumZmpz / fnormal;
    } else {
        zmp.fill (0.);
    }
    return zmp;
}

int test()
{
    std::cout << "Starting" << std::endl;

  /// sampling period
    const double dt=5e-3;

    stateObservation::flexibilityEstimation::ModelBaseEKFFlexEstimatorIMU est;

  /// Measurement noise covariance
    Matrix Cov;
    Cov.resize(6,6);
    double unitCov1 = 1e6;
    double unitCov2 = 1e6;
    Cov <<  unitCov1,0,0,0,0,0,
            0,unitCov1,0,0,0,0,
            0,0,unitCov1,0,0,0,
            0,0,0,unitCov2,0,0,
            0,0,0,0,unitCov2,0,
            0,0,0,0,0,unitCov2;

  /// Initializations
    // Dimensions
    const unsigned kinit=0;
    const unsigned kmax=42761;
    const unsigned measurementSize=6;

    const unsigned stateSize=18;
    unsigned contactNbr = 2;
    // State initialization => not used here because it is set in model-base-ekf-flex-estimator-imu


    const unsigned inputSize=est.getInputSize();

     // Input initialization
    Vector u0=Vector::Zero(inputSize-6*contactNbr,1);
    u0 <<   0.0135672,
            0.001536,
            0.80771,
            -2.50425e-06,
            -1.03787e-08,
            5.4317e-08,
            -2.50434e-06,
            -1.03944e-08,
            5.45321e-08,
            48.1348,
            46.9498,
            1.76068,
            -0.0863332,
            -0.59487,
            -0.0402246,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            -0.098,
            -1.21619e-10,
            1.1174,
            3.06752e-22,
            -1.06094e-20,
            7.75345e-22,
            -2.84609e-06,
            -1.18496e-08,
            -4.52691e-18,
            2.95535e-20,
            -1.0346e-18,
            7.58731e-20,
            -0.000284609,
            -1.18496e-06,
            -4.52691e-16;


    est.setSamplingPeriod(dt);

    est.setInput(u0);
    est.setMeasurementInput(u0);


    est.setKfe(40000*Matrix3::Identity());
    est.setKfv(600*Matrix3::Identity());
    est.setKte(400*Matrix3::Identity());
    est.setKtv(10*Matrix3::Identity());

   /// Definitions of input vectors
     // Measurement
     IndexedMatrixArray y;
     std::cout << "Loading measurements file" << std::endl;
     y.readFromFile("data/source_measurement.dat",1,measurementSize);
     // Input
     IndexedMatrixArray u;
     std::cout << "Loading input file" << std::endl;
     u.readFromFile("data/source_input.dat",1,inputSize);
     //state
     IndexedMatrixArray xRef;
     std::cout << "Loading reference state file" << std::endl;
     xRef.readFromFile("data/source_state.dat",stateSize,1);


     //zmp estimated ref
     IndexedMatrixArray zmpEstimatedRef;
     std::cout << "Loading reference zmpEstimated file" << std::endl;
     zmpEstimatedRef.readFromFile("data/source_zmpestimated.dat",3,1);
     //zmp
     IndexedMatrixArray zmpRef;
     std::cout << "Loading reference zmp file" << std::endl;
     zmpRef.readFromFile("data/source_zmp.dat",3,1);
     //forcesAndMoments
     IndexedMatrixArray forcesAndMomentsRef;
     std::cout << "Loading reference forcesAndMoments file" << std::endl;
     forcesAndMomentsRef.readFromFile("data/source_forcesAndMoments.dat",12,1);
     //forceLLEG
     IndexedMatrixArray forceLLEGRef;
     std::cout << "Loading reference forceLLEG file" << std::endl;
     forceLLEGRef.readFromFile("data/source_forceLLEG.dat",6,1);
     //forceRLEG
     IndexedMatrixArray forceRLEGRef;
     std::cout << "Loading reference forceRLEG file" << std::endl;
     forceRLEGRef.readFromFile("data/source_forceRLEG.dat",6,1);
     //forceLLEG estimated
     IndexedMatrixArray forceLLEGEst;
     std::cout << "Loading estimated forceLLEG file" << std::endl;
     forceLLEGEst.readFromFile("data/source_forcesSupport1.dat",6,1);
     //forceRLEG estimated
     IndexedMatrixArray forceRLEGEst;
     std::cout << "Loading estimated forceRLEG file" << std::endl;
     forceRLEGEst.readFromFile("data/source_forcesSupport2.dat",6,1);

   /// Definition of ouptut vectors
     // State: what we want
     IndexedMatrixArray x_output;
     // Measurement
     IndexedMatrixArray y_output;
     // Input
     IndexedMatrixArray u_output;

     // ZMP estimated
     IndexedMatrixArray zmpEstimated_output;
     IndexedMatrixArray zmpReEstimated_output;
     // ZMP
     IndexedMatrixArray zmp_output;
     // force LLEG
     IndexedMatrixArray forceLLEG_output;
     // force RLEG
     IndexedMatrixArray forceRLEG_output;
     // forcesAndMoments
     IndexedMatrixArray forcesAndMoments_output;

     IndexedMatrixArray deltax_output;

     est.setMeasurementNoiseCovariance(Cov);
     est.setContactsNumber(contactNbr);

    Vector flexibility;
    flexibility.resize(18);
    Vector xdifference(flexibility);

    Vector forcesAndMoments_source;
    forcesAndMoments_source.resize(12);

    Vector forcesAndMomentsReference;
    forcesAndMomentsReference.resize(12);

    stateObservation::Vector forcesAndMoments;
    forcesAndMoments.resize(12);
    Vector forcesAndMomentsReal;
    forcesAndMomentsReal.resize(12);
    Vector zmpEstimated, zmpReEstimated;
    zmpEstimated.resize(3), zmpReEstimated.resize(3);
    Vector zmp;
    zmp.resize(3);

    est.setContactModel(stateObservation::flexibilityEstimation::
                ModelBaseEKFFlexEstimatorIMU::contactModel::elasticContact);

    Matrix sensorPosition1, sensorPosition2, sensorPosition1Est, sensorPosition2Est;
    sensorPosition1.resize(4,4);
    sensorPosition2.resize(4,4);
    sensorPosition1Est.resize(4,4);
    sensorPosition2Est.resize(4,4);

    std::cout << "Beginning reconstruction "<<std::endl;
    for (unsigned k=kinit+2;k<kmax;++k)
    {
        est.setMeasurement(y[k].transpose());
        est.setMeasurementInput(u[k].transpose());

        flexibility = est.getFlexibilityVector();

        x_output.setValue(flexibility,k);
        y_output.setValue(y[k],k);
        u_output.setValue(u[k],k);

        forcesAndMomentsReference=forcesAndMomentsRef[k];
        forcesAndMoments = est.getForcesAndMoments();
//        forcesAndMoments.block(0,0,6,1)=forceRLEGEst[k];
//        forcesAndMoments.block(6,0,6,1)=forceLLEGEst[k];

        // Compute forces and ZMP estimated in this code
        sensorPosition1=kine::vector6ToHomogeneousMatrix((u[k].block(0,42,1,6)).transpose());
        sensorPosition2=kine::vector6ToHomogeneousMatrix((u[k].block(0,48,1,6)).transpose());
        zmpReEstimated=computeZmp (forcesAndMoments, sensorPosition1, sensorPosition2, contactNbr);

//        // Compute forces and ZMP estimated
//        sensorPosition1Est=kine::vector6ToHomogeneousMatrix((u[k].block(0,42,1,6)).transpose());
//        sensorPosition2Est=kine::vector6ToHomogeneousMatrix((u[k].block(0,48,1,6)).transpose());
//        sensorPosition1Est(0,3)=0;
//        sensorPosition2Est(0,3)=0;
//        zmpEstimated=computeZmp (forcesAndMoments, sensorPosition1Est, sensorPosition2Est, contactNbr);

//        // Compute ZMP real
//        forcesAndMoments_source.block(0,0,6,1)=forceRLEGRef[k];
//        forcesAndMoments_source.block(6,0,6,1)=forceLLEGRef[k];
//        zmp=computeZmp (forcesAndMoments_source, sensorPosition1, sensorPosition2, contactNbr);

        if(zmpEstimated(0)>100)
        {
            zmpEstimated(0)=100;
        }
        if(zmpEstimated(1)>100)
        {
            zmpEstimated(1)=100;
        }
        if(zmpEstimated(2)>100)
        {
            zmpEstimated(2)=100;
        }
        if(zmpEstimated(0)<-100)
        {
            zmpEstimated(0)=-100;
        }
        if(zmpEstimated(1)<-100)
        {
            zmpEstimated(1)=-100;
        }
        if(zmpEstimated(2)<-100)
        {
            zmpEstimated(2)=-100;
        }

        forcesAndMoments_output.setValue(forcesAndMoments,k);
        zmpReEstimated_output.setValue(zmpReEstimated,k);
        //zmp_output.setValue(zmp,k);

    }

    std::cout << "Completed "<<std::endl;

    x_output.writeInFile("state.dat");
    std::cout << "coucou" << std::endl;
    y_output.writeInFile("measurement.dat");
    std::cout << "coucou" << std::endl;
    u_output.writeInFile("input.dat");
    std::cout << "coucou" << std::endl;
    forcesAndMoments_output.writeInFile("forcesAndMoments.dat");
    std::cout << "coucou" << std::endl;
    zmpReEstimated_output.writeInFile("zmpReEstimated.dat");
    std::cout << "coucou" << std::endl;
    //zmp_output.writeInFile("zmp.dat");

    return 0;

}

int main()
{

    return test();

}






