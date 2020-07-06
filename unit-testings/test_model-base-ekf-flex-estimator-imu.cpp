#include <iostream>
#include <fstream>

//#include <state-observation/noise/gaussian-white-noise.hpp>
//#include <state-observation/examples/offline-ekf-flexibility-estimation.hpp>
//#include <state-observation/dynamical-system/dynamical-system-simulator.hpp>
//#include <state-observation/tools/miscellaneous-algorithms.hpp>

#include <state-observation/flexibility-estimation/model-base-ekf-flex-estimator-imu.hpp>
#include <state-observation/config.h>
#include <time.h>


using namespace stateObservation;

typedef kine::indexes<kine::rotationVector> indexes;

int test()
{
    std::cout << "Starting" << std::endl;

  /// sampling period
    const double dt=5e-3;

  /// Measurement noise covariance
    Matrix Cov;
    Cov.resize(6,6);
    double unitCov = 1e-2;
    Cov <<  unitCov,0,0,0,0,0,
            0,unitCov,0,0,0,0,
            0,0,unitCov,0,0,0,
            0,0,0,unitCov,0,0,
            0,0,0,0,unitCov,0,
            0,0,0,0,0,unitCov;

  /// Initializations
    // Dimensions
    const size_t kinit=0;
    const size_t kmax=1400;
    const unsigned measurementSize=6;
    (void) measurementSize; ///avoid warning
    const unsigned inputSize=60;
    const unsigned stateSize=18;
    (void) stateSize;
    unsigned contactNbr = 2;
    // State initialization => not used here because it is set in model-base-ekf-flex-estimator-imu

     // Input initialization
    Vector u0=Vector::Zero(inputSize-6*contactNbr,1);
    u0 <<  0.0135673,
             0.001536,
             0.80771,
             -2.63605e-06,
             -1.09258e-08,
             5.71759e-08,
             2.71345,
             0.3072,
             161.542,
             48.1348,
             46.9498,
             1.76068,
             -0.0863332,
             -0.594871,
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
             -6.23712e-11,
             1.1174,
             1.58984e-22,
             -5.43636e-21,
             3.9598e-22,
             -2.99589e-06,
             -1.24742e-08,
             -4.7647e-18,
             3.17968e-20,
             -1.08727e-18,
             7.91959e-20,
             -0.000299589,
             -1.24742e-06,
             -4.7647e-16,
             0,
             0,
             0,
             0,
             0,
             0; ///external forces and moments

    stateObservation::flexibilityEstimation::ModelBaseEKFFlexEstimatorIMU est;
    est.setSamplingPeriod(dt);
    est.setInput(u0);
    est.setMeasurementInput(u0);

    est.setRobotMass(56.8679920);
    stateObservation::Vector3 v1; v1.setOnes();
    v1=1000*v1;
    est.setTorquesLimit(v1);
    stateObservation::Vector3 v2; v2.setOnes();
    v2=1000*v2;
    est.setForcesLimit(v2);
    est.setKfe(40000*Matrix3::Identity());
    est.setKte(600*Matrix3::Identity());
    est.setKfv(600*Matrix3::Identity());
    est.setKtv(60*Matrix3::Identity());

   /// Definitions of input vectors
     // Measurement
    IndexedVectorArray y;
    std::cout << "Loading measurements file" << std::endl;
    y.readVectorsFromFile("inputFiles/source_measurement.dat");
    std::cout << "Done, size " << y.size()<<  std::endl;
     // Input
    IndexedVectorArray u;
     std::cout << "Loading input file" << std::endl;
    u.readVectorsFromFile("inputFiles/source_input.dat");
    std::cout << "Done, size " << u.size()<<  std::endl;
      //state
    IndexedVectorArray xRef;
    std::cout << "Loading reference state file" << std::endl;
    xRef.readVectorsFromFile("inputFiles/source_state.dat");
    std::cout << "Done, size " << xRef.size()<<  std::endl;


   /// Definition of ouptut vectors
     // State: what we want
    IndexedVectorArray x_output;
     // Measurement
    IndexedVectorArray y_output;
     // Input
    IndexedVectorArray u_output;

    IndexedVectorArray deltax_output;


    est.setMeasurementNoiseCovariance(Cov);

    est.setContactsNumber(contactNbr);


    Vector flexibility;
    flexibility.resize(18);
    Vector xdifference(flexibility);

    tools::SimplestStopwatch stopwatch;
    IndexedVectorArray computationTime_output;
    double computationTime_sum=0;
    Vector computeTime;
    computeTime.resize(1);

    Vector errorsum=Vector::Zero(12);

    est.setContactModel(stateObservation::flexibilityEstimation::
                ModelBaseEKFFlexEstimatorIMU::contactModel::elasticContact);

    // for (TimeIndex k=u.getFirstIndex(); k<u.getNextIndex();++k)
    // {
    //   Vector uk(est.getInputSize());
    //   uk.head<42>() = u[k].head<42>();
    //   uk.segment<6>(42) << 0,0,0,0,0,0;
    //   for (unsigned  i =0; i<contactNbr ;++i)
    //   {
    //     uk.segment<6>(42+12*i)=u[k].segment<6>(42+6*i);
    //     uk.segment<6>(42+12*i+6)<<0,0,0,0,0,0;
    //   }

    //   u[k]=uk;
    // }

    std::cout << "Beginning reconstruction "<<std::endl;
    for (size_t k=kinit+2;k<kmax;++k)
    {
        est.setMeasurement(y[k]);
        est.setMeasurementInput(u[k]);

        stopwatch.start();

        flexibility = est.getFlexibilityVector();

        computeTime[0]=stopwatch.stop();

        xdifference =flexibility.head<12>()-xRef[k].head<12>();

        x_output.setValue(flexibility,k);
        y_output.setValue(y[k],k);
        u_output.setValue(u[k],k);
        deltax_output.setValue(xdifference,k);

        computationTime_output.setValue(computeTime,k);
        computationTime_sum+=computeTime[0];
    }

    std::cout << "Completed "<<std::endl;

    computeTime[0]=computationTime_sum/(kmax-kinit-2);
    computationTime_output.setValue(computeTime,kmax);
    computationTime_output.writeInFile("computationTime.dat");

    /// this code is useful to generate new input files
    // x_output.writeInFile("unit-testings/result-state.dat");
    // y_output.writeInFile("unit-testings/result-measurement.dat");
    // u_output.writeInFile("unit-testings/result-input.dat");

    errorsum = errorsum/(kmax-kinit-2);

    Vector error(4);

    typedef flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::state State;

    error(0) = sqrt((errorsum(State::pos) + errorsum(State::pos+1) + errorsum(State::pos+2)));
    error(1) = sqrt((errorsum(State::linVel) + errorsum(State::linVel+1) + errorsum(State::linVel+2)));
    error(2) = sqrt((errorsum(State::ori) + errorsum(State::ori+1) + errorsum(State::ori+2)));
    error(3) = sqrt((errorsum(State::angVel) + errorsum(State::angVel+1) + errorsum(State::angVel+2)));

    std::cout << "Mean computation time " << computeTime[0] <<std::endl;

    std::cout << "Mean error " << error.transpose() <<std::endl;

    double posgain = 40000;
    double linvelgain = 600;

    double origain = 10;
    double angvelgain = 1;


    double syntherror = posgain*error(0)+ linvelgain*error(1) + origain*error(2) + angvelgain*error(3);

    std::cout << "Synthetic error " << posgain*error(0) << " + " << linvelgain*error(1)
                          << " + " << origain*error(2) << " + " << angvelgain*error(3)
                          << " = " << syntherror << std::endl;


    if (syntherror == syntherror)
    {
      if (syntherror>0.2)
      {
        std::cout << "Failed : error is too big !!"<< std::endl <<"The end" << std::endl;
        return 1;
      }
    }
    else
    {
      if (syntherror>0.2)
      {
        std::cout << "Failed : NaN !!"<< std::endl <<"The end" << std::endl;
        return 2;
      }
    }
    
    if (Release_BUILD)
    {
        if (computeTime[0] > 2e5)
        {
            std::cout << "Failed : Computation time is too long !!" << std::endl << "The end" << std::endl;
            return 2;
        }
    }

    std::cout << "Succeed !!"<< std::endl <<"The end" << std::endl;
    return 0;
}

int main()
{

    return test();

}
