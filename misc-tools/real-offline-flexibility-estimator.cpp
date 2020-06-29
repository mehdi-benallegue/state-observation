#include <iostream>
#include <fstream>

#include <state-observation/noise/gaussian-white-noise.hpp>
#include <state-observation/examples/offline-ekf-flexibility-estimation.hpp>
#include <state-observation/dynamical-system/dynamical-system-simulator.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>
#include <state-observation/dynamical-system/imu-dynamical-system.hpp>


using namespace stateObservation;

typedef kine::indexes<kine::rotationVector> indexes;

///sampling period
const double dt=5e-3;

int testDerivator()
{
    /// The number of samples
    const unsigned kmax=3000;

    ///sampling period
    const double dt=1e-3;

    ///Sizes of the states for the state, the measurement, and the input vector
    const unsigned stateSize=18;
    const unsigned measurementSize=6;
    //const unsigned inputSize=6;

    ///The array containing all the states, the measurements and the inputs
    IndexedVectorArray x;
    IndexedVectorArray y;
    IndexedVectorArray u;

    ///The covariance matrix of the process noise and the measurement noise
    Matrix q;
    Matrix r;

    {
        ///simulation of the signal
        /// the IMU dynamical system functor
        IMUDynamicalSystem imu;

        ///The process noise initialization
        Matrix q1=Matrix::Identity(stateSize,stateSize)*0.00;
        GaussianWhiteNoise processNoise(imu.getStateSize());
        processNoise.setStandardDeviation(q1);
        imu.setProcessNoise( & processNoise );
        q=q1*q1.transpose();

        ///The measurement noise initialization
        Matrix r1=Matrix::Identity(measurementSize,measurementSize)*0.0;
        GaussianWhiteNoise MeasurementNoise(imu.getMeasurementSize());
        MeasurementNoise.setStandardDeviation(r1);
        imu.setMeasurementNoise( & MeasurementNoise );
        r=r1*r1.transpose();

        ///the simulator initalization
        DynamicalSystemSimulator sim;
        sim.setDynamicsFunctor(&imu);

        ///initialization of the state vector
        Vector x0=Vector::Zero(stateSize,1);
        sim.setState(x0,0);

        ///construction of the input
        /// the input is constant over 10 time samples
        for (unsigned i=0;i<kmax/10;++i)
        {
            Vector uk=Vector::Zero(imu.getInputSize(),1);

            uk[0]=0.4 * sin(M_PI/10*i);
            uk[1]=0.6 * sin(M_PI/12*i);
            uk[2]=0.2 * sin(M_PI/5*i);

            uk[3]=10  * sin(M_PI/12*i);
            uk[4]=0.07  * sin(M_PI/15*i);
            uk[5]=0.05 * sin(M_PI/5*i);

            ///filling the 10 time samples of the constant input
            for (int j=0;j<10;++j)
            {
                u.setValue(uk,i*10+j);
            }

            ///give the input to the simulator
            ///we only need to give one value and the
            ///simulator takes automatically the appropriate value
            sim.setInput(uk,10*i);

        }

        ///set the sampling perdiod to the functor
        imu.setSamplingPeriod(dt);

        ///launched the simulation to the time kmax+1
        sim.simulateDynamicsTo(kmax+1);

        ///extract the array of measurements and states
        y = sim.getMeasurementArray(1,kmax);
        x = sim.getStateArray(1,kmax);
    }

    IndexedVectorArray dta;

    for (TimeIndex i=x.getFirstIndex() ; i<x.getNextIndex() ; ++i)
    {
        Vector xi = Vector::Zero(6,1);

        xi.head(3) = Vector(x[i]).segment(indexes::pos,3);
        xi.tail(3) = Vector(x[i]).segment(indexes::ori,3);

        dta.setValue(xi,i);
    }

    IndexedVectorArray state = kine::reconstructStateTrajectory(dta,dt);

    ///file of output
    std::ofstream f;
    f.open("trajectory.dat");

    for (TimeIndex i=x.getFirstIndex() ; i<x.getNextIndex() ; ++i)
    {
        //f<<dta[i].transpose()<<"\t#####\t\t"<<Vector(Vector(x[i]).segment(9,3)).transpose()<<"\t#####\t\t"<<(x[i]-state[i]).transpose()<<std::endl;
        f<<(x[i]-state[i])<<std::endl<<std::endl;
    }

    return 0;


}

IndexedVectorArray getMeasurements(const char * accelerometerSignal, const  char * gyrometerSignal)
{
    std::ifstream facc;
    std::ifstream fgyr;

    facc.open(accelerometerSignal);
    fgyr.open(gyrometerSignal);

    Vector3 mAcc;
    Vector3 mGyr;

    Vector yk=Vector::Zero(6,1);

    IndexedVectorArray y;

    bool continuation=true;

    while (continuation)
    {
        unsigned k1;
        unsigned k2;
        facc >> k1;
        fgyr >> k2;
        if (facc.eof()||facc.eof()||k1!=k2)
            continuation=false;

        if (continuation)
        {
            facc >> mAcc[0]>> mAcc[1]>> mAcc[2];
            fgyr >> mGyr[0]>> mGyr[1]>> mGyr[2];

            yk.head(3)=mAcc;
            yk.tail(3)=mGyr;
            y.setValue(yk,k1);
        }
    }

    return y;
}

IndexedVectorArray getTrajectory(const char * PositionOrientation)
{
    std::ifstream f;

    f.open(PositionOrientation);

    Vector6 configuration;

    IndexedVectorArray up;

    bool continuation=true;

    while (continuation)
    {
        unsigned k;

        f >> k;

        if (f.eof())
            continuation=false;

        if (continuation)
        {
            f >> configuration[0]>> configuration[1]>> configuration[2]
              >> configuration[3]>> configuration[4]>> configuration[5];

            up.setValue(configuration,k);
        }
    }

    IndexedVectorArray state = kine::reconstructStateTrajectory(up,dt);

    return state;
}



int test (const IndexedVectorArray & y, const IndexedVectorArray & u)
{
    /// The number of samples
    const unsigned stateSize = 18;



    ///the initalization of an estimation of the initial state
    Vector xh0=Vector::Zero(stateSize,1);

    std::vector<Vector3, Eigen::aligned_allocator<Vector3> > contactPositions;

    contactPositions.push_back(Matrix::Zero(3,0));

    stateObservation::IndexedVectorArray xh=
        stateObservation::examples::offlineEKFFlexibilityEstimation
        (y,u,xh0,1,contactPositions,dt);

    ///file of output
    std::ofstream f;
    f.open("trajectory.dat");

    double error;

    ///the reconstruction of the state
    for (TimeIndex i=y.getFirstIndex();i<y.getNextIndex();++i)
    {
        ///display part, useless
        Vector3 g;
        {
            g = Vector(y[i]).head(3);
            g.normalize();
        }

        Vector3 gh;
        {
            Matrix3 Rh;

            Vector3 orientationV=Vector(xh[i]).segment(indexes::ori,3);
            double angle=orientationV.norm();
            if (angle > cst::epsilonAngle)
                Rh = AngleAxis(angle, orientationV/angle).toRotationMatrix();
            else
                Rh = Matrix3::Identity();
            gh=Rh.transpose()*Vector3::UnitZ();
            gh.normalize();
        }

        error = acos(double(g.transpose()*gh)) * 180 / M_PI;


        f << i<< " \t "<< error << " \t\t\t "
          << g.transpose() << " \t\t\t " << gh.transpose() << " \t\t\t "
          << xh[i].transpose() << std::endl;
    }

    std::cout << "Error " << error << ", test: " ;

    if (error > 2)
    {
        std::cout << "FAILED !!!!!!!";
        return 1;
    }
    else
    {
        std::cout << "SUCCEEDED !!!!!!!";
        return 0;
    }
}

int main()
{

    IndexedVectorArray y =
            getMeasurements("dg_HRP2LAAS-accelerometer.dat",
                "dg_HRP2LAAS-gyrometer.dat");

    IndexedVectorArray u= getTrajectory("IMUTrajectory.dat");



    return test(y,u);

}
