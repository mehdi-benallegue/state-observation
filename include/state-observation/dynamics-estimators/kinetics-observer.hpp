/**
 * \file      kinetics-observer-base.hpp
 * \author    Mehdi Benallegue
 * \date      2018
 * \brief     Unified Kinetics estimator
 *
 * \details
 *
 *
 */

#ifndef FLEXIBILITYESTIMATION_EKFFLEXIBILITYESTIMATORBASE_H
#define FLEXIBILITYESTIMATION_EKFFLEXIBILITYESTIMATORBASE_H

#include <map>
#include <set>

#include <boost/utility.hpp>

#include <state-observation/tools/definitions.hpp>
#include <state-observation/observer/extended-kalman-filter.hpp>
#include <state-observation/flexibility-estimation/flexibility-estimator-base.hpp>
#include <state-observation/tools/rigid-body-kinematics.hpp>
#include <state-observation/dynamical-system/dynamical-system-functor-base.hpp>
#include <state-observation/sensors-simulation/accelerometer-gyrometer.hpp>
#include <state-observation/noise/noise-base.hpp>


///////////SIZE OF VECTORS
static const int stateSizeBase = 3+4+3+3+3+3+3;
static const int stateSizePerContact = 3+4+3+3;
static const int stateTangentSizeBase = 3+3+3+3+3+3+3;
static const int stateTangentSizePerContact = 3+3+3+3;

static const int sizeAcceleroSignal = 3;
static const int sizeGyroSignal = 3;
static const int sizeIMUSignal = sizeAcceleroSignal+sizeGyroSignal;

static const int sizeFTSignal = 6;
static const int sizePoseSignal = 7;
static const int sizePoseSignalTangent = 6;


////////////DEFAULT VALUES FOR COVARIANCE //////
static const double positionInitVariance = !e-4;
static const double orientationInitVariance = 1e-4;
static const double linVelInitVariance = 1e-6;
static const double angVelInitVariance = 1e-6;
static const double forceInitVariance = 1e100;
static const double torqueInitVariance = 1e100;

static const double positionProcessVariance = 1e-8;
static const double orientationProcessVariance = 1e-8;
static const double linVelProcessVariance = 1e-8;
static const double angVelProcessVariance = 1e-8;
static const double forceProcessVariance = 1e-8;
static const double TorqueProcessVariance = 1e-8;
static const double TorqueProcessVarianceYaw = 1e-8;

static const double acceleroVariance = 1e-4;
static const double gyroVariance = 1e-8;
static const double forceSensorVariance = 1e-8;
static const double torqueSensorVariance = 1e-10;
static const double positionSensorVariance = 1e-4;
static const double orientationSensorVariance = 1e-3;

///default derivation steps
static const double defaultdx=1e-6;






namespace stateObservation
{

   /**
    * \class  EKFFlexibilityEstimatorBase
    * \brief  This class is the base class of the flexibility estimators that
    *         use an extended Kalman Filter. Several methods require to be overloaded
    *         to derive an implementation from this base class.
    *
    */

    class KineticsObserver: 
        protected DynamicalSystemFunctorBase
    {
    public:
        typedef kine::Kinematics Kinematics;

        /// The constructor.
        ///  \li maxContacts : maximum number of contacts,
        ///  \li dx gives the derivation step for a finite differences derivation method
        KineticsObserver(int maxContacts=4);


        ///virtual destructor
        virtual ~KineticsObserver();


        ///Gets the state size
        virtual unsigned getStateSize() const;

        ///Gets the measurement size
        virtual unsigned getMeasurementSize() const;

        /// ///////////////////////////////////////////////////////////
        /// Getting, setting the current time and running the estimation
        /// //////////////////////////////////////////////////////////

        /// gets the sampling time
        virtual double getSamplingTime() const;

        /// sets the sampling time
        virtual void setSamplingTime(double) ;

        /// this function triggers the estimation itself
        void update();

        /// ////////////////////////////////
        /// Getting and setting the state
        /// ////////////////////////////////

        /// Gets an estimation of the state in the form of a state vector $\hat{x_{k+1}}$
        virtual Vector getStateVector() const;

        /// Get the Kinematics of the observed frams
        virtual Kinematics getKinematics() const;

        /// Get the kinematics of a given frame
        virtual Kinematics getKinematics(const Kinematics &) const;

        ///get the contact force provided by the estimator 
        /// which is different from a contact sensor measurement
        virtual Vector3 getContactForce(int contactNbr);
        virtual Vector3 getContactMoment(int contactNbr);
        virtual Vector6 getContactForceandMoment(int contactNbr);


        virtual Vector getExternalForces() const;



        ///This function allows to estimate the acceleration even if
        /// withAccelerationEstimation is not set
        virtual void estimateAcceleration();

        ///Sets a value of the state x_k provided from another source
        /// can be used for initialization of the estimator or use 
        /// setStateVector
        virtual void setStateVector(Vector3,bool resetCovariance=true);


        ///Sets a value for the kinematics part of the state
        virtual void setKinematics(const Kinematics &, bool resetForces=true,
                                   bool resetCovariance=true);

        /// /////////////////////////////////////////////////////
        /// Setting and getting the state of the estimation
        /// ////////////////////////////////////////////////////
        void setWithExternalForces(bool b = true);

        void setWithAccelerationEstimation(bool b = true);

        /// ///////////////////////////////////////////////
        /// Getting and setting input data and measurements
        /// /////////////////////////////////////////////

        /// sets the measurement of the IMU (gyrometer, accelerometer, kinematics and number of the IMU)
        /// accelero and gyrometer are the measurement
        /// localKine gets the kinematics of the IMU expressed in the observed frame
        /// the best is to provide the position, the orientation, 
        /// the angular and linear velocities and the linear acceleration 
        /// Nevertheless if velocities or accelerations are not available they will be 
        /// automatically computed through finite differences
        /// the acceleroCov and the gyroCov are the covariance matrices of the sensors
        virtual int setIMU(const Vector3 & accelero, const  Vector3 & gyrometer, const Kinematics &localKine, int num=-1);
        virtual int setIMU(const Vector3 & accelero, const  Vector3 & gyrometer, const Matrix3& acceleroCov, 
                                                        const Matrix3 gyroCov, const Kinematics &localKine, int num=-1);
        

        
        virtual void setIMUDefaultCovarianceMatrix(const Matrix3& acceleroCov, const Matrix3 gyroCov);


        ////////// Contact stters, these are MANDATORY for every contact at every iteration ///
        /// if the contact is equipped with FT sensor call setContactFTSensor 
        /// otherwise calls etContactWithoutSensor

        /// sets the measurements of a force or torque sensor of the contact numbered contactNumber
        /// localKine sets the kinematics of the contact expressed in the observed frame
        /// the best is to provide the position , the orientation, 
        /// the angular and the linear velocities. 
        virtual void setContactFTSensor(const Vector3& force, const Vector3& torque, const Kinematics &localKine, int contactNumber);
        virtual void setContactFTSensor(const Vector3& force, const Vector3& torque, const Matrix6 ForcetorqueCovMatrix, 
                                                                                    const Kinematics &localKine, int contactNumber);

        virtual void setContactFTSensorDefaultCovarianceMatrix(const Matrix6 & forceTorqueSensorCovMat);
        
        /// sets the the kinematics of the contact expressed in the observed frame
        /// the best is to provide the position, the orientation, the angular and the linear velocities. 
        virtual void setContactWithNoSensor(const Kinematics &localKine, int contactNumber);
                                                               
        /// Set a measurement of the pose. The input is the Measured kinematics 
        /// namely position and orientation.
        virtual void setPoseSensor(const Kinematics &);
        virtual void setPoseSensor(const Kinematics &, const Matrix6 CovarianceMatrix);

        virtual void setPoseSensorDefaultCovarianceMatrix(const Matrix & , int contactNumber);


        ///if only force of torque is available, set the unavailable value to zero 
        virtual void setExternalForceTorqueMeasurement(const Vector3 force = Vector3::Zero(), const Vector3 torque = Vector3::Zero());


        ////////////// Contact management ///////////////////////////////
        /// Set a new contact
        /// -pose is the initial guess on the position of the contact. Only position 
        /// and orientation are enough
        /// -initialCovarianceMatrix is the covariance matrix expressing the 
        ///  uncertainty of the initial guess (if no initial guess is available 
        ///  give a rough position with a high initial covariance matrix)
        /// -processCovarianceMatrix is the covariance matrix expressing the 
        ///  rate at which the contact slides (set to zero for no sliding)
        /// -linear, angular stiffness and damping set the flexibility model
        ///  of the contact
        /// set contactNumber to -1 in order to set the number automatically
        /// returns the number of this contact
        virtual int addContact(const Kinematics & pose, 
                            const Matrix6 & initialCovarianceMatrix, const Matrix6 & processCovarianceMatrix, 
                            const Matrix3 & linearStiffness, const Matrix & linearDamping, 
                            const Matrix3 & angularStiffness, const Matrix & angularDamping, 
                            int contactNumner=-1);
        virtual int addContact(const Kinematics & pose, 
                            const Matrix6 & processCovarianceMatrix, 
                            const Matrix3 & linearStiffness, const Matrix & linearDamping, 
                            const Matrix3 & angularStiffness, const Matrix & angularDamping, 
                            int contactNumner=-1);
        virtual int removeContact(int contactnbr);

        virtual void clearContacts();

        virtual int getNumberOfContact() const;

        std::vector<int> getListOfContact() const;

        

        ///Sets the covariance matrix of the flexibility Guess
        virtual void setStateCovariance(const Matrix & P);

        ///Gets the covariance matrix of the flexibility
        virtual Matrix getStateCovariance() const;

        ///Sets the covariance matrices for the process noises
        /// \li Q process noise
        virtual void setProcessNoiseCovariance(const Matrix & Q);

        ///gets the covariance matrices
        virtual Vector getMeasurementVector();

        /// Gets a const reference on the extended Kalman filter
        virtual const stateObservation::ExtendedKalmanFilter & getEKF() const;

        /// Gets a reference on the extended Kalman filter
        /// modifying this object may lead to instabilities
        virtual stateObservation::ExtendedKalmanFilter & getEKF();

        ///Resets the covariance matrices to their original values
        virtual void resetCovarianceMatrices();

        ///to reset all the sensor inputs and provided contact positions
        void resetIteration();


protected:
        ///////////// DYNAMICAL SYSTEM IMPLEMENTATION
        virtual Vector stateDynamics(const Vector &x, const Vector &u, TimeIndex k);

        virtual Vector measureDynamics(const Vector &x, const Vector &u, TimeIndex k);

        ///Sets a noise which disturbs the state dynamics
        virtual void setProcessNoise(NoiseBase *);

        ///Removes the process noise
        virtual void resetProcessNoise();
        ///Gets the process noise
        virtual NoiseBase *getProcessNoise() const;

        ///Sets a noise which disturbs the measurements
        virtual void setMeasurementNoise(NoiseBase *);
        ///Removes the measurement noise
        virtual void resetMeasurementNoise();
        ///Gets a pointer on the measurement noise
        virtual NoiseBase *getMeasurementNoise() const;

        ///Set the period of the time discretization
        virtual void setSamplingPeriod(double dt);

        virtual Matrix getAMatrix(const Vector &xh);
        virtual Matrix getCMatrix(const Vector &xp);

        ///Gets the input size
        virtual unsigned getInputSize() const;

        static void stateSum(const  Vector& stateVector, const Vector& tangentVector, Vector& sum);

        static void stateDifference(const Vector& stateVector1, const Vector& stateVector2, Vector& difference);

        static void measureDifference(const Vector& measureVector1, const Vector& measureVector2, Vector& difference);
        /////////////////////////////////////// 

        virtual void setFiniteDifferenceStep(const Vector & dx);
        virtual void useFiniteDifferencesJacobians(bool b=true);



    protected:

        struct Sensor
        {
            Sensor(int signalSize):size(signalSize), time(0) {}
            virtual ~Sensor(){}
            int index;
            int size;
            int time;
            
            inline Vector extractFromVector(const Vector & v){return v.segment(size,index);}
        };

        struct IMU:
        public Sensor
        {
            virtual ~IMU(){}
            IMU():Sensor(sizeIMUSignal){}
            Kinematics kinematics;
            Vector6 acceleroGyro;
            Matrix3 covMatrixAccelero;
            Matrix3 covMatrixGyro;
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW        
        };

        typedef std::map<int, IMU, std::less<int>, 
             Eigen::aligned_allocator<std::pair<const int, IMU> > > MapIMU;
        typedef MapIMU::iterator MapIMUIterator;
        typedef MapIMU::const_iterator MapIMUConstIterator;
        typedef std::pair<int, IMU> PairIMU;

        
        struct Contact:
        public Sensor
        {
            Contact():Sensor(sizeFTSignal),withRealSensor(false){}
            virtual ~Contact(){}
            Kinematics kinematics;
            Vector6 forceTorque;
            Matrix6 covMatrix;
            bool withRealSensor;

            static int numberOfRealSensors;


            EIGEN_MAKE_ALIGNED_OPERATOR_NEW        
        };

        struct stateContainer
        {
            Kinematics kine;

            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        };

        typedef std::map <int, Contact, std::less<int>, 
                Eigen::aligned_allocator<std::pair<const int, Contact> > > MapContact;
        typedef MapContact::iterator MapContactIterator;
        typedef MapContact::const_iterator MapContactConstIterator;
        typedef std::pair<int, Contact> PairContact;

        struct AbsolutePoseSensor:
        public Sensor
        {
            AbsolutePoseSensor():Sensor(sizePoseSignal),isSet(false){}

            Kinematics pose;
            Matrix6 covMatrix;
            bool isSet;
        };

        AbsolutePoseSensor absPoseSensor_;
        MapIMU imuSensors_;
        MapContact contacts_;

        int stateSize_; 
        int stateTangentSize_;
        int measurementSize_; 

        double dt_;

        Vector stateVector_;
        Vector stateVectorDx_;
        Vector oldStateVector_; 

        Vector measurementVector_;
        Matrix measurementCovMatrix_;    

        Matrix3 acceleroDefaultCovMat_;
        Matrix3 gyroCovDefaultMat_;
        Matrix6 contactFTSensorDefaultCovMat_;
        Matrix6 poseSensorCovMat_;
        
        Matrix3 estPositionCovMat;
        Matrix3 estOrientationCovMat;
        Matrix3 estLinVelCovMat;
        Matrix3 estAngVelCovMat;
        

        bool finiteDifferencesJacobians_;

        stateObservation::ExtendedKalmanFilter ekf_;

        TimeIndex k_est;
        TimeIndex k_data;

        ///function to call before all the measurements
        ///detects if there is a new estimation beginning and then
        ///calls the reset of the iteration
        virtual void startNewIteration_();

    private:

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    };

}

#endif // FLEXIBILITYESTIMATION_EKFFLEXIBILITYESTIMATORBASE_H
