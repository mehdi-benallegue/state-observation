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
        unsigned getStateSize() const;

        ///Gets the measurement size
        unsigned getMeasurementSize() const;

        /// ///////////////////////////////////////////////////////////
        /// Getting, setting the current time and running the estimation
        /// //////////////////////////////////////////////////////////

        /// gets the sampling time
        double getSamplingTime() const;

        /// sets the sampling time
        void setSamplingTime(double) ;

        /// sets the mass of the robot
        void setMass(double);

        /// this function triggers the estimation itself
        Vector update();

        /// ////////////////////////////////
        /// Getting and setting the state
        /// ////////////////////////////////

        /// Gets an estimation of the state in the form of a state vector $\hat{x_{k+1}}$
        Vector getStateVector() const;

        /// Get the Kinematics of the observed frams
        Kinematics getKinematics() const;

        /// Get the kinematics of a given frame
        Kinematics getKinematics(const Kinematics & localKinematics) const;
        Kinematics getKinematics(const Kinematics & localKinematics);

        ///get the contact force provided by the estimator 
        /// which is different from a contact sensor measurement
        Vector6 getContactWrench(int contactNbr) const;
        Kinematics getContactPosition(int contactNbr) const;

        /// gets the external unmodeled forces
        Vector6 getUnmodeledWrench() const;

        ///This function allows to estimate the acceleration 
        /// it returns a Kinemactis
        Kinematics estimateAccelerations();

        ///Sets a value for the kinematics part of the state
        /// if resetForces is set to true the forces are set to zero
        ///
        void setStateKinematics(const Kinematics &, bool resetContactWrenches=true,
                                   bool resetCovariance=true);

        void setGyroBias(const Vector3 &,  bool resetCovariance=true);

        void setStateUnmodeledWrench(const Vector6 &, bool resetCovariance=true);
        
        ///Sets a value of the state x_k provided from another source
        /// can be used for initialization of the estimator 
        void setStateVector(const Vector &,bool resetCovariance=true);

        /// Add known external forces and moments which are not due to contact
        /// they must be expressed in the same frame as the kinematic root
        void setAdditionalWrench(const Vector6& );

        /// Getters for the indexes of the state Vector
        inline unsigned kineIndex() const;
        inline unsigned posIndex() const;
        inline unsigned oriIndex() const;
        inline unsigned linVelIndex() const;
        inline unsigned angVelIndex() const;
        inline unsigned gyroBiasIndex() const;
        inline unsigned unmodeledWrenchIndex() const;
        inline unsigned unmodeledForceIndex() const;
        inline unsigned unmodeledTorqueIndex() const;
        inline unsigned contactKineIndex(int contactNbr) const;
        inline unsigned contactPosIndex(int contactNbr) const;
        inline unsigned contactOriIndex(int contactNbr) const;
        inline unsigned contactForceIndex(int contactNbr) const;
        inline unsigned contactTorqueIndex(int contactNbr) const;
        inline unsigned contactWrenchIndex(int contactNbr) const;        

        /// /////////////////////////////////////////////////////
        /// Setting and getting the state of the estimation
        /// ////////////////////////////////////////////////////
        void setWithUnmodeledWrench(bool b = true);

        void setWithAccelerationEstimation(bool b = true);

        void setWithGyroBias(bool b = true);

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
        int setIMU(const Vector3 & accelero, const  Vector3 & gyrometer, const Kinematics &localKine, int num=-1);
        int setIMU(const Vector3 & accelero, const  Vector3 & gyrometer, const Matrix3& acceleroCov, 
                                                        const Matrix3 gyroCov, const Kinematics &localKine, int num=-1);
                
        void setIMUDefaultCovarianceMatrix(const Matrix3& acceleroCov, const Matrix3 gyroCov);

        ////////// Contact stters, these are MANDATORY for every contact at every iteration ///
        /// if the contact is equipped with wrench sensor call setContactWrenchSensor 
        /// otherwise calls etContactWithoutSensor

        /// sets the measurements of a force or torque sensor of the contact numbered contactNumber
        /// localKine sets the kinematics of the contact expressed in the observed frame
        /// the best is to provide the position , the orientation, 
        /// the angular and the linear velocities. 
        void setContactWrenchSensor(const Vector3& force, const Vector3& torque, const Kinematics &localKine, int contactNumber);
        void setContactWrenchSensor(const Vector3& force, const Vector3& torque, const Matrix6 ForcetorqueCovMatrix, 
                                                                                    const Kinematics &localKine, int contactNumber);

        void setContactWrenchSensorDefaultCovarianceMatrix(const Matrix6 & forceTorqueSensorCovMat);
        
        /// sets the the kinematics of the contact expressed in the observed frame
        /// the best is to provide the position, the orientation, the angular and the linear velocities. 
        void setContactWithNoSensor(const Kinematics &localKine, int contactNumber);
                                                               
        /// Set a measurement of the pose. The input is the Measured kinematics 
        /// namely position and orientation.
        void setPoseSensor(const Kinematics &);
        void setPoseSensor(const Kinematics &, const Matrix6 CovarianceMatrix);
 
        void setPoseSensorDefaultCovarianceMatrix(const Matrix & , int contactNumber);

        ///if only force of torque is available, set the unavailable value to zero 
        virtual void setUnmodeledForceTorqueMeasurement(const Vector3 force = Vector3::Zero(), const Vector3 torque = Vector3::Zero());

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
        int addContact(const Kinematics & pose, 
                            const Matrix6 & initialCovarianceMatrix, const Matrix6 & processCovarianceMatrix, 
                            const Matrix3 & linearStiffness,  const Matrix3 & linearDamping, 
                            const Matrix3 & angularStiffness, const Matrix3 & angularDamping, 
                            int contactNumner=-1);
        /// version with default stiffness and damping
        /// use when the contact parameters are known
        int addContact(const Kinematics & pose, 
                            const Matrix6 & initialCovarianceMatrix, const Matrix6 & processCovarianceMatrix, 
                            int contactNumner=-1);
        /// version when the contact position is perfectly known
        int addContact(const Kinematics & pose, 
                            const Matrix3 & linearStiffness,  const Matrix3 & linearDamping, 
                            const Matrix3 & angularStiffness, const Matrix3 & angularDamping, 
                            int contactNumner=-1);
        /// version when the position is perfectly known but not the stiffness and damping
        int addContact(const Kinematics & pose, int contactNumner=-1);

        int removeContact(int contactnbr);

        void clearContacts();

        int getNumberOfContact() const;

        std::vector<int> getListOfContacts() const;

        ///Sets the covariance matrix of the flexibility Guess
        void setStateCovariance(const Matrix & P);

        void setKinematicsCovariance(const Matrix & );
        void setGyroBiasCovariance(const Matrix3 & );
        void setUnmodeledWrenchCovariance(const Matrix3 & );

        ///Gets the covariance matrix of the flexibility
        Matrix getStateCovariance() const;

        ///Sets the covariance matrices for the process noises
        /// \li Q process noise
        void setProcessNoiseCovariance(const Matrix & Q);

        ///gets the covariance matrices
        Vector getMeasurementVector();

        /// Gets a const reference on the extended Kalman filter
        const stateObservation::ExtendedKalmanFilter & getEKF() const;

        /// Gets a reference on the extended Kalman filter
        /// modifying this object may lead to instabilities
        stateObservation::ExtendedKalmanFilter & getEKF();

        ///Resets the covariance matrices to their original values
        void resetStateCovarianceMat();
        void resetStateKinematicsCovMat();
        void resetStateGyroBiasCovMat();
        void resetStateUnmodeledWrenchCovMat();       

        void resetProcessCovarianceMat();
        void resetProcessKinematicsCovMat();
        void resetProcessGyroBiasCovMat();
        void resetProcessUnmodeledWrenchCovMat();

        void resetSensorsCovMat();

        void resetAllCovairanceMatrices();

        ///removes all the contacts !!!        
        void resetContacts();

        ///to reset all the sensor inputs and provided contact positions but 
        void resetInputs();

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

        ////////////////////////////
        virtual void setFiniteDifferenceStep(const Vector & dx);
        virtual void useFiniteDifferencesJacobians(bool b=true);

        Vector stateNaNCorrection_();

        Vector6 computeAccelerations_();

        ///updates stateKine_ from the stateVector            
        void updateKine_();

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
            Contact():Sensor(sizeWrench),withRealSensor(false),stateIndex(-1){}
            virtual ~Contact(){}
            Kinematics kinematics;
            Vector6 forceTorque;
            Matrix6 covMatrix;
            bool withRealSensor;
            int stateIndex;

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
            AbsolutePoseSensor():Sensor(sizePose),isSet(false){}

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

        Kinematics stateKinematics_;

        Vector stateVector_;
        Vector stateVectorDx_;
        Vector oldStateVector_; 

        Vector6 additionalWrench_;

        Vector measurementVector_;
        Matrix measurementCovMatrix_;    

        Matrix3 estPositionCovMat;
        Matrix3 estOrientationCovMat;
        Matrix3 estLinVelCovMat;
        Matrix3 estAngVelCovMat;
        
        stateObservation::ExtendedKalmanFilter ekf_;
        bool finiteDifferencesJacobians_;
        bool withGyroBias_;
        bool withUnmodeledWrench_;

        TimeIndex k_est;
        TimeIndex k_data;

        double mass_;

        double dt_;

        /// resets one block on the diagonal of the state Covariance Matrix
        /// i.e. sets value of a square block on the diagonal of the covMat
        /// and sets to zero all the values related to their lines and columns
        template <int blockSize>
        static void setBlockStateCovariance(Matrix & covMat, const Matrix & covBlock, 
        int blockIndex, int matrixSize);

        ///function to call before all the measurements
        ///detects if there is a new estimation beginning and then
        ///calls the reset of the iteration
        void startNewIteration_();

    private:

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        ///////////SIZE OF VECTORS
        static const unsigned sizeAcceleroSignal = 3;
        static const unsigned sizeGyroSignal = 3;
        static const unsigned sizeIMUSignal = sizeAcceleroSignal+sizeGyroSignal;

        static const unsigned sizePos = 3;
        static const unsigned sizeOri = 4;
        static const unsigned sizeOriTangent = 3;
        static const unsigned sizeLinVel = sizePos;
        static const unsigned sizeAngVel = sizeOriTangent;
        static const unsigned sizeGyroBias = sizeGyroSignal;

        static const unsigned sizeForce = 3;
        static const unsigned sizeTorque = 3;

        static const unsigned sizeWrench = sizeForce + sizeTorque;


        static const unsigned sizeStateKine =sizePos+
                                             sizeOri+
                                             sizeLinVel+
                                             sizeAngVel;
        static const unsigned sizeStateBase = sizeStateKine+
                                              sizeGyroBias+
                                              sizeForce+
                                              sizeTorque;
        static const unsigned sizeStatePerContact = sizePos+
                                                    sizeOri+
                                                    sizeForce+
                                                    sizeTorque;
        static const unsigned sizeStateTangentBase = sizePos+
                                                    sizeOriTangent+
                                                    sizeLinVel+
                                                    sizeAngVel+
                                                    sizeGyroBias+
                                                    sizeForce+
                                                    sizeTorque;
        static const unsigned sizeStateTangentPerContact = sizePos+
                                                        sizeOriTangent+
                                                        sizeForce+
                                                        sizeTorque;

        static const unsigned sizePose = sizePos+
                                        sizeOri;
    
        static const unsigned sizeContactKine = sizePose;

        static const unsigned sizePoseTangent = sizePos+
                                                sizeOriTangent;


        static const Kinematics::Flags::Byte flagsStateKine =  Kinematics::Flags::position |
                                                               Kinematics::Flags::orientation |
                                                               Kinematics::Flags::linVel |
                                                               Kinematics::Flags::angVel;

        static const Kinematics::Flags::Byte flagsContactKine =  Kinematics::Flags::position |
                                                                 Kinematics::Flags::orientation;

        static const Kinematics::Flags::Byte flagsKineSensor = Kinematics::Flags::position |
                                                                Kinematics::Flags::orientation;

    ////////////DEFAULT VALUES //////
        static const double defaultMass;

        static const double statePoseInitVarianceDefault ;
        static const double stateOriInitVarianceDefault ;
        static const double stateLinVelInitVarianceDefault ;
        static const double stateAngVelInitVarianceDefault ;
        static const double gyroBiasInitVarianceDefault ;
        static const double unmodeledWrenchInitVarianceDefault ;
        static const double contactForceInitVarianceDefault ;
        static const double contactTorqueInitVarianceDefault ;

        static const double statePoseProcessVarianceDefault ;
        static const double stateOriProcessVarianceDefault ;
        static const double stateLinVelProcessVarianceDefault ;
        static const double stateAngVelProcessVarianceDefault ;
        static const double gyroBiasProcessVarianceDefault ;
        static const double unmodeledWrenchProcessVarianceDefault ;
        static const double contactForceProcessVarianceDefault ;
        static const double contactTorqueProcessVarianceDefault ;
        
        static const double acceleroVarianceDefault ;
        static const double gyroVarianceDefault ;
        static const double forceSensorVarianceDefault ;
        static const double torqueSensorVarianceDefault ;
        static const double positionSensorVarianceDefault ;
        static const double orientationSensorVarianceDefault;

        static const double linearStiffnessDefault;
        static const double angularStiffnessDefault;
        static const double linearDampingDefault;
        static const double angularDampingDefault;
        

    protected:
        /// Default Stiffness and damping
        Matrix3 linearStiffnessMatDefault_;
        Matrix3 angularStiffnessMatDefault_;
        Matrix3 linearDampingMatDefault_;
        Matrix3 angularDampingMatDefault_;

        ////////////Sensor Covariance mnatrices
        Matrix3 acceleroCovMatDefault_;
        Matrix3 gyroCovMatDefault_;
        Matrix6 contactWrenchSensorCovMatDefault_;
        Matrix6 poseSensorCovMatDefault_;

        Matrix3 statePosInitCovMat_;
        Matrix3 stateOriInitCovMat_;
        Matrix3 stateLinVelInitCovMat_;
        Matrix3 stateAngVelInitCovMat_;
        Matrix3 gyroBiasInitCovMat_;
        Matrix6 unmodeledWrenchInitCovMat_;
        Matrix6 contactWrenchInitCovMat_;

        Matrix3 statePosProcessCovMat_;
        Matrix3 stateOriProcessCovMat_;
        Matrix3 stateLinVelProcessCovMat_;
        Matrix3 stateAngVelProcessCovMat_;
        Matrix3 gyroBiasProcessCovMat_;
        Matrix6 unmodeledWrenchProcessCovMat_;
        Matrix6 contactWrenchProcessCovMat_;

        Matrix12 stateKineMatricsInitCovMat_;
        Matrix12 stateKineMatricsProcessCovMat_;

        ///default derivation steps
        static const double defaultdx;

    };

#include <state-observation/dynamics-estimators/kinetics-observer.hxx>

}

#endif // FLEXIBILITYESTIMATION_EKFFLEXIBILITYESTIMATORBASE_H
