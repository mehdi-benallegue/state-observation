/**
 * \file      kalman-filter-base.hpp
 * \author    Mehdi Benallegue
 * \date       2012
 * \brief      Defines the base class of a Kalman filter
 *
 *             It mostly implements the equations of Kalman filtering
 *             It is suitablle by derivation to be used incases of Linear,
 *             linearized and extended Kalman filtering. It may be
 *             derived to unscented Kalman filtering, but non-straighforwardly
 *             because the state vector is modified.
 *
 *             x_{k+1}=f(x_k,u_k)+v_k
 *
 *             y_k=h(x_k,u_k)+w_k
 *
 * \details
 *
 *
 */



#ifndef KALMANFILTERBASEHPP
#define KALMANFILTERBASEHPP

#include <state-observation/api.h>
#include <state-observation/observer/zero-delay-observer.hpp>
#include <state-observation/tools/state-vector-arithmetics.hpp>


namespace stateObservation
{


/**
     * \class  KalmanFilterBase
     * \brief
     *        It mostly implements the equations of Kalman filtering
     *        It is suitablle by derivation to be used incases of Linear,
     *        linearized and extended Kalman filtering. It may be
     *        derived to unscented Kalman filtering, but non-straighforwardly
     *        because the state vector is modified. This class requires
     *        to be derived to overload the update routine and the measurements
     *        simulation routine.
     *
     *             x_{k+1}=f(x_k,u_k)+v_k
     *
     *             y_k=h(x_k,u_k)+w_k
     *
     *
     *
     */
    class STATE_OBSERVATION_DLLAPI KalmanFilterBase: public ZeroDelayObserver,
                            protected StateVectorArithmetics
    {
    public:

        /// The type of the jacobian df/dx
        typedef Matrix Amatrix;

        /// The type of the jacobian dh/dx
        typedef Matrix Cmatrix;

        /// The type of the covariance matrix of the process noise v
        typedef Matrix Qmatrix;

        /// The type of the covariance matrix of the measurement noise w
        typedef Matrix Rmatrix;

        /// The type of the covariance matrix of the state estimation error.
        typedef Matrix Pmatrix;

        typedef Eigen::LLT<Pmatrix> LLTPMatrix;

        ///StateVector is the type of state tangent vector
        typedef Vector StateVectorTan;

        ///MeasureVector is the type of measurement tanegnt vector
        typedef Vector MeasureVectorTan;


        /// Default constructor
        KalmanFilterBase();

        /// The constructor
        ///  \li n : size of the state vector
        ///  \li m : size of the measurements vector
        ///  \li p : size of the input vector
        KalmanFilterBase(Index n,Index m,Index p=0);

        /// The constructor to use in case the dimension of the state space
        /// is smaller that its vector representation. For example
        /// The state could be made of rotations matrices (3x3 matrix: size = 9)
        /// Or quaternions (size =4) while they lie in a 3D space.
        /// In general the representation is in a Lie group and the representation
        /// of state derivatives are expressed in a Lie algebra.
        /// Use setSumFunction for Kalman update (mandatory)
        ///
        /// The update can then be done using exponential maps.
        ///  \li n : size of the state vector representation
        ///  \li nt : dimension of the tangent space to the state space
        ///  \li m : size of the measurements vector
        ///  \li mt : dimension of the tangent space to the measurement space
        ///  \li p : size of the input vector
        KalmanFilterBase(Index n, Index nt, Index m, Index mt, Index p);

        /// Set the value of the jacobian df/dx
        virtual void setA(const Amatrix& A);

        virtual Matrix getA() const;

        /// Clear the jacobian df/dx
        virtual void clearA();


        /// Set the value of the Jacobian dh/dx
        virtual void setC(const Cmatrix& C);

        virtual Matrix getC() const;

        /// Clear the jacobian dh/dx
        virtual void clearC();


        /// Set the measurement noise covariance matrix
        virtual void setR(const Rmatrix& R);
        inline void setMeasurementCovariance(const Rmatrix& R){setR(R);}

        /// Set the measurement noise covariance matrix
        virtual Matrix getR() const ;
        inline Matrix getMeasurementCovariance() const {return getR();}

        /// Clear the measurement noise covariance matrix
        virtual void clearR();


        /// Set the process noise covariance matrix
        virtual void setQ(const Qmatrix& Q);
        inline void setProcessCovariance(const Qmatrix& Q){setQ(Q);}

        /// Set the process noise covariance matrix
        virtual Matrix getQ() const;
        inline Matrix getProcessCovariance() const {return getQ();}


        /// Clear the process noise covariance matrix
        virtual void clearQ();


        /// Set the covariance matrix of the current time state estimation error
        virtual void setStateCovariance(const Pmatrix& P);

        /// Clear the covariace matrix of the current time state estimation
        /// error
        virtual void clearStateCovariance();

        /// Get the covariance matrix of the current time state estimation
        virtual Pmatrix getStateCovariance() const;

        /// Resets all the observer
        virtual void reset();


        /// Get a matrix having the size of the A matrix having "c" values
        Amatrix getAmatrixConstant(double c) const;

        /// Get a matrix having the size of the A matrix having random values
        Amatrix getAmatrixRandom() const;

        /// Get a matrix having the size of the A matrix having zero values
        Amatrix getAmatrixZero() const;

        /// Get an identity matrix having the size of the A matrix
        Amatrix getAmatrixIdentity() const;

        ///checks whether or not a matrix has the dimensions of the A matrix
        bool checkAmatrix(const Amatrix & ) const;


        /// Get a matrix having the size of the C matrix having "c" values
        Cmatrix getCmatrixConstant(double c) const;

        /// Get a matrix having the size of the C matrix having random values
        Cmatrix getCmatrixRandom() const;

        /// Get a matrix having the size of the C matrix having zero values
        Cmatrix getCmatrixZero() const;

        ///checks whether or not a matrix has the dimensions of the C matrix
        bool checkCmatrix(const Cmatrix &) const;


        /// Get a matrix having the size of the Q matrix having "c" values
        Qmatrix getQmatrixConstant(double c) const;

        /// Get a matrix having the size of the Q matrix having random values
        Qmatrix getQmatrixRandom() const;

        /// Get a matrix having the size of the Q matrix having zero values
        Qmatrix getQmatrixZero() const;

        /// Get an identity matrix having the size of the Q matrix
        Qmatrix getQmatrixIdentity() const;

        ///checks whether or not a matrix has the dimensions of the Q matrix
        bool checkQmatrix(const Qmatrix &) const;


        /// Get a matrix having the size of the R matrix having "c" values
        Rmatrix getRmatrixConstant(double c) const;

        /// Get a matrix having the size of the R matrix having random values
        Rmatrix getRmatrixRandom() const;

        /// Get a matrix having the size of the R matrix having zero values
        Rmatrix getRmatrixZero() const;

        /// Get an identity matrix having the size of the R matrix
        Rmatrix getRmatrixIdentity() const;

        ///checks whether or not a matrix has the dimensions of the R matrix
        bool checkRmatrix(const Rmatrix &) const;


        /// Get a matrix having the size of the P matrix having "c" values
        Pmatrix getPmatrixConstant(double c) const;

        /// Get a matrix having the size of the P matrix having random values
        Pmatrix getPmatrixRandom() const;

        /// Get a matrix having the size of the P matrix having zero values
        Pmatrix getPmatrixZero() const;

        /// Get an identity matrix having the size of the P matrix
        Pmatrix getPmatrixIdentity() const;

        /// Checks whether or not a matrix has the dimensions of the P matrix
        bool checkPmatrix(const Pmatrix & ) const;

        ///Gives a vector of state tangent vector size having duplicated "c" value
        virtual StateVectorTan stateTangentVectorConstant( double c ) const;

        ///Gives a vector of state  tangent vector size having random values
        virtual StateVectorTan stateTangentVectorRandom() const;

        ///Gives a vector of state  tangent vector size having zero values
        virtual StateVectorTan stateTangentVectorZero() const;

        ///Tells whether or not the vector has the dimensions of a state tangent vector
        virtual bool checkStateTangentVector(const StateVectorTan & v ) const;

        ///Gives a vector of measurement tangent vector size having duplicated "c" value
        virtual MeasureVectorTan measureTangentVectorConstant( double c ) const;

        ///Gives a vector of measurement tangent vector size having random values
        virtual MeasureVectorTan measureTangentVectorRandom() const;

        ///Gives a vector of measurement tangent vector size having zero values
        virtual MeasureVectorTan measureTangentVectorZero() const;

        ///Tells whether or not the vector has the dimensions of a measurement tangent vector
        virtual bool checkMeasureTangentVector(const MeasureVectorTan &) const;

        /// Changes the dimension of the state vector:
        ///resets the internal container for the state vector and
        ///the containers for the matrices A, C, Q, P
        virtual void setStateSize(Index n);

        /// Changes the dimension of the state vector:
        /// n is the dimension of the state representation
        /// and nt is the dimension of the tangent vector representation
        ///resets the internal container for the state vector and
        ///the containers for the matrices A, C, Q, P
        virtual void setStateSize(Index n, Index nt);

        /// Changes the dimension of the measurement vector:
        ///resets the internal container for the measurement vectors and
        ///the containers for the matrices C, R
        virtual void setMeasureSize(Index m);

        /// Changes the dimension of the measurement vector:
        /// m is the size of the measurementVector
        /// mt is the dimension of the measurement tangent space
        ///resets the internal container for the measurement vectors and
        ///the containers for the matrices C, R
        virtual void setMeasureSize(Index m, Index mt);

        /// Get simulation of the measurement y_k using the state estimation
        virtual MeasureVector getSimulatedMeasurement(TimeIndex k);

        ///Get the last vector of innovation of the Kalman filter
        virtual StateVector getInnovation();

        /// A function that gives the prediction (this is NOT the estimation of the state),
        /// for the estimation call getEstimateState method
        /// it is only an execution of the state synamics with the current state
        /// estimation and the current input value
        inline StateVector updateStatePrediction();

        ///update the predicted state, enables to precompute the predicted measurementŔ
        ///triggers also Vector updateStatePrediction()
        ///returns the measurement prediction
        ///definition in the bottom of this file
        inline MeasureVector updateStateAndMeasurementPrediction();

        ///get the last predicted state
        StateVector getLastPrediction() const;

        ///get the last predicted measurement
        MeasureVector getLastPredictedMeasurement() const;

        ///get the last Kalman gain matrix
        Matrix getLastGain() const;

        ///set update functions for sum and difference for the state vector
        /// (used for the case of multiplicative Kalman filter)
        void setStateArithmetics(StateVectorArithmetics * arith);

    protected:

        /// the size of tangent space of the state space
        Index nt_;

        /// the size of tangent space of the measurement space
        Index mt_;

        /// The type of Kalman gain matrix
        typedef Matrix Kmatrix;

        /// The Kalman filter loop
        virtual StateVector oneStepEstimation_();

        /// The abstract method to overload to implement f(x,u)
        virtual StateVector prediction_(TimeIndex k)=0;

        /// The abstract method to overload to implement h(x,u)
        virtual MeasureVector simulateSensor_(const StateVector& x, TimeIndex k)=0;

        /// Predicts the sensor measurement,
        virtual MeasureVector predictSensor_(TimeIndex k);

        /// Containers for the jacobian matrix of the process
        Matrix a_;

        /// Containers for the jacobian matrix of the measurement
        Matrix c_;

        /// Container for the process noise covariance matrice
        Matrix q_;

        /// Container for the measurement noise covariance matrice
        Matrix r_;

        /// Container for the covariance matrix of the estimation error
        Matrix pr_;

        /// container for the prediction
        IndexedVector xbar_;

        /// container for the prediction of the sensor
        IndexedVector ybar_;

        ///Vector containing the inovation of the Kalman filter
        Vector innovation_;

        struct optimizationContainer
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
            Matrix pbar;
            Vector xhat;
            Vector inoMeas;
            Matrix inoMeasCov;
            Matrix inoMeasCovInverse;
            LLTPMatrix inoMeasCovLLT;
            Matrix kGain;
            Matrix t;
        } oc_;

        StateVectorArithmetics * arithm_;

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    /*inline*/ Vector KalmanFilterBase::updateStatePrediction()
    {
        prediction_(this->x_.getTime()+1);
        return xbar_();
    }

    /*inline*/ Vector KalmanFilterBase::updateStateAndMeasurementPrediction()
    {
        updateStatePrediction();
        predictSensor_(this->x_.getTime()+1);
        return ybar_();
    }


}

#endif //KALMANFILTERBASEHPP
