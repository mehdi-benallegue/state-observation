/*
 * dynamical-system.hpp
 *
 *  Created on: 19 mai 2014
 *      Author: alexis
 */

#ifndef DYNAMICAL_SYSTEM_HPP_
#define DYNAMICAL_SYSTEM_HPP_


#include <vector>

#include <state-observation/api.h>
#include <state-observation/dynamical-system/dynamical-system-functor-base.hpp>
#include <state-observation/noise/noise-base.hpp>
#include <state-observation/sensors-simulation/accelerometer-gyrometer.hpp>
#include <state-observation/tools/rigid-body-kinematics.hpp>
#include <state-observation/tools/hrp2.hpp>

#include <Eigen/Cholesky>


namespace stateObservation
{

  namespace flexibilityEstimation
  {

    /**
    * \class  DynamicalSystem
    * \brief  This class describes the dynamics of a robot's flexibility
    *         this dynamics is the simplest possible system, the flexibility
    *         is expressed as a rotation against the contact positions with no
    *         other hypothesis than that the contact points are at constant position
    *
    */
    class STATE_OBSERVATION_DLLAPI	IMUElasticLocalFrameDynamicalSystem :
      public stateObservation::DynamicalSystemFunctorBase
    {
    public:
      struct input
      {
        ///indexes of the different components of a vector of the input state
        static const unsigned posCom = 0;
        static const unsigned velCom = 3;
        static const unsigned accCom = 6;
        static const unsigned inertia = 9;
        static const unsigned angMoment = 15;
        static const unsigned dotInertia = 18;
        static const unsigned dotAngMoment = 24;
        static const unsigned posIMU = 27;
        static const unsigned oriIMU = 30;
        static const unsigned linVelIMU = 33;
        static const unsigned angVelIMU = 36;
        static const unsigned linAccIMU = 39;
        static const unsigned additionalForces = 42;
        static const unsigned contacts = 48;

        static const unsigned sizeBase = 48;
      };

      struct state
      {
        static const unsigned pos = 0;
        static const unsigned ori = 3;
        static const unsigned linVel = 6;
        static const unsigned angVel = 9;
        static const unsigned fc = 12;
        static const unsigned unmodeledForces = 24;
        static const unsigned comBias = 30;
        static const unsigned drift = 32;

        static const unsigned size =35;
      };


      struct contactModel
      {
        ///indexes of the different components of a vector of the input state
        static const unsigned elasticContact= 1;
        static const unsigned pendulum= 2;
        static const unsigned none= 0;

      };


      typedef Eigen::LLT<Matrix3> LLTMatrix3;

      ///constructor
      explicit IMUElasticLocalFrameDynamicalSystem(double dt);

      ///virtual destructor
      virtual ~IMUElasticLocalFrameDynamicalSystem();

      void test();

      // Get the contact wrench
      void computeContactWrench
              (const Matrix3& orientation, const Vector3& position,
               const IndexedVectorArray& contactPosV, const IndexedVectorArray& contactOriV,
               const Vector& fc, const Vector& tc, const Vector3 & fm, const Vector3& tm,
               const Vector3& addForce, const Vector3 & addMoment);

      stateObservation::Vector computeAccelerations(const Vector& x,const Vector& u);

      // computation of the acceleration linear
      virtual void computeAccelerations
      (const Vector3& positionCom, const Vector3& velocityCom,
       const Vector3& accelerationCom, const Vector3& AngMomentum,
       const Vector3& dotAngMomentum,
       const Matrix3& Inertia, const Matrix3& dotInertia,
       const IndexedVectorArray& contactPos,
       const IndexedVectorArray& contactOri,
       const Vector3& position, const Vector3& linVelocity,
       Vector3& linearAcceleration,  const Vector3& oriVector ,
       const Matrix3& orientation, const Vector3& angularVel,
       Vector3& angularAcceleration,
       const Vector& fc, const Vector& tc,
       const Vector3 & fm, const Vector3& tm,
       const Vector3 & addForces, const Vector3& addMoments);

      ///Description of the state dynamics
      virtual stateObservation::Vector stateDynamics
      (const stateObservation::Vector& x, const stateObservation::Vector& u,
       TimeIndex k);

      ///compute the jacobien of the state dynamics at the last computed value
      stateObservation::Matrix stateDynamicsJacobian();

      ///compute the jacobien of the state dynamics at a given state
      stateObservation::Matrix stateDynamicsJacobian(const stateObservation::Vector& x, const stateObservation::Vector& u,
       TimeIndex k);

      ///sets the finite differences derivation step vector
      void setFDstep(const stateObservation::Vector & dx);


      ///Description of the sensor's dynamics
      virtual stateObservation::Vector measureDynamics
      (const stateObservation::Vector& x, const stateObservation::Vector& u,
       TimeIndex k);

      ///compute the Jacobien of the measurements dynamics at the last computed value
      stateObservation::Matrix measureDynamicsJacobian();

      ///compute the Jacobien of the measurements dynamics at a given state value
      stateObservation::Matrix measureDynamicsJacobian(const stateObservation::Vector& x, const stateObservation::Vector& u,
       TimeIndex k);

      ///Sets a noise which disturbs the state dynamics
      virtual void setProcessNoise( stateObservation::NoiseBase * );

      ///Removes the process noise
      virtual void resetProcessNoise();

      ///Gets the process noise
      virtual stateObservation::NoiseBase * getProcessNoise() const;

      ///Sets a noise which disturbs the measurements
      virtual void setMeasurementNoise( stateObservation::NoiseBase * );

      ///Removes the measurement noise
      virtual void resetMeasurementNoise();

      ///Gets a pointer on the measurement noise
      virtual stateObservation::NoiseBase * getMeasurementNoise() const;

      ///Set the period of the time discretization
      virtual void setSamplingPeriod(double dt);

      ///Gets the state size
      virtual Index getStateSize() const;

      ///Gets the input size
      virtual Index getInputSize() const;

      ///Sets the input size
      virtual void setInputSize(Index i);

      ///Gets the contact number
      ///virtual

      ///Gets the contacts position

      ///Gets the measurement size
      virtual Index getMeasurementSize() const;

      ///Sets the number of contacts
      virtual void setContactsNumber(unsigned);

      virtual void setPe(stateObservation::Vector3 Pe)
      {
          pe=Pe;
      }

      ///Gets the nimber of contacts
      inline unsigned getContactsNumber(void) const
      {
         return  nbContacts_;
      }

      virtual void setContactModel(unsigned nb);

      virtual void setPrinted(bool b)
      {
          printed_ = b;
      }

      virtual void computeElastContactForcesAndMoments
      (const IndexedVectorArray& contactPosArray,
       const IndexedVectorArray& contactOriArray,
       const IndexedVectorArray& contactVelArray,
       const IndexedVectorArray& contactAngVelArray,
       const Vector3& position, const Vector3& linVelocity,
       const Vector3& oriVector, const Matrix3& orientation,
       const Vector3& angVel,
       Vector& fc, Vector& tc);

      virtual void computeElastPendulumForcesAndMoments
      (const IndexedVectorArray& PrArray,
       const Vector3& position, const Vector3& linVelocity,
       const Vector3& oriVector, const Matrix3& orientation,
       const Vector3& angVel,
       Vector& forces, Vector& moments);

      void computeForcesAndMoments
      (const IndexedVectorArray& position1,
       const IndexedVectorArray& position2,
       const IndexedVectorArray& velocity1,
       const IndexedVectorArray& velocity2,
       const Vector3& position, const Vector3& linVelocity,
       const Vector3& oriVector, const Matrix3& orientation,
       const Vector3& angVel,
       Vector& fc, Vector& tc);

      virtual void computeForcesAndMoments
      (const Vector& x,
       const Vector& u);


      virtual Vector getForcesAndMoments();

      virtual Vector getForcesAndMoments (const Vector& x,
       const Vector& u);

      virtual Vector getMomentaDotFromForces(const Vector& x, const Vector& u);
      virtual Vector getMomentaDotFromKinematics(const Vector& x, const Vector& u);

      virtual void iterateDynamicsEuler
      (const Vector3& positionCom, const Vector3& velocityCom,
       const Vector3& accelerationCom, const Vector3& AngMomentum,
       const Vector3& dotAngMomentum,
       const Matrix3& Inertia, const Matrix3& dotInertia,
       const IndexedVectorArray& contactPos,
       const IndexedVectorArray& contactOri,
       Vector3& position, Vector3& linVelocity, Vector& fc1,
       Vector3 &oriVector, Vector3& angularVel, Vector& fc2,
       const Vector3 & fm, const Vector3& tm,
       const Vector3 & addForces, const Vector3& addMoments,
       double dt
      );

      virtual void iterateDynamicsRK4
      (const Vector3& positionCom, const Vector3& velocityCom,
       const Vector3& accelerationCom, const Vector3& AngMomentum,
       const Vector3& dotAngMomentum,
       const Matrix3& Inertia, const Matrix3& dotInertia,
       const IndexedVectorArray& contactPos,
       const IndexedVectorArray& contactOri,
       Vector3& position, Vector3& linVelocity, Vector& fc1,
       Vector3 &oriVector, Vector3& angularVel, Vector& fc2,
       const Vector3 & fm, const Vector3& tm,
       const Vector3 & addForces, const Vector3& addMoments,
       double dt
      );

      virtual void setWithForceMeasurements(bool b);
      virtual bool getWithForceMeasurements() const;
      virtual void setWithComBias(bool b);
      virtual bool getWithComBias() const;
      virtual void setWithAbsolutePosition(bool b);
      virtual bool getWithAbsolutePosition() const;
      void setWithUnmodeledForces(bool b);


      virtual void setKfe(const Matrix3 & m);
      virtual void setKfv(const Matrix3 & m);
      virtual void setKte(const Matrix3 & m);
      virtual void setKtv(const Matrix3 & m);

      virtual void setKfeRopes(const Matrix3 & m);
      virtual void setKfvRopes(const Matrix3 & m);
      virtual void setKteRopes(const Matrix3 & m);
      virtual void setKtvRopes(const Matrix3 & m);

      virtual Matrix getKfe() const;
      virtual Matrix getKfv() const;
      virtual Matrix getKte() const;
      virtual Matrix getKtv() const;

      virtual void setRobotMass(double d);

      virtual double getRobotMass() const;

    protected:

      bool printed_;

      stateObservation::AccelerometerGyrometer sensor_;

      stateObservation::NoiseBase * processNoise_;

      void updateMeasurementSize_();

      double dt_;

      double robotMass_;
      double robotMassInv_;

      Matrix3& computeRotation_(const Vector3 & x, int i);

      static const Index stateSize_=state::size;
      Index inputSize_;
      static const Index measurementSizeBase_=6;
      unsigned nbContacts_;
      unsigned contactModel_;

      Vector fc_;
      Vector tc_;

      Vector dx_;

      Vector xk1_;
      Vector xk_;
      Vector uk_;

      Vector xk_fory_;
      Vector yk_;
      Vector uk_fory_;


      Index measurementSize_;

      std::vector <Vector3,Eigen::aligned_allocator<Vector3> > contactPositions_;

      Matrix3 Kfe_, Kte_, Kfv_, Ktv_;
      Matrix3 KfeRopes_, KteRopes_, KfvRopes_, KtvRopes_;

      TimeIndex kcurrent_;

      bool withForceMeasurements_;
      bool withComBias_;
      bool withAbsolutePos_;
      bool withUnmodeledForces_;

      stateObservation::Vector3 pe;

      double marginalStabilityFactor_;
      //a scaling factor a=1-epsilon to avoid the natural marginal stability of
      //the dynamics x_{k+1}=x_k we replace it with x_{k+1}=a*x_k

      unsigned index_;

      struct Optimization
      {

        Vector6 momentaDot;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        Vector3 positionFlex;
        Vector3 velocityFlex;
        Vector3 accelerationFlex;
        Vector3 orientationFlexV;
        Vector3 angularVelocityFlex;
        Vector3 angularAccelerationFlex;
        Vector3 positionComBias;

        Matrix3 rFlex;
        Matrix3 rFlexT;

        Vector3 drift;
        Vector3 pdrift;
        Matrix3 rdrift;

        double cy,sy;

        AngleAxis orientationAA;

        Vector xk1;
        Vector xk;
        Vector xk1dx;

        Vector xdx;

        Vector xk_fory;
        Vector yk;
        Vector ykdy;

        TimeIndex k_fory;

        Matrix3 orinertia;

        LLTMatrix3 invinertia;

        Matrix3 rtotal;
        Vector3 ptotal;
        AngleAxis aatotal;
        Vector3 oritotal;

        Matrix Jx;
        Matrix Jy;

        Matrix3 rimu;
        Vector3 imuAcc;
        Vector3 imuOmega;
        Vector sensorState;

        Vector3 positionCom;
        Vector3 velocityCom;
        Vector3 accelerationCom;
        Vector3 AngMomentum;
        Vector3 dotAngMomentum;

        Vector3 positionControl;
        Vector3 velocityControl;
        Vector3 accelerationControl;
        Vector3 orientationControlV;
        Vector3 angularVelocityControl;

        Matrix3 rControl;

        IndexedVectorArray contactPosV;
        IndexedVectorArray contactOriV;
        IndexedVectorArray contactVelArray;
        IndexedVectorArray contactAngVelArray;

        Matrix3 inertia;
        Matrix3 dotInertia;

        IndexedVectorArray efforts;

        Vector3 f, fi;
        Vector3 t;

        //unmodelled and unmeasured forces
        Vector3 fm;
        Vector3 tm;

        Vector3 linearAcceleration;
        Vector3 angularAcceleration;

        Vector3 vf;
        Vector3 vt;

        Vector3 crosstempV;
        Matrix3 crosstempM;

        //elastic contact forces and moments
        Matrix3 Rci; //rotation of contact i
        Matrix3 Rcit;//transpose of previous
        Vector3 contactPos; //
        Vector3 contactVel;
        Vector3 RciContactPos;
        Vector3 globalContactPos;
        Matrix3 Rt;

        Vector3 forcei;
        Vector3 momenti;

        //additional forces and moments
        Vector3 addForce;
        Vector3 addMoment;

        Matrix3 skewV;
        Matrix3 skewV2;
        Matrix3 skewVR;
        Matrix3 skewV2R;
        Matrix3 RIRT;
        Vector3 wx2Rc;
        Vector3 _2wxRv;
        Vector3 Ra;
        Vector3 Rc;
        Vector3 Rcp;

        //optimization of orientation transformation between vector3 to rotation matrix

        Matrix3 curRotation0;
        Vector3 orientationVector0;
        Matrix3 curRotation1;
        Vector3 orientationVector1;
        Matrix3 curRotation2;
        Vector3 orientationVector2;
        Matrix3 curRotation3;
        Vector3 orientationVector3;

        Optimization()
          :
          curRotation0(Matrix3::Identity()),
          orientationVector0(Vector3::Zero()),
          curRotation1(Matrix3::Identity()),
          orientationVector1(Vector3::Zero()),
          curRotation2(Matrix3::Identity()),
          orientationVector2(Vector3::Zero()),
          curRotation3(Matrix3::Identity()),
          orientationVector3(Vector3::Zero())
        {}

        inline Vector3& orientationVector(int i)
        {
          if (i==0)
            return orientationVector0;
          if (i==1)
            return orientationVector1;
          if (i==2)
            return orientationVector2;

          return orientationVector3;
        }

        inline Matrix3& curRotation(int i)
        {
          if (i==0)
            return curRotation0;
          if (i==1)
            return curRotation1;
          if (i==2)
            return curRotation2;

          return curRotation3;
        }

      } op_;

    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    };
  }
}



#endif /* DYNAMICAL_SYSTEM_HPP_ */
