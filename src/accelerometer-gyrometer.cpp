#include <state-observation/sensors-simulation/accelerometer-gyrometer.hpp>

namespace stateObservation
{
    AccelerometerGyrometer::AccelerometerGyrometer():
    r_(Matrix3::Zero()),
    acc_(Vector3::Zero()),
    omega_(Vector3::Zero()),
    output_(Vector::Zero(measurementSize_,1)),
    currentStateSize_(stateSize_)
    {
#ifdef STATEOBSERVATION_VERBOUS_CONSTRUCTORS
      std::cout<<std::endl<<"AccelerometerGyrometer Constructor"<<std::endl;
#endif //STATEOBSERVATION_VERBOUS_CONSTRUCTOR

      matrixMode_=false;
    }


    size_t AccelerometerGyrometer::getStateSize_() const
    {
        return currentStateSize_;
    }

    size_t AccelerometerGyrometer::getMeasurementSize_() const
    {
        return measurementSize_;
    }

    Vector AccelerometerGyrometer::computeNoiselessMeasurement_()
    {
      if (!matrixMode_)
      {
        Quaternion q(state_.head<4>());

        r_=q.toRotationMatrix();
        acc_= state_.segment(4,3);
        omega_= state_.tail(3);
      }
      else
      {
        r_=Eigen::Map<Matrix3>(&state_[0]);
        acc_= state_.segment<3>(9);
        omega_ = state_.tail<3>();
      }

      output_.head<3>()=accelerationMeasure(acc_,r_);
      output_.tail<3>()=rotationVelocityMeasure(omega_, r_);

      return output_;
    }

    void AccelerometerGyrometer::setMatrixMode(bool matrixMode)
    {
      matrixMode_=matrixMode;
      if (!matrixMode_)
        currentStateSize_= stateSize_;
      else
        currentStateSize_= stateSizeMatrix_;
    }


}
