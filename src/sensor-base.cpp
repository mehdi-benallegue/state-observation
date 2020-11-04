#include <state-observation/sensors-simulation/sensor-base.hpp>

namespace stateObservation
{
    SensorBase::SensorBase()
        :noise_(0x0)
    {
    }

    bool SensorBase::checkStateVector(const Vector & state) const
    {
        return (state.rows()==getStateSize() && state.cols()==1);
    }

    void SensorBase::setNoise(NoiseBase * noise)
    {
        noise_=noise;
    }

    void SensorBase::resetNoise()
    {
        noise_=0x0;
    }

    Vector SensorBase::stateVectorZero()const
    {
        return Vector::Zero(getStateSize(),1);
    }

    NoiseBase * SensorBase::getNoise() const
    {
        return noise_;
    }

}
