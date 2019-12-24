#include <state-observation/dynamical-system/dynamical-system-functor-base.hpp>

namespace stateObservation
{

    DynamicalSystemFunctorBase::DynamicalSystemFunctorBase()
    {
#ifdef STATEOBSERVATION_VERBOUS_CONSTRUCTORS
       std::cout<<std::endl<<"DynamicalSystemFunctorBase Constructor"<<std::endl;
#endif //STATEOBSERVATION_VERBOUS_CONSTRUCTORS
        //ctor
    }

    DynamicalSystemFunctorBase::~DynamicalSystemFunctorBase()
    {
        //dtor
    }

    bool DynamicalSystemFunctorBase::checkStateVector(const Vector & v)
    {
        return (unsigned(v.rows())==getStateSize() && v.cols()==1);
    }

    bool DynamicalSystemFunctorBase::checkInputvector(const Vector & v)
    {
        return (unsigned(v.rows())==getInputSize() && unsigned(v.cols())==1);
    }

    void StateArithmetics::stateSum(const  Vector& stateVector, const Vector& tangentVector, Vector& sum)
    {
      detail::defaultSum(stateVector,tangentVector,sum);
    }

    void StateArithmetics::stateDifference(const Vector& stateVector1, const Vector& stateVector2, Vector& difference)
    {
      detail::defaultDifference(stateVector1,stateVector2,difference)  ;
    }

    void StateArithmetics::measurementDifference(const Vector& measureVector1, const Vector& measureVector2, Vector& difference)
    {
      detail::defaultDifference(measureVector1,measureVector2,difference);
    }

}
