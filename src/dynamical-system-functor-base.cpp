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


}
