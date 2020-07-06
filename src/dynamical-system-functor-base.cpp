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
        return (size_t(v.rows())==getStateSize() && v.cols()==1);
    }

    bool DynamicalSystemFunctorBase::checkInputvector(const Vector & v)
    {
        return (size_t(v.rows())==getInputSize() && size_t(v.cols())==1);
    }


}
