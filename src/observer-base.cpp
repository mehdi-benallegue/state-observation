#include  <state-observation/observer/observer-base.hpp>

namespace stateObservation
{



    void ObserverBase::reset()
    {
        clearStates();
        clearMeasurements();
        clearInputs();
    }

    ObserverBase::ObserverBase()
    {
        n_=m_=p_=0;
    }

    ObserverBase::ObserverBase(Index n, Index m, Index p):
            n_(n),m_(m),p_(p)
    {

    }



    ObserverBase::StateVector ObserverBase::stateVectorConstant( double c ) const
    {
        return StateVector::Constant(n_,1,c);
    }

    ObserverBase::StateVector ObserverBase::stateVectorRandom() const
    {
        return StateVector::Random(n_,1);
    }

    ObserverBase::StateVector ObserverBase::stateVectorZero() const
    {
        return StateVector::Zero(n_,1);
    }

    bool ObserverBase::checkStateVector(const StateVector & v) const
    {
        return (v.rows()==n_ && v.cols()==1);
    }


    ObserverBase::MeasureVector ObserverBase::measureVectorConstant( double c ) const
    {
        return MeasureVector::Constant(m_,1,c);
    }

    ObserverBase::MeasureVector ObserverBase::measureVectorRandom() const
    {
        return MeasureVector::Random(m_,1);
    }

    ObserverBase::MeasureVector ObserverBase::measureVectorZero() const
    {
        return MeasureVector::Zero(m_,1);
    }

    bool ObserverBase::checkMeasureVector(const MeasureVector & v) const
    {
        return (v.rows()==m_ && v.cols()==1);
    }


    ObserverBase::InputVector ObserverBase::inputVectorConstant( double c ) const
    {
        return InputVector::Constant(p_,1,c);
    }

    ObserverBase::InputVector ObserverBase::inputVectorRandom() const
    {
        return InputVector::Random(p_,1);
    }

    ObserverBase::InputVector ObserverBase::inputVectorZero() const
    {
        return InputVector::Zero(p_,1);
    }

    bool ObserverBase::checkInputVector(const InputVector & v) const
    {
        return (v.rows()==p_ && v.cols()==1);
    }


    void ObserverBase::setStateSize(Index n)
    {
        n_=n;
    }

    Index ObserverBase::getStateSize() const
    {
        return n_;
    }

    void ObserverBase::setMeasureSize(Index m)
    {
        m_=m;
    }

    Index ObserverBase::getMeasureSize() const
    {
        return m_;
    }


    void ObserverBase::setInputSize(Index p)
    {
        p_=p;
    }

    Index ObserverBase::getInputSize() const
    {
        return p_;
    }

}
