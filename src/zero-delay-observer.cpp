#include <state-observation/observer/zero-delay-observer.hpp>

namespace stateObservation
{

    void ZeroDelayObserver::setState
                            (const ObserverBase::StateVector& x_k,TimeIndex k)
    {
        BOOST_ASSERT(checkStateVector(x_k)
                            && "The size of the state vector is incorrect");

        x_.set(x_k,k);
        while (y_.size()>0 && y_.getFirstIndex()<=k)
        {
            y_.popFront();
        }

        if (p_>0)
            while (u_.size()>0 && u_.getFirstIndex()<k)
            {
                u_.popFront();
            }
    }

    void ZeroDelayObserver::clearStates()
    {
        x_.reset();
    }

    void ZeroDelayObserver::setMeasurement
                    (const ObserverBase::MeasureVector& y_k,TimeIndex k)
    {

        BOOST_ASSERT(checkMeasureVector(y_k)
                && "The size of the measure vector is incorrect");
        if (y_.size()>0)
            BOOST_ASSERT ((y_.getNextIndex()==k || y_.checkIndex(k))
                && "ERROR: The time is set incorrectly for \
                                the measurements (order or gap)");
        else
            BOOST_ASSERT ( (!x_.isSet() || x_.getTime()==k-1)
                && "ERROR: The time is set incorrectly for the measurements \
                                (must be [current_time+1])");

        y_.setValue(y_k,k);
    }

    void ZeroDelayObserver::clearMeasurements()
    {
        y_.reset();
    }

    void ZeroDelayObserver::setInput
                    (const ObserverBase::InputVector& u_k,TimeIndex k)
    {
        if (p_>0)
        {
            BOOST_ASSERT(checkInputVector(u_k)
                        && "The size of the input vector is incorrect");

            if (u_.size()>0)
                BOOST_ASSERT ((u_.getNextIndex()==k || u_.checkIndex(k))
                        && "ERROR: The time is set incorrectly \
                                for the inputs (order or gap)");
            else
            {
                BOOST_ASSERT ( (!x_.isSet() || x_.getTime()==k || x_.getTime()==k-1)
                    && "ERROR: The time is set incorrectly for the \
                          inputs (must be [current_time] or [current_time+1])");
            }

            u_.setValue(u_k,k);
        }
    }

    void ZeroDelayObserver::clearInputs()
    {
        if (p_>0)
            u_.reset();
    }

    ObserverBase::StateVector
    ZeroDelayObserver::getEstimatedState(TimeIndex k)
    {
        TimeIndex k0=x_.getTime();

        BOOST_ASSERT(k0<=k
                && "ERROR: The observer cannot estimate previous states");

        for (TimeIndex i=k0;i<k;++i)
        {
            oneStepEstimation_();
            if (y_.getFirstIndex()<k)
                y_.popFront();

            if (p_>0)
                if (u_.getFirstIndex()<k)
                    u_.popFront();
        }

        return x_();
    }


    TimeIndex ZeroDelayObserver::getCurrentTime()const
    {
        return x_.getTime();
    }

    Vector ZeroDelayObserver::getInput(TimeIndex k) const
    {
        return u_[k];
    }

    TimeSize ZeroDelayObserver::getInputsNumber()const
    {
        return u_.size();
    }

    TimeIndex ZeroDelayObserver::getInputTime()const
    {
        if (u_.size()>0)
        {
            return u_.getLastIndex();
        }
        else
        {
            return 0;
        }
    }

    Vector ZeroDelayObserver::getMeasurement(TimeIndex k) const
    {
        return y_[k];
    }

    TimeIndex ZeroDelayObserver::getMeasurementTime()const
    {
        BOOST_ASSERT(y_.size()>0
                && "ERROR: There is no measurements registered (past measurements are erased)");
        return y_.getLastIndex();
    }

    TimeSize ZeroDelayObserver::getMeasurementsNumber()const
    {
        return TimeSize(y_.size());
    }


    void ZeroDelayObserver::setStateSize(Index n)
    {
        if (n!=n_)
        {
            ObserverBase::setStateSize(n);
            clearStates();
        }
    }

    void ZeroDelayObserver::setMeasureSize(Index m)
    {
        if (m!=m_)
        {
            ObserverBase::setMeasureSize(m);
            clearMeasurements();
        }
    }

    void ZeroDelayObserver::setInputSize(Index p)
    {
        if (p!=p_)
        {
            ObserverBase::setInputSize(p);
            clearInputs();
        }
    }
}
