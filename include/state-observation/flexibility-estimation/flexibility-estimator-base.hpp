/**
 * \file     flexibility-estimator-base.hpp
 * \author   Mehdi Benallegue
 * \date     2013
 * \brief    Definitions of base class for flexibility estimator.
 *
 * \details
 *
 *
 */

#ifndef FLEXIBILITY_ESTIMATION_FLEXIBILITY_ESTIMATOR_BASE_HPP
#define FLEXIBILITY_ESTIMATION_FLEXIBILITY_ESTIMATOR_BASE_HPP

#include <state-observation/api.h>
#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{
namespace flexibilityEstimation
{
    /**
    * \class  FlexibilityEstimatorBase
    * \brief  This class is the base class of the flexibility estimators.
    *
    */
    class STATE_OBSERVATION_DLLAPI FlexibilityEstimatorBase
    {
    public:
        /// virtual destructor
        virtual ~FlexibilityEstimatorBase(){}

        ///The constructor
        explicit FlexibilityEstimatorBase();

        ///Sets a value of the flexibility x_k provided from another source
        /// can be used for initialization of the estimator
        virtual void setFlexibilityGuess(const Matrix &)=0;

        /// Sets the value of the next sensor measurement y_{k+1}
        virtual void setMeasurement(const Vector &)=0;

        /// Gets an estimation of the flexibility in the form of a homogeneous matrix
        virtual Matrix4 getFlexibility()=0;



    protected:



    };

} //namespace flexibilityEstimation
}//namespace stateobservation

#endif //FLEXIBILITY_ESTIMATION_FLEXIBILITY_ESTIMATOR_BASE_HPP
