/**
 * \file      linear-acceleration.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief     Implements the accelerometer algorithm
 *
 * \details
 *
 *
 */



#ifndef SENSORALGORITHMSLINEARACCELERATIONHPP
#define SENSORALGORITHMSLINEARACCELERATIONHPP

#include <state-observation/api.h>
#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{
    namespace algorithm
    {
        /**
         * \class LinearAcceleration
         * \brief Implements the measurements given by an accelerometer.
         *
         */

        class STATE_OBSERVATION_DLLAPI LinearAcceleration
        {
        public:
            ///virtual destructor
            virtual ~LinearAcceleration(){}

            ///The acceleration measurement in the local frame represented by the orientation Matrix
            Vector3 accelerationMeasure(const Vector3 & acceleration, const Matrix3 & orientation) const;

        protected:

        };
    }

}

#endif //SENSORALGORITHMSLINEARACCELERATIONHPP
