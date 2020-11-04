/**
 * \file      magnetic-field.hpp
 * \author    Joseph Mirabel
 * \date      2015
 * \brief     Implements the magnetic field algorithm
 *
 * \details
 *
 *
 */



#ifndef SENSORALGORITHMSMAGNETICFIELDHPP
#define SENSORALGORITHMSMAGNETICFIELDHPP

#include <state-observation/api.h>
#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{
    namespace algorithm
    {
        /**
         * \class MagneticField
         * \brief Implements the measurements given by an magnetometer.
         *
         */

        class STATE_OBSERVATION_DLLAPI MagneticField
        {
        public:

            MagneticField();

            ///virtual destructor
            virtual ~MagneticField(){}

            ///The magnetic field measurement in the local frame represented by the orientation Matrix
            Vector3 magneticFieldMeasure(const Matrix3 & orientation) const;

        private:
            Vector3 earthLocalMagneticField_;
        };
    }

}

#endif //SENSORALGORITHMSLINEARACCELERATIONHPP
