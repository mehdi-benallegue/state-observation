/**
 * \file      probability-law-simulation.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief
 *
 *
 *
 */


#ifndef SENSORSSIMULATIONPROBABILITYLAWSIMULATIONHPP
#define SENSORSSIMULATIONPROBABILITYLAWSIMULATIONHPP

#include <boost/random.hpp>

#include <state-observation/tools/definitions.hpp>


namespace stateObservation
{
    namespace tools
    {
        class ProbabilityLawSimulation
        {
        public:

            ///gets White Gaussian Noise
            ///having a given bias and standard deviation(std)
            static Matrix getGaussianVector( const Matrix & std, const Matrix & bias,
                unsigned rows, unsigned cols=1);

        protected:
            static boost::lagged_fibonacci1279 gen_;

        };

    }
}



#endif //SENSORSSIMULATIONPROBABILITYLAWSIMULATIONHPP
