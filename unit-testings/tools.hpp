#ifndef OBSERVATIONTOOLSHPP
#define OBSERVATIONTOOLSHPP

#include <Eigen/Core>
#include <boost/random.hpp>

namespace stateObserver
{
    namespace unitTesting
    {
        class Tools
        {
        public:

            //add White Gaussian Noise to a vector
            //having a given bias and standard deviation(std)
            static Eigen::MatrixXd getWGNoise( const Eigen::MatrixXd& std, const Eigen::MatrixXd& bias, size_t rows, size_t cols=1);

        protected:
            static boost::lagged_fibonacci1279 gen_;

        };

#       include "tools.hxx"
    }
}



#endif //OBSERVATIONTOOLSHPP
