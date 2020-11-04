#ifndef OBSERVATIONTOOLSHPP
#define OBSERVATIONTOOLSHPP

#include <boost/random.hpp>
#include <Eigen/Core>

namespace stateObserver
{
namespace unitTesting
{
class Tools
{
public:
  // add White Gaussian Noise to a vector
  // having a given bias and standard deviation(std)
  static Eigen::MatrixXd getWGNoise(const Eigen::MatrixXd & std,
                                    const Eigen::MatrixXd & bias,
                                    Index rows,
                                    Index cols = 1);

protected:
  static boost::lagged_fibonacci1279 gen_;
};

#include "tools.hxx"
} // namespace unitTesting
} // namespace stateObserver

#endif // OBSERVATIONTOOLSHPP
