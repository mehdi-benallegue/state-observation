#include <state-observation/tools/probability-law-simulation.hpp>

namespace stateObservation
{
namespace tools
{

std::random_device ProbabilityLawSimulation::rd_;
std::mt19937 ProbabilityLawSimulation::gen_{ProbabilityLawSimulation::rd_()};

double ProbabilityLawSimulation::getGaussianScalar(double std, double bias)
{
  std::normal_distribution<> g(bias, std);
  return g(gen_);
}

Matrix ProbabilityLawSimulation::getGaussianVector(const Matrix & std, const Matrix & bias, Index rows, Index cols)
{
  std::normal_distribution<> g(0, 1);
  Matrix ret = Matrix::Zero(rows, cols);
  for(Index i = 0; i < rows; ++i)
  {
    for(Index j = 0; j < cols; ++j) ret(i, j) = g(gen_);
  }
  ret = std * ret + bias;

  return ret;
}
} // namespace tools

} // namespace stateObservation
