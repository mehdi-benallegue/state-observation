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

double ProbabilityLawSimulation::getUniformScalar(double min, double max)
{
  std::uniform_real_distribution<double> g(min, max);
  return g(gen_);
}
} // namespace tools

} // namespace stateObservation
