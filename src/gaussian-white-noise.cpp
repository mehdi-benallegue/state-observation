#include <boost/assert.hpp>
#include <Eigen/Cholesky>
#include <state-observation/noise/gaussian-white-noise.hpp>
#include <state-observation/tools/probability-law-simulation.hpp>

namespace stateObservation
{

GaussianWhiteNoise::GaussianWhiteNoise(Index dimension)
: dim_(dimension), std_(Matrix::Identity(dimension, dimension)), bias_(Vector::Zero(dimension, 1)),
  sum_(detail::defaultSum)
{
}

GaussianWhiteNoise::GaussianWhiteNoise() : dim_(0), sum_(detail::defaultSum) {}

Vector GaussianWhiteNoise::getNoisy(const Vector & v)
{
  checkVector_(v);

  sum_(v, tools::ProbabilityLawSimulation::getGaussianVector(std_, bias_, dim_), noisy_);

  return noisy_;
}

void GaussianWhiteNoise::setStandardDeviation(const Matrix & std)
{
  checkMatrix_(std);
  std_ = std;
}

void GaussianWhiteNoise::setCovarianceMatrix(const Matrix & cov)
{
  checkMatrix_(cov);
  Matrix L(cov.llt().matrixL());
  std_ = L;
}

void GaussianWhiteNoise::setBias(const Vector & bias)
{
  checkVector_(bias);
  bias_ = bias;
}

Index GaussianWhiteNoise::getDimension() const
{
  return dim_;
}

void GaussianWhiteNoise::setDimension(Index dim)
{
  dim_ = dim;
  bias_ = Vector::Zero(dim, 1);
  std_ = Matrix::Identity(dim, dim);
}

void GaussianWhiteNoise::checkMatrix_(const Matrix & m) const
{
  (void)m; // avoid warning
  BOOST_ASSERT(m.rows() == dim_ && m.cols() == dim_ && "ERROR: Matrix incorrecly dimemsioned");
}

void GaussianWhiteNoise::checkVector_(const Vector & v) const
{
  (void)v; // avoid warning
  BOOST_ASSERT(v.rows() == dim_ && v.cols() == 1 && "ERROR: Vector incorrecly dimemsioned");
}

void GaussianWhiteNoise::setSumFunction(void (*sum)(const Vector & stateVector,
                                                    const Vector & tangentVector,
                                                    Vector & result))
{
  sum_ = sum;
}
} // namespace stateObservation
