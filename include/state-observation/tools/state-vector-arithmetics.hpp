#ifndef STATE_VECTOR_ARITHMETICS_HPP
#define STATE_VECTOR_ARITHMETICS_HPP

#include <state-observation/api.h>
#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{
namespace detail
{
void STATE_OBSERVATION_DLLAPI defaultSum(const Vector & stateVector, const Vector & tangentVector, Vector & result);

void STATE_OBSERVATION_DLLAPI defaultDifference(const Vector & stateVector1,
                                                const Vector & stateVector2,
                                                Vector & difference);

} // namespace detail

/**
 * \class  StateArithmetics
 * \brief
 *        This class is used to customize the way the difference between measurements,
 *        the state update function and the differentiation are performed.
 *         default is the usual natual arithmetics. overload any ohter one
 *
 */
class STATE_OBSERVATION_DLLAPI StateVectorArithmetics
{
public:
  virtual void stateSum(const Vector & stateVector, const Vector & tangentVector, Vector & sum);

  virtual void stateDifference(const Vector & stateVector1, const Vector & stateVector2, Vector & difference);

  virtual void measurementDifference(const Vector & measureVector1, const Vector & measureVector2, Vector & difference);
};
} // namespace stateObservation

#endif // STATE_VECTOR_ARITHMETICS_HPP
