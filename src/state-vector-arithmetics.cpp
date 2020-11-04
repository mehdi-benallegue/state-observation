#include <state-observation/tools/state-vector-arithmetics.hpp>

namespace stateObservation
{
void detail::defaultSum(const Vector &stateVector, const Vector &tangentVector, Vector &result)
{
    result.noalias() = stateVector + tangentVector;
}

void detail::defaultDifference(const Vector &stateVector1, const Vector &stateVector2, Vector &difference)
{
    difference.noalias() = stateVector1 - stateVector2;
}

void StateVectorArithmetics::stateSum(const Vector &stateVector, const Vector &tangentVector, Vector &sum)
{
    detail::defaultSum(stateVector, tangentVector, sum);
}

void StateVectorArithmetics::stateDifference(const Vector &stateVector1, const Vector &stateVector2, Vector &difference)
{
    detail::defaultDifference(stateVector1, stateVector2, difference);
}

void StateVectorArithmetics::measurementDifference(const Vector &measureVector1, const Vector &measureVector2, Vector &difference)
{
    detail::defaultDifference(measureVector1, measureVector2, difference);
}

} // namespace stateObservation
