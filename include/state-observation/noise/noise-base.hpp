/**
 * \file      noise-base.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief
 *
 *
 *
 */

#ifndef SENSORSIMULATIONNOISEBASEHPP
#define SENSORSIMULATIONNOISEBASEHPP

#include <state-observation/api.h>
#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{

/**
 * \class  NoiseBase
 * \brief
 *
 * \details
 *
 */

class STATE_OBSERVATION_DLLAPI NoiseBase
{
public:
  /// Virtual destructor
  virtual ~NoiseBase() {}

  /// The method to overload to produce the noisy version of a given vector
  virtual Vector getNoisy(const Vector &) = 0;

protected:
};

} // namespace stateObservation

#endif // SENSORSIMULATIONNOISEBASEHPP
