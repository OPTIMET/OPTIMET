#ifndef COUPLING_H_
#define COUPLING_H_

#include "Tools.h"
#include "Spherical.h"
#include "Types.h"

namespace optimet {
/**
 * The Coupling class implements the calculation of the A and B coupling
 * coefficients.
 * @warning Do not use without initialization!
 */
struct Coupling {
  Matrix<t_complex> diagonal;    /**< The A_nmkl coefficients. */
  Matrix<t_complex> offdiagonal; /**< The B_nmkl coefficients. */

  /**
   * Initialization constructor for the Coupling class.
   * @param relR_ the relative spherical vector between two spheres.
   * @param waveK_ the complex wave number.
   * @param regular_ the regular flag.
   * @param nMax_ the maximum value of the n iterator.
   */
  Coupling(Spherical<t_real> relR_, t_complex waveK_, t_uint nMax_, bool regular_ = true);
};
}

#endif /* COUPLING_H_ */
