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
class Coupling {
public:
  t_uint nMax;            /**< The maximum value of the n iterator. */
  Spherical<t_real> relR; /**< Relative spherical vector between two spheres. */
  t_complex waveK;        /**< The complex wave number. */
  bool regular; /**< Specifies if we should be using the regular (J) functions.
                   1-yes, 0-no. */

  Matrix<t_complex> dataApq; /**< The A_nmkl coefficients. */
  Matrix<t_complex> dataBpq; /**< The B_nmkl coefficients. */

  /**
   * Initialization constructor for the Coupling class.
   * @param relR_ the relative spherical vector between two spheres.
   * @param waveK_ the complex wave number.
   * @param regular_ the regular flag.
   * @param nMax_ the maximum value of the n iterator.
   */
  Coupling(Spherical<t_real> relR_, t_complex waveK_, bool regular_,
           t_uint nMax_)
      // yep, regular is not regular because optimet
      : nMax(nMax_),
        relR(relR_),
        waveK(waveK_),
        regular(not regular_),
        dataApq(Matrix<t_complex>::Zero(Tools::iteratorMax(nMax),
                                        Tools::iteratorMax(nMax))),
        dataBpq(Matrix<t_complex>::Zero(Tools::iteratorMax(nMax),
                                        Tools::iteratorMax(nMax))) {}

  /**
   * Default destructor for the Coupling class.
   */
  virtual ~Coupling() {}

  /**
   * Populate the dataApq and dataBpq variables.
   * @return 0 if successful, 1 otherwise.
   */
  int populate();
};
}

#endif /* COUPLING_H_ */
