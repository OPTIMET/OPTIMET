/** @file constants.h
 * Header file containing constants.
 * Header file with common constants that can be included wherever needed.
 * Constants are defined in constants.cpp as extern to prevent double
 * memory allocation. All constants should be abbreviated with cons.
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <complex>

extern const double consPi;       /**< The PI constant. */
extern const double consEpsilon0; /**< The permittivity of vacuum in <SI>. */
extern const double consMu0;      /**< The permeability of vacuum in <SI>. */

extern const double consC; /**< The speed of light in vacuum in <SI>. */

extern const double consFrnmTom; /**< Conversion from nm to m in <SI>. */

extern const std::complex<double> consC1;  /**< The complex number 1. */
extern const std::complex<double> consCi;  /**< The complex number i. */
extern const std::complex<double> consCmi; /**< The complex number -i. */
extern const std::complex<double> consCm1; /**< The complex number -1. */
extern const std::complex<double> consC0;  /**< The complex number 0. */

extern const double errEpsilon; /**< The predefined epsilon error value for
                                   comparison of floats. */

#endif /* CONSTANTS_H_ */
