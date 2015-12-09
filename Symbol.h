/**
 * The Symbol class contains special symbol routines.
 * Symbol is released under the GSL. To view a copy
 * of the licence, look in the documentation.
 */

#ifndef OPTIMET_SYMBOL_H
#define OPTIMET_SYMBOL_H

#include "Scatterer.h"
#include "ElectroMagnetic.h"

namespace optimet {
namespace symbol {

/**
 * Calculates the Wigner 3j symbol: \n
 * (j1 j2 j3 \n
 *  m1 m2 m3).
 * @param j1 coefficient.
 * @param j2 coefficient.
 * @param j3 coefficient.
 * @param m1 coefficient.
 * @param m2 coefficient.
 * @param m3 coefficient.
 * @return the Wigner 3j symbol.
 */
double Wigner3j(int j1, int j2, int j3, int m1, int m2, int m3);

std::complex<double> up_mn(int m, int n, int nMax,
                           const std::complex<double> & cmn_1,
                           const std::complex<double> & dmn_1, double omega,
                           const Scatterer & object,
                           const ElectroMagnetic & bground);

std::complex<double> vp_mn(int m, int n, int nMax,
                           const std::complex<double> & cmn_1,
                           const std::complex<double> & dmn_1, double omega,
                           const Scatterer & object,
                           const ElectroMagnetic & bground);

std::complex<double> upp_mn(int m, int n, int nMax,
                            const std::complex<double> & cmn_1,
                            const std::complex<double> & dmn_1, double omega,
                            const Scatterer & object);

};

} // namespace symbol
} // namespace optimet

#endif /* OPTIMET_SYMBOL_H */
