#include "AJ_AuxFuns.h"

#include "AuxCoefficients.h"
#include "Bessel.h"
#include "CompoundIterator.h"
#include "constants.h"

#include <cmath>
#include <iostream>
#include <cstdlib>

namespace optimet {

// computes all corresponding Spherical Harmonics values given (nMax, m)
std::vector<std::complex<double>> compute_Yn_m(const Spherical<double> &R,
                                               const int nMax, const int m) {
  std::vector<std::complex<double>> Ynm(nMax + 1);

  std::vector<double> dn(nMax + 1);

  std::vector<double> Wigner(nMax + 1), dWigner(nMax + 1);

  AuxCoefficients auxCoefficients;
  auxCoefficients.compute_dn(nMax, dn.data());
  auxCoefficients.VIGdVIG(nMax, m, R, Wigner.data(), dWigner.data());

  const double dm =
      std::pow(-1.0, static_cast<double>(m)); // Legendre to Wigner function
  const std::complex<double> exp_imphi(
      std::cos(static_cast<double>(m) * R.phi),
      std::sin(static_cast<double>(m) * R.phi));

  for (int i = 0; i <= nMax; i++) {
    // obtain Spherical hormonics - Ynm
    const int n = i;
    double d_n = static_cast<double>(i);

    const double d_temp = dm * dn[i] * std::sqrt(d_n * (d_n + 1.0));

    if (m == 0 && n == 0)
      Ynm[i] = std::sqrt(1.0 / (4.0 * consPi));
    else
      Ynm[i] = d_temp * exp_imphi * Wigner[i];
  }

  return Ynm;
}

// computes all corresponding Spherical Harmonics values given (nMax) -
// m=[-nMax, nMax]
std::vector<std::complex<double>> compute_Yp(const Spherical<double> &R,
                                             const int nMax) {

  std::vector<std::complex<double>> dataYp(CompoundIterator::max(nMax) + 1);

  // (n,m)=(0,0)
  std::vector<std::complex<double>> Yn_m = compute_Yn_m(R, nMax, 0);
  dataYp[CompoundIterator::max(nMax)] = Yn_m[0];

  // all other elements compute and map into CompoundIterator format
  for (CompoundIterator q(nMax, nMax); q < q.max(nMax); ++q) {
    for (int n = std::abs(q.second); n <= nMax; n++) {
      if (n != 0) {
        Yn_m = compute_Yn_m(R, nMax, q.second);
        CompoundIterator p(n, q.second);
        dataYp[static_cast<int>(p)] = Yn_m[n];
      }
    }
  }

  return dataYp;
}

} // namespace optimet
