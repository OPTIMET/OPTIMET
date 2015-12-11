#ifndef OPTIMET_BESSEL_H
#define OPTIMET_BESSEL_H

#include <vector>
#include <tuple>
#include <complex>

#include "constants.h"

extern "C" {
int zbesj_(const double *, const double *, const double *, const long int *,
           const long int *, double *, double *, long int *, long int *);
int zbesh_(const double *, const double *, const double *, const long int *,
           const long int *, const long int *, double *, double *, long int *,
           long int *);
}

namespace optimet {

enum BESSEL_TYPE { Bessel = 0, Hankel1 = 1, Hankel2 = 2 };

/*!
 * The bessel function implements the Spherical Bessel and Hankel functions
 * and their derivatives, calculated from the zeroth order up to the maximum
 * order.
 *
 * \tparam BesselType   the type of function:
 *                        \c 0 - Bessel,
 *                        \c 1 - Hankel (first kind),
 *                        \c 2 - Hankel (second kind)
 * \tparam ScalingType  the scaling type:
 *                        \c 0 - unscaled,
 *                        \c 1 - scaled
 *
 * \param [in] z          the argument for the Bessel function
 * \param [in] max_order  the maximum order of functions to calculate
 *
 * \return a tuple containing the values of the spherical bessel and hankel
 *           functions in the first element, and their derivatives in the second
 *
 * \warning The zeroth order derivative (never used) is not accurate!
 */
template <BESSEL_TYPE BesselType, bool ScalingType>
std::tuple<std::vector<std::complex<double>>, std::vector<std::complex<double>>>
bessel(const std::complex<double> &z, long int max_order) {
  static constexpr double order = 0.5;
  static constexpr long int scaling_type = ScalingType + 1;

  std::vector<std::complex<double>> data(max_order + 1);
  std::vector<std::complex<double>> ddata(max_order + 1);

  if (std::abs(z) <= errEpsilon) {
    for (int i = 0; i <= max_order; i++)
      data[i] = ddata[i] = std::complex<double>(0.0, 0.0);
  } else {
    const double zr = z.real();
    const double zi = z.imag();

    // Return vectors for real and imaginary parts
    // (+1 for zeroth order, +1 for derivative)
    std::vector<double> cyr(max_order + 2);
    std::vector<double> cyi(max_order + 2);

    long int zeroUnderflow, ierr;

    if (BesselType == Bessel)
      // Calculate the Bessel function of the first kind
      zbesj_(&zr, &zi, &order, &scaling_type, &cyr.size(), cyr.data(),
             cyi.data(), &zeroUnderflow, &ierr);
    else
      // Calculate the Hankel function of the first or second kind
      zbesh_(&zr, &zi, &order, &scaling_type, &bessel_type, &cyr.size(),
             cyr.data(), cyi.data(), &zeroUnderflow, &ierr);

    if (ierr != 0)
      throw std::range_error("Error computing Bessel/Hankel functions");

    // Assemble the direct functions
    const double r = std::sqrt(consPi / (2.0 * z));
    for (int i = 0; i <= max_order; i++)
      data[i] = r * std::complex<double>(cyr[i], cyi[i]);

    // Assemble the derivative functions
    for (int i = 0; i < max_order; i++)
      ddata[i] = consCm1 * data[i + 1] + ((double)i / z) * data[i];

    // The last derivative
    ddata[maxOrder] = consCm1 * r * std::complex<double>(cyr[max_order + 1],
                                                         cyi[max_order + 1]) +
                      ((double)max_order / z) * data[max_order];
  }

  return std::make_tuple(data, ddata);
}

} // namespace optimet

#endif /* OPTIMET_BESSEL_H */
