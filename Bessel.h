#ifndef OPTIMET_BESSEL_H
#define OPTIMET_BESSEL_H

#include <vector>
#include <tuple>
#include <complex>
#include <cmath>

#include "constants.h"

extern "C"
{
  int zbesj_(double *, double *, double *, long int *, long int *, double *, double *, long int *, long int *);
  int zbesh_(double *, double *, double *, long int *, long int *, long int *, double *, double *, long int *, long int *);
}

namespace optimet {

/*!
 * The Bessel function implements the Spherical Bessel and Hankel functions
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
 * \tparam MaxOrder     the maximum order of functions to calculate
 *
 * \param [in] z  the argument for the Bessel function
 *
 * \return a tuple containing the values of the spherical bessel and hankel
 *           functions in the first element, and their derivatives in the second
 *
 * \warning The zeroth order derivative (never used) is not accurate!
 */
template <long int BesselType, long int ScalingType, long int MaxOrder>
std::tuple<std::vector<std::complex<double>>, std::vector<std::complex<double>>>
bessel(const std::complex<double> & z) {

  static_assert(BesselType >= 0 && BesselType <= 2, "Wrong besselType");
  static_assert(MaxOrder >= 1, "Maximum order required for Bessel smaller than "
                               "1!");

  const double zr = z.real();
  const double zi = z.imag();
  const double order = 0.5;

  // Return vectors for real and imaginary parts
  // (+1 for zeroth order, +1 for derivative)
  std::vector<double> cyr(MaxOrder + 2);
  std::vector<double> cyi(MaxOrder + 2);

  long int zeroUnderflow, ierr;

  if (besselType == 0)
    // Calculate the Bessel function of the first kind
    zbesj_(&zr, &zi, &order, &scaling_type, &cyr.size(),
           cyr.data(), cyi.data(), &zeroUnderflow, &ierr);
  else
    // Calculate the Hankel function of the first or second kind
    zbesh_(&zr, &zi, &order, &scaling_type, &bessel_type, &cyr.size(),
           cyr.data(), cyi.data(), &zeroUnderflow, &ierr);

  if (ierr != 0)
    throw std::range_error("Error computing Bessel/Hankel functions");

  std::vector<std::complex<double>> data(max_order + 1);
  std::vector<std::complex<double>> ddata(max_order + 1);

  // Assemble the direct functions
  for(int i = 0; i <= MaxOrder; i++) {
    if (std::abs(z) <= errEpsilon)
      data[i] = std::complex<double>(0.0, 0.0);
    else
      data[i] = std::sqrt(consPi / (2.0 * z)) *
                std::complex<double>(cyr[i], cyi[i]);
  }

  //The last derivative
  ddata[maxOrder] = consCm1 * std::sqrt(consPi / (2.0 * z)) *
                    std::complex<double>(cyr[MaxOrder + 1],
                                         cyi[MaxOrder + 1]) +
                    ((double)MaxOrder / z) * data[MaxOrder];

  for (int i = 0; i < MaxOrder; i++)
    ddata[i] = consCm1 * data[i + 1] + ((double)i / z) * data[i];

  return std::make_tuple(data, ddata);
}

} // namespace optimet

#endif /* OPTIMET_BESSEL_H */
