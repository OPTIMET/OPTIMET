// (C) University College London 2017
// This file is part of Optimet, licensed under the terms of the GNU Public License
//
// Optimet is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Optimet is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Optimet. If not, see <http://www.gnu.org/licenses/>.

#ifndef OPTIMET_BESSEL_H
#define OPTIMET_BESSEL_H

#include "constants.h"
#include <complex>
#include <iostream>
#include <tuple>
#include <vector>

extern "C" {
int zbesj_(const double *, const double *, const double *, const long int *, const long int *,
           double *, double *, long int *, long int *);
int zbesh_(const double *, const double *, const double *, const long int *, const long int *,
           const long int *, double *, double *, long int *, long int *);
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
template <BESSEL_TYPE BesselType, bool Scaling = false>
std::tuple<std::vector<std::complex<double>>, std::vector<std::complex<double>>>
bessel(const std::complex<double> &z, long int max_order) {
  // Calling FORTRAN functions from C/C++ expects the arguments to be pointers
  // to int/real which means they must be rvalues
  const double order = 0.5;
  const long int scaling = Scaling ? 2 : 1;
  const long int bessel_type = BesselType;

  std::vector<std::complex<double>> data(max_order + 1);
  std::vector<std::complex<double>> ddata(max_order + 1);

  if(std::abs(z) <= errEpsilon) {
    for(int i = 0; i <= max_order; i++)
      data[i] = ddata[i] = std::complex<double>(0.0, 0.0);
    if(BesselType == Bessel)
      data[0] = std::complex<double>(1, 0);
  } else {
    const double zr = z.real();
    const double zi = z.imag();

    // Return vectors for real and imaginary parts
    // (+1 for zeroth order, +1 for derivative)
    const long int size = max_order + 2;
    std::vector<double> cyr(size);
    std::vector<double> cyi(size);

    long int zeroUnderflow, ierr;

    if(BesselType == Bessel)
      // Calculate the Bessel function of the first kind
      zbesj_(&zr, &zi, &order, &scaling, &size, cyr.data(), cyi.data(), &zeroUnderflow, &ierr);
    else
      // Calculate the Hankel function of the first or second kind
      zbesh_(&zr, &zi, &order, &scaling, &bessel_type, &size, cyr.data(), cyi.data(),
             &zeroUnderflow, &ierr);

    switch(ierr) {
    case 0:
      break;
    case 1:
      throw std::runtime_error("Incorrect input to Henkel/Bessel functions");
      break;
    case 2:
      throw std::runtime_error("Overflow when computing to Henkel/Bessel functions");
      break;
    case 3:
      throw std::runtime_error(
          "Loss of significance: result is less than half of machine accuracy");
      break;
    case 4:
      throw std::runtime_error("Complete loss of significance");
      break;
    case 5:
      throw std::runtime_error("Algorithmic conditions not met");
      break;
    default:
      throw std::runtime_error("Unknown error when computing Henkel/Bessel functions");
      break;
    }

    // Assemble the direct functions
    const std::complex<double> r = std::sqrt(consPi / (2.0 * z));
    for(int i = 0; i <= max_order; i++)
      data[i] = r * std::complex<double>(cyr[i], cyi[i]);

    // Assemble the derivative functions
    for(int i = 0; i < max_order; i++)
      ddata[i] = consCm1 * data[i + 1] + ((double)i / z) * data[i];

    // The last derivative
    ddata[max_order] = consCm1 * r * std::complex<double>(cyr[max_order + 1], cyi[max_order + 1]) +
                       ((double)max_order / z) * data[max_order];
  }

   if ( (z.imag() == 0.0) && (BesselType == Bessel) ){   
   // corrects the error if the imaginary part of the argument is exactly zero
   for (int i = 0; i <= max_order; i++){
   data[i]=data[i].real();
   ddata[i]=ddata[i].real();
   }
  
  }

  return std::make_tuple(data, ddata);
}

template <BESSEL_TYPE BesselType, bool Scaling = false>
std::vector<std::complex<double>>
bessel3der(const std::complex<double> &z, long int max_order) {

  // Calling FORTRAN functions from C/C++ expects the arguments to be pointers
  // to int/real which means they must be rvalues
  const double order = 0.5;
  const long int scaling = Scaling ? 2 : 1;
  const long int bessel_type = BesselType;
  
  

  std::vector<std::complex<double>> data(max_order + 1);
  std::vector<std::complex<double>> ddata(max_order + 1);
  std::vector<std::complex<double>> dddata(max_order + 1);
  
  if(std::abs(z) <= errEpsilon) {
    for(int i = 0; i <= max_order; i++){
      data[i] = ddata[i] = std::complex<double>(0.0, 0.0);
      dddata[i] = std::complex<double>(0.0, 0.0);
      }
    if(BesselType == Bessel)
      
      data[0] = std::complex<double>(1, 0);
  } else {
    const double zr = z.real();
    const double zi = z.imag();

    // Return vectors for real and imaginary parts
    // (+1 for zeroth order, +1 for derivative)
    const long int size = max_order + 2;
    std::vector<double> cyr(size);
    std::vector<double> cyi(size);

     long int zeroUnderflow, ierr;
     
  
    if(BesselType == Bessel)
      // Calculate the Bessel function of the first kind
      zbesj_(&zr, &zi, &order, &scaling, &size, cyr.data(), cyi.data(), &zeroUnderflow, &ierr);
    else
      // Calculate the Hankel function of the first or second kind
      zbesh_(&zr, &zi, &order, &scaling, &bessel_type, &size, cyr.data(), cyi.data(),
             &zeroUnderflow, &ierr);
             
    switch(ierr) {
    case 0:
      break;
    case 1:
      throw std::runtime_error("Incorrect input to Henkel/Bessel functions");
      break;
    case 2:
      throw std::runtime_error("Overflow when computing to Henkel/Bessel functions");
      break;
    case 3:
      throw std::runtime_error(
          "Loss of significance: result is less than half of machine accuracy");
      break;
    case 4:
      throw std::runtime_error("Complete loss of significance");
      break;
    case 5:
      throw std::runtime_error("Algorithmic conditions not met");
      break;
    default:
      throw std::runtime_error("Unknown error when computing Henkel/Bessel functions");
      break;
    }
    
     

    // Assemble the direct functions
    const std::complex<double> r = std::sqrt(consPi / (2.0 * z));
    for(int i = 0; i <= max_order; i++)
      data[i] = r * std::complex<double>(cyr[i], cyi[i]);

    // Assemble the first derivative functions
    for(int i = 0; i < max_order; i++)
      ddata[i] = consCm1 * data[i + 1] + ((double)i / z) * data[i];

    // The last derivative
    ddata[max_order] = consCm1 * r * std::complex<double>(cyr[max_order + 1], cyi[max_order + 1]) +
                       ((double)max_order / z) * data[max_order];
                       
   // Assemble the second derivative functions
    for(int i = 1; i <= max_order; i++)
      dddata[i] = consCm1 * ((double)i / z) * ddata[i] + ((double)i / std::pow(z , 2.0)) * data[i] + ddata[i - 1];
                       
                       
  }

  if ( (z.imag() == 0.0) && (BesselType == Bessel) ){
   // corrects the error if the imaginary part of the argument is exactly zero
   for (int i = 1; i <= max_order; i++){
   dddata[i]=dddata[i].real();
   }
  
  }

  return dddata;
}



inline std::tuple<std::vector<std::complex<double>>, std::vector<std::complex<double>>>
bessel(const std::complex<double> &z, enum BESSEL_TYPE besselType, bool scale, long int nMax) {
  switch(besselType) {
  case Bessel:
    return (scale) ? bessel<Bessel, true>(z, nMax) : bessel<Bessel, false>(z, nMax);
  case Hankel1:
    return (scale) ? bessel<Hankel1, true>(z, nMax) : bessel<Hankel1, false>(z, nMax);
  case Hankel2:
    return (scale) ? bessel<Hankel2, true>(z, nMax) : bessel<Hankel2, false>(z, nMax);
  }
}

} // namespace optimet

#endif /* OPTIMET_BESSEL_H */
