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

#ifndef EXCITATION_H_
#define EXCITATION_H_

#include "Spherical.h"
#include "Types.h"
#include "constants.h"
#include <complex>
#include <memory>


namespace optimet {
/**
 * The Excitation class implements the incoming wave coefficients.
 * Possible Excitation types are:
 *  plane wave -- FF
 * SH sources - Polarization densities surface and bulk
 */
class Excitation {
public:
  //! The incoming wave values as (0, Einc_the, Einc_phi) in [V/m].
  SphericalP<std::complex<double>> Einc;
  Spherical<double> vKInc; /**< The incoming wavevector angular values as (R,
                              Kinc_the, Kinc_pphi). */
  long nMax;               /**< The maximum value for the n iterator. */
  unsigned long type;     
  bool SH_cond; // condition for the existence SH sources
  Vector<t_complex> dataIncAp; /**< The incoming wave a_n^m coefficients. */
  Vector<t_complex> dataIncBp; /**< The incoming wave b_n^m coefficients. */

  std::complex<double> waveK, bgcoef; /**< The incoming wave wavenumber. */
  
  
  /**
   * Initialization constructor for the Excitation class.
   * @param type_ the type of the wave.
   * @param Einc_ the incoming wave values.
   * @param waveKInc_ the incoming wavevector values.
   * @param nMax_ the maximum value of the n iterator.
   */
  Excitation(unsigned long type_, SphericalP<std::complex<double>> Einc_, bool SH_cond,
             Spherical<double> vKInc_, int nMax_, std::complex<double> bgcoeff);

  /**
   * Default destructor for the Excitation class.
   */
  virtual ~Excitation() {}

  /**
   * Update method for the Excitation class.
   * @param type_ the type of the wave.
   * @param Einc_ the incoming wave values.
   * @param waveKInc_ the incoming wavevector values.
   * @param nMax_ the maximum value of the n iterator.
   */
  void update(unsigned long type_, SphericalP<std::complex<double>> Einc_,
              Spherical<double> vKInc_, int nMax_);

  /**
   * Populate the a_n^m and b_n^m coefficients for the incoming wave.
   * Stored in dataIncAp and dataIncBp with CompoundIterator indexing.
   * @return 0 if successful, 1 otherwise.
   */
  int populate();

  /**
   * Returns the Incoming Local matrix converting the a,b coefficients into
   * local ones.
   * @param point_ the new origin of the coordinate system.
   * @param Inc_local_ the incoming local matrix.
   * @param nMax_ the maximum value of the n iterator.
   * @return 0 if successful, 1 otherwise.
   */
  int getIncLocal(Spherical<double> point_, std::complex<double> *Inc_local_,
                  int nMax_) const;
                                  
  /**
   * Updates the wavelength of the current excitation object to a new value.
   * @param lambda_ the new value of the wavelength.
   */
  void updateWavelength(double lambda_);

  //! The incoming wave wavelength
  t_real lambda() const { return 2 * constant::pi / std::real(vKInc.rrr); }
  //! The incoming wave frequency number (calculated)
  optimet::t_real wavenumber() const { return std::real(vKInc.rrr); }
  //! The incoming wave frequency (calculated)
  optimet::t_real omega() const { return constant::c * wavenumber(); }

};
}
#endif /* EXCITATION_H_ */
