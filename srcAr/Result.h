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

#ifndef RESULT_H_
#define RESULT_H_

#include "CompoundIterator.h"
#include "Excitation.h"
#include "Geometry.h"
#include "OutputGrid.h"
#include "Spherical.h"
#include "SphericalP.h"
#include <vector>
#ifdef OPTIMET_MPI
#include <mpi.h>
#endif
#include <complex>

namespace optimet {
/**
 * The Result class is used to provide post simulation
 * output functions including field profiles, absorption
 * and scattering cross sections, and scattering
 * coefficients.
 */

class Result {
private:
  std::shared_ptr<Geometry> geometry;     /**< Pointer to the Geometry. */
  std::shared_ptr<Excitation> excitation; /**< Pointer to the Excitation. */
  std::complex<double> waveK;             /**< The std::complex wave number. */
  Result *result_FF;                      /**< The Fundamental Frequency results vector. */
  Result *result_SH;                      /**< The Second Harmonic results vector. */
  //! Maximum nMax
  optimet::t_uint nMax;
  optimet::t_uint nMaxS;
public:
  Vector<t_complex> scatter_coef;   /**< The scattering coefficients. */
  Vector<t_complex> internal_coef;  /**< The internal coefficients. */
  Vector<t_complex> scatter_coef_SH;   /**< The scattering coefficients. second harmonic */
  Vector<t_complex> internal_coef_SH;  /**< The internal coefficients. second harmonic */
  std::vector<double *> CLGcoeff; // Coefficients needed for SH source calculations
  /**
   * Initialization constructor for the Result class.
   * Fundamental Frequency version.
   * @param geometry_ the pointer to the geometry.
   * @param excitation_ the pointer to the excitation.
   */
  Result(std::shared_ptr<Geometry> geometry_, std::shared_ptr<Excitation> xcitation_);

  /**
   * Initialization constructor for the Result class.
   * Second Harmonic version.
   * @param geometry_ the pointer to the geometry.
   * @param excitation_ the pointer to the excitation.
   * @param result_FF_ the pointer to the Fundamental Frequency results.
   */
  Result(std::shared_ptr<Geometry> geometry_, std::shared_ptr<Excitation> excitation_,
         Result *result_FF_);

  /**
   * Default destructor for the Result class.
   */
  virtual ~Result(){};

  /**
   * Initialization method for the Result class.
   * Fundamental Frequency version.
   * @param geometry_ the pointer to the geometry.
   * @param excitation_ the pointer to the excitation.
   */
  void
  init(std::shared_ptr<Geometry> geometry_, std::shared_ptr<Excitation> excitation_);

  /**
   * Update method for the Result class.
   * @param geometry_ the pointer to the geometry.
   * @param excitation_ the pointer to the excitation.
   */
  void
  update(std::shared_ptr<Geometry> geometry_, std::shared_ptr<Excitation> excitation_);

  /**
   * Initialization constructor for the Result class.
   * Second Harmonic version.
   * @param geometry_ the pointer to the geometry.
   * @param excitation_ the pointer to the excitation.
   * @param result_FF_ the pointer to the Fundamental Frequency results.
   */
  void init(std::shared_ptr<Geometry> geometry_, std::shared_ptr<Excitation> excitation_,
            Result *result_FF_);


  /**
   * Returns the E and H fields at a given point.
   * @param R_ the coordinates of the point.
   * @param EField_ SphericalP vector that will store the E field.
   * @param HField_ SphericalP vector that will store the H field.
   * @param projection_ defines spherical (1) or cartesian (0) projection.
   */
  void getEHFields(Spherical<double> R_, SphericalP<std::complex<double>> &EField_FF,
       SphericalP<std::complex<double>> &HField_FF, SphericalP<std::complex<double>> &EField_SH,
       SphericalP<std::complex<double>> &HField_SH, bool projection_, std::complex<double> *coeffXmn, std::complex<double> *coeffXpl) const;
                   
  #ifdef OPTIMET_MPI
  //Returns the Extinction Cross Section for Fundamental Frequency.
  double getExtinctionCrossSection(int gran1, int gran2);
  
  // Scattering cross section for FF
  double getScatteringCrossSection(int gran1, int gran2);

  // Returns the Scattering Cross Section for SH Frequency.
  double getScatteringCrossSection_SH(int gran1, int gran2);  
 #endif

  /**
   * Populate a grid with E and H fields.
   * @param oEGrid_FF the OutputGrid object for the E fields, fundamental freq.
   * @param oHGrid_FF the OutputGrid object for the H fields, fundamental freq.
   * @param oEGrid_SH the OutputGrid object for the E fields, SH freq.
   * @param oHGrid_SH the OutputGrid object for the H fields, SH freq.
   * @param projection_ defines spherical (1) or cartesian (0) projection.
   * @return 0 if succesful, 1 otherwise.
   */
#ifdef OPTIMET_MPI
 int setFields(std::vector<double> &Rr, std::vector<double> &Rthe, 
std::vector<double> &Rphi, bool projection_, std::vector<double *> CLGcoeff);
#endif
};
}
#endif /* RESULT_H_ */
