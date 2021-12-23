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

#ifndef ELECTROMAGNETIC_H_
#define ELECTROMAGNETIC_H_

#include <complex>
#include <vector>

/**
 * @brief The ElectroMagnetic class implements EM properties of scatterers.
 * The ElectroMagnetic class can be used to create a collection of related
 * electromagnetic properties for scatterers which are either directly
 * initialized or calculated from existing data.
 */

class ElectroMagnetic {
private:
  //
public:
  /**
   * Default constructor for ElectroMagnetic.
   * Initializes epsilon, mu, epsilon_r and mu_r with vacuum values.
   */
  ElectroMagnetic();

  /**
   * Initialization constructor for ElectroMagnetic.
   * Initializes epsilon_r and mu_r with RELATIVE values by calling init_r.
   * @param epsilon_r_ the complex relative value for permittivity.
   * @param mu_r_ the complex relative value for permeability.
   * other parameters are values for the second harmonic case
   * @see init_r()
   */
  ElectroMagnetic(std::complex<double> epsilon_r_, std::complex<double> mu_r_, std::complex<double> epsilon_r_SH, std::complex<double> ksippp, std::complex<double> ksiparppar, std::complex<double> gamma);

  
  //Default destructor for ElectroMagnetic.
  virtual ~ElectroMagnetic();
  // Model variables
  int modelType; /**< Mode being used. Currently supports 0-fixed, 3-Hydrodynamic_Gold_Johnson, 4-Silicon_Schinke */
  
  double lambda;
 

// Schinke data for lambda in range 0.5 - 1.45um
   std::vector<double> refInd, Extcoeff;

  // Fundamental Frequency variables
  std::complex<double> epsilon;   /**< The absolute electric permittivity. */
  std::complex<double> mu;        /**< The absolute magnetic permeability. */
  std::complex<double> epsilon_r; /**< The relative electric permittivity. */
  std::complex<double> mu_r;      /**< The relative magnetic permeability. */

  // Second Harmonic variables, Rudnick-Stern parameters
  std::complex<double> a_SH; 
  std::complex<double> b_SH; 
  std::complex<double> d_SH; 
  std::complex<double>
      epsilon_SH; /**< The absolute second harmonic electric permittivity. */
  std::complex<double>
      mu_SH; /**< The absolute second harmonic magnetic permeability. */
  std::complex<double>
      epsilon_r_SH; /**< The relative second harmonic electric permittivity. */
  std::complex<double>
      mu_r_SH; /**< The relative second harmonic magnetic permeability. */
  
  // surface tensor components    
  std::complex<double> ksippp;    
  std::complex<double> ksiparppar;

  // bulk tensor component
  std::complex<double> gamma;  
  
  /**
   * Initialization function for ElectroMagnetic.
   * Initializes epsilon_r and mu_r with RELATIVE values. Calculate absolute
   * values. Set non-linear values to Bachelier values.
   * @param epsilon_r_ the complex relative value for permittivity.
   * @param mu_r_ the complex relative value for permeability.
   * @see ElectroMagnetic()
   * @see init()
   */
  void init_r(std::complex<double> epsilon_r_, std::complex<double> mu_r_, std::complex<double> epsilon_r_SH, std::complex<double> ksippp, std::complex<double> ksiparppar, std::complex<double> gamma);

                        
 void initHydrodynamicModel_r(std::complex<double> a_SH,
                        std::complex<double> b_SH,
                        std::complex<double> d_SH,
                        std::complex<double> mu_r_);                       

  void initSiliconModel_r(std::complex<double> mu_r_);
  
  void populateHydrodynamicModel();

  void populateSiliconModel();

  /**
   * Updates the ElectroMagnetic object to a new wavelength.
   * @param lambda_ the new wavelength.
   */
  void update(double lambda_);
};

#endif /* ELECTROMAGNETIC_H_ */
