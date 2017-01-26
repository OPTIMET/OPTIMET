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

/**
 * @brief The ElectroMagnetic class implements EM properties for scatterers.
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
   * @see init_r()
   */
  ElectroMagnetic(std::complex<double> epsilon_r_, std::complex<double> mu_r_);

  /**
   * Default destructor for ElectroMagnetic.
   */
  virtual ~ElectroMagnetic();
  // Model variables
  int modelType; /**< Mode being used. Currently supports 0 - fixed, 1 -
                    Sellmeier, 2 - Drude. */
  double B1;
  double C1;
  double B2;
  double C2;
  double B3;
  double C3;
  double B4;
  double C4;
  double B5;
  double C5;
  double lambda;
  std::complex<double> plasma_freq;
  std::complex<double> damping_freq;

  // Fundamental Frequency variables
  std::complex<double> epsilon;   /**< The absolute electric permittivity. */
  std::complex<double> mu;        /**< The absolute magnetic permeability. */
  std::complex<double> epsilon_r; /**< The relative electric permittivity. */
  std::complex<double> mu_r;      /**< The relative magnetic permeability. */

  // Second Harmonic variables
  std::complex<double> a_SH; /**< The a Non-Linear coefficient (cfn13arx). */
  std::complex<double> b_SH; /**< The b Non-Linear coefficient (cfn13arx). */
  std::complex<double> d_SH; /**< The d Non-Linear coefficient (cfn13arx). */
  std::complex<double>
      epsilon_SH; /**< The absolute second harmonic electric permittivity. */
  std::complex<double>
      mu_SH; /**< The absolute second harmonic magnetic permeability. */
  std::complex<double>
      epsilon_r_SH; /**< The relative second harmonic electric permittivity. */
  std::complex<double>
      mu_r_SH; /**< The relative second harmonic magnetic permeability. */

  /**
   * Initialization function for ElectroMagnetic.
   * Initializes epsilon and mu with ABSOLUTE values. Calculates relative
   * values. Set non-linear values to Bachelier values.
   * @param epsilon_ the complex absolute value for permittivity.
   * @param mu_ the complex absolute value for permeability.
   * @see ElectroMagnetic()
   * @see init_r()
   */
  void init(std::complex<double> epsilon_, std::complex<double> mu_);

  /**
   * Initialization function for ElectroMagnetic.
   * Initializes epsilon_r and mu_r with RELATIVE values. Calculate absolute
   * values. Set non-linear values to Bachelier values.
   * @param epsilon_r_ the complex relative value for permittivity.
   * @param mu_r_ the complex relative value for permeability.
   * @see ElectroMagnetic()
   * @see init()
   */
  void init_r(std::complex<double> epsilon_r_, std::complex<double> mu_r_);

  /**
   * Initialization function for Electromagnetic using the Sellmeier equation.
   * DOES NOT CALCULATE THE VALUES! USE UPDATE FOR THAT!
   * @param B1_ the Sellmeier B1 parameter.
   * @param C1_ the Sellmeier C1 parameter.
   * @param B2_ the Sellmeier B2 parameter.
   * @param C2_ the Sellmeier C2 parameter.
   * @param B3_ the Sellmeier B3 parameter.
   * @param C3_ the Sellmeier C3 parameter.
   * @param mu_r_ the magnetic permeabilitty (relative).
   * @param lambda_ the simulation wavelength.
   */
  void initSellmeier_r(double B1_, double C1_, double B2_, double C2_,
                       double B3_, double C3_, double B4_, double C4_,
                       double B5_, double C5_, std::complex<double> mu_r_);

  void initDrudeModel_r(std::complex<double> plasma_frequency_,
                        std::complex<double> damping_frequency_,
                        std::complex<double> mu_r_);

  /**
   * Populates the Sellmeier coefficients.
   */
  void populateSellmeier();

  void populateDrudeModel();

  /**
   * Updates the ElectroMagnetic object to a new wavelength.
   * @param lambda_ the new wavelength.
   */
  void update(double lambda_);
};

#endif /* ELECTROMAGNETIC_H_ */
