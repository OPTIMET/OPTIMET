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

#ifndef COUPLING_H_
#define COUPLING_H_

#include "Tools.h"
#include "Types.h"

namespace optimet {
/**
 * The Coupling class implements the calculation of the A and B coupling
 * coefficients.
 * @warning Do not use without initialization!
 */
struct Coupling {
  Matrix<t_complex> diagonal;    /**< The A_nmkl coefficients. */
  Matrix<t_complex> offdiagonal; /**< The B_nmkl coefficients. */

  /**
   * Initialization constructor for the Coupling class.
   * @param relR_ the relative spherical vector between two spheres.
   * @param waveK_ the complex wave number.
   * @param regular_ the regular flag.
   * @param nMax_ the maximum value of the n iterator.
   */
  Coupling(Spherical<t_real> relR_, t_complex waveK_, t_uint nMax_, bool regular_ = true);
};
}

#endif /* COUPLING_H_ */
