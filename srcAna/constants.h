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

/** @file constants.h
 * Header file containing constants.
 * Header file with common constants that can be included wherever needed.
 * Constants are defined in constants.cpp as extern to prevent double
 * memory allocation. All constants should be abbreviated with cons.
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include "Types.h"

extern const optimet::t_real consPi;       /**< The PI constant. */
extern const optimet::t_real consEpsilon0; /**< The permittivity of vacuum in <SI>. */
extern const optimet::t_real consMu0;      /**< The permeability of vacuum in <SI>. */

extern const optimet::t_real consC; /**< The speed of light in vacuum in <SI>. */

extern const optimet::t_real consFrnmTom; /**< Conversion from nm to m in <SI>. */

extern const optimet::t_complex consC1;  /**< The complex number 1. */
extern const optimet::t_complex consCi;  /**< The complex number i. */
extern const optimet::t_complex consCmi; /**< The complex number -i. */
extern const optimet::t_complex consCm1; /**< The complex number -1. */
extern const optimet::t_complex consC0;  /**< The complex number 0. */

extern const optimet::t_real errEpsilon; /**< The predefined epsilon error value for
                                   comparison of floats. */
namespace optimet {
namespace constant {
extern const t_real pi;
extern const t_real epsilon0; /**< The permittivity of vacuum in <SI>. */
extern const t_real mu0;      /**< The permeability of vacuum in <SI>. */
extern const t_real c;        /**< The speed of light in vacuum in <SI>. */

extern const t_real from_nm_to_m; /**< Conversion from nm to m in <SI>. */

extern const t_complex i; /**< The complex number 1. */

extern const t_real tolerance;
}
}

#endif /* CONSTANTS_H_ */
