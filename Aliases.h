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

#ifndef ALIASES_H_
#define ALIASES_H_

/** @file Aliases.h
 * Header file containing common aliases used in the code.
 * Aliases are implemented as preprocessor defines and start with O3D.
 * Alias groups:
 *  1 - solver methods
 *  8 - geometry types
 *  9 - output grids
 */

/** @brief Cartesian Regular grid. */
#define O3DCartesianRegular 911

/** @brief Spiral structure. */
#define O3DGeometrySpiral 311

/** @brief Mischenko 1996 solver (Direct). */
#define O3DSolverDirect 111

/** @brief Stout 2002 solver (Indirect). */
#define O3DSolverIndirect 112

#endif /* ALIASES_H_ */
