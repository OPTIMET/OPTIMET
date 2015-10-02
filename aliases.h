#ifndef ALIASES_H_
#define ALIASES_H_

/** @file aliases.h
 * Header file containing common aliases used in the code.
 * Aliases are implemented as preprocessor defines and start with O3D.
 * Alias groups:
 *  1 - solver methods
 *  8 - geometry types
 * 	9 - output grids
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
