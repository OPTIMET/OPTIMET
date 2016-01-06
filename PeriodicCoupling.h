/*
 * Periodic_Couplin.h
 *
 *  Created on: June 27, 2013
 *      Author: uceealj
 */

#ifndef Periodic_Coupling_H_
#define Periodic_Coupling_H_

#include "Types.h"
#include "Spherical.h"
#include "Cartesian.h"
#include <complex>

namespace optimet {
void compute_AlBe_nmlk(Spherical<t_real> R, t_complex waveK, t_int BHreg,
                       t_int n_max, t_complex ****AlBe_nmlk);
}
#endif /*Periodic_Coupling_H_ */
