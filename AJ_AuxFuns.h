#ifndef AJ_AuxFuns_H_
#define AJ_AuxFuns_H_

#include "Spherical.h"
#include "CompoundIterator.h"
#include "Spherical.h"
#include "Bessel.h"
#include "Tools.h"
#include "constants.h"
#include "gsl/gsl_sf_gamma.h"
#include <math.h>
#include <iostream>

#include "Scatterer.h"
#include "Tools.h"
#include "Bessel.h"
#include "CompoundIterator.h"
#include "Tools.h"
#include "Symbol.h"
#include "Excitation.h"
#include "Algebra.h"
#include <iostream>
#include <cmath>

int compute_YJn_m(Spherical<double> R, complex<double> waveK, int BHreg, int nMax, int m, complex<double> *YJnm);
int compute_YJp(Spherical<double> R, complex<double> waveK, int BHreg, int nMax, complex<double> *dataYJp);

int compute_Yn_m(Spherical<double> R, complex<double> waveK, int nMax, int m, complex<double> *Ynm);
int compute_Yp(Spherical<double> R, complex<double> waveK, int nMax, complex<double> *dataYp);

//int getTLocal(double omega_, Scatterer particle, int nMax_, complex<double> ** T_local_);

#endif /* AJ_AuxFuns_H_ */
