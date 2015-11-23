#ifndef AJ_AuxFuns_H_
#define AJ_AuxFuns_H_

#include "Spherical.h"

#include <complex>

int compute_YJn_m(Spherical<double> R, std::complex<double> waveK, int BHreg,
  int nMax, int m, std::complex<double> *YJnm);
int compute_YJp(Spherical<double> R, std::complex<double> waveK, int BHreg,
  int nMax, std::complex<double> *dataYJp);

int compute_Yn_m(Spherical<double> R, std::complex<double> waveK, int nMax,
  int m, std::complex<double> *Ynm);
int compute_Yp(Spherical<double> R, std::complex<double> waveK, int nMax,
  std::complex<double> *dataYp);

#endif /* AJ_AuxFuns_H_ */
