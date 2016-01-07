#ifndef AJ_AuxFuns_H_
#define AJ_AuxFuns_H_

#include "Spherical.h"

#include <vector>
#include <complex>

namespace optimet {

std::vector<std::complex<double>> compute_Yp(const Spherical<double> &R,
                                             const int nMax);
} // namespace optimet

#endif /* AJ_AuxFuns_H_ */
