#include "constants.h"

extern const double consPi = 3.14159265358979323846;
extern const double consC = 299792458;
extern const double consMu0 = 4.0 * consPi * 1e-7;
extern const double consEpsilon0 = 1.0 / (consMu0 * consC * consC);

extern const double consFrnmTom = 1e-9;

extern const std::complex<double> consC1 = std::complex<double>(1.0, 0.0);
extern const std::complex<double> consCi = std::complex<double>(0.0, 1.0);
extern const std::complex<double> consCmi = std::complex<double>(0.0, -1.0);
extern const std::complex<double> consCm1 = std::complex<double>(-1.0, 0.0);
extern const std::complex<double> consC0 = std::complex<double>(0.0, 0.0);

extern const double errEpsilon = 1e-10;

namespace optimet {
namespace constant {
extern const double pi = consPi;
extern const t_real epsilon0 = consEpsilon0;
extern const t_real mu0 = consMu0;
extern const t_real c = consC;

extern const t_real from_nm_to_m = consFrnmTom;

extern const t_complex i{0, 1};

extern const t_real tolerance = errEpsilon;
}
}
