#include "AuxCoefficients.h"

#include "constants.h"
#include "Tools.h"
#include "Legendre.h"
#include "CompoundIterator.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "gsl/gsl_sf_gamma.h"

namespace optimet {

std::vector<double> AuxCoefficients::compute_dn(int nMax) {
  std::vector<double> dn(nMax + 1);

  dn[0] = -1000; // this should return infinity!
  for (int i = 1; i <= nMax; i++)
    dn[i] = std::sqrt((2.0 * (double)i + 1.0) /
                      (4.0 * consPi * (double)(i * (i + 1))));

  return dn;
}

std::vector<SphericalP<std::complex<double>>>
AuxCoefficients::compute_Pn(int nMax, const std::vector<double> &Wigner) {
  std::vector<SphericalP<std::complex<double>>> Pn(nMax + 1);

  for (int i = 0; i <= nMax; i++) {
    Pn[i].rrr = Wigner[i];
    Pn[i].the = 0.;
    Pn[i].phi = 0.;
  }

  return Pn;
}

std::vector<SphericalP<std::complex<double>>>
AuxCoefficients::compute_Cn(int nMax, int m, const Spherical<double> &R,
                            const std::vector<double> &Wigner,
                            const std::vector<double> &dWigner) {
  std::vector<SphericalP<std::complex<double>>> Cn(nMax + 1);

  for (int i = 0; i <= nMax; i++) {
    double A;
    if (m == 0)
      A = 0.0;
    else if (std::abs(R.the) < 1e-10 ||
             (std::abs(R.the) - consPi + 1e-10) > 0.0)
      A = double(m) / std::cos(R.the) * dWigner[i];
    else
      A = double(m) / std::sin(R.the) * Wigner[i];

    Cn[i].rrr = std::complex<double>(0.0, 0.0);
    Cn[i].the = std::complex<double>(0.0, A);
    Cn[i].phi = std::complex<double>(-dWigner[i], 0.0);
  }

  return Cn;
}

std::vector<SphericalP<std::complex<double>>>
AuxCoefficients::compute_Bn(int nMax, int m, const Spherical<double> &R,
                            const std::vector<double> &Wigner,
                            const std::vector<double> &dWigner) {
  std::vector<SphericalP<std::complex<double>>> Bn(nMax + 1);

  for (int i = 0; i <= nMax; i++) {
    double A;
    if (m == 0)
      A = 0.0;
    else if (std::abs(R.the) < 1e-10 ||
             (std::abs(R.the) - consPi + 1e-10) > 0.0)
      A = double(m) / std::cos(R.the) * dWigner[i];
    else
      A = double(m) / std::sin(R.the) * Wigner[i];

    Bn[i].rrr = std::complex<double>(0.0, 0.0);
    Bn[i].the = std::complex<double>(dWigner[i], 0.0);
    Bn[i].phi = std::complex<double>(0.0, A);
  }

  return Bn;
}

std::vector<SphericalP<std::complex<double>>> AuxCoefficients::compute_Mn(
    int nMax, int m, const Spherical<double> &R, std::complex<double> waveK,
    const std::vector<double> &dn,
    const std::vector<SphericalP<std::complex<double>>> &Cn,
    BESSEL_TYPE besselType) {
  std::vector<SphericalP<std::complex<double>>> Mn(nMax + 1);

  std::vector<std::complex<double>> data, ddata;
  std::tie(data, ddata) = optimet::bessel(R.rrr * waveK, besselType, 0, nMax);

  const double dm = std::pow(-1.0, double(m)); // Legendre to Wigner function
  const std::complex<double> exp_imphi(std::cos(double(m) * R.phi),
                                       std::sin(double(m) * R.phi));

  for (int i = 0; i <= nMax; i++) {
    const std::complex<double> c_temp = dm * dn[i] * exp_imphi;

    Mn[i].rrr = c_temp * data[i] * Cn[i].rrr;
    Mn[i].the = c_temp * data[i] * Cn[i].the;
    Mn[i].phi = c_temp * data[i] * Cn[i].phi;
  }

  return Mn;
}

std::vector<SphericalP<std::complex<double>>> AuxCoefficients::compute_Nn(
    int nMax, int m, const Spherical<double> &R, std::complex<double> waveK,
    const std::vector<double> &dn,
    const std::vector<SphericalP<std::complex<double>>> &Pn,
    const std::vector<SphericalP<std::complex<double>>> &Bn,
    BESSEL_TYPE besselType) {
  std::vector<SphericalP<std::complex<double>>> Nn(nMax + 1);

  const std::complex<double> Kr = waveK * R.rrr;

  std::vector<std::complex<double>> data, ddata;
  std::tie(data, ddata) = optimet::bessel(R.rrr * waveK, besselType, 0, nMax);

  const double dm = std::pow(-1.0, double(m)); // Legendre to Wigner function
  const std::complex<double> exp_imphi(std::cos(double(m) * R.phi),
                                       std::sin(double(m) * R.phi));

  Nn[0].rrr = std::complex<double>(0.0, 0.0);
  Nn[0].the = std::complex<double>(0.0, 0.0);
  Nn[0].phi = std::complex<double>(0.0, 0.0);

  for (int i = 1; i <= nMax; i++) {
    const double d_n = double(i);
    Nn[i].rrr = (1.0 / Kr) * dm * dn[i] *
                ((d_n * (d_n + 1.0) * data[i] * Pn[i].rrr) +
                 ((Kr * ddata[i] + data[i]) * Bn[i].rrr)) *
                exp_imphi;
    Nn[i].the = (1.0 / Kr) * dm * dn[i] *
                ((d_n * (d_n + 1.0) * data[i] * Pn[i].the) +
                 ((Kr * ddata[i] + data[i]) * Bn[i].the)) *
                exp_imphi;
    Nn[i].phi = (1.0 / Kr) * dm * dn[i] *
                ((d_n * (d_n + 1.0) * data[i] * Pn[i].phi) +
                 ((Kr * ddata[i] + data[i]) * Bn[i].phi)) *
                exp_imphi;
  }

  return Nn;
}

std::tuple<std::vector<double>, std::vector<double>>
AuxCoefficients::VIGdVIG(int nMax, int m, const Spherical<double> &R) {

  std::vector<double> Wigner(nMax + 1, 0.0), dWigner(nMax + 1, 0.0);

  // I and II - determine n_min and obtain vig_d_n_min
  const bool check_m_negative = (m < 0);
  // vector spherical function index : n_min=max(|l|,|m|), where l == 0
  const int n_min = m = std::abs(m);
  // wigner function argument : [0<=the<=PI]
  const double vig_the = (check_m_negative) ? consPi - R.the : R.the;
  // wigner function auxiliary variable   : x=cos(the)
  const double vig_x = std::cos(vig_the); // calculate in radians

  // III - calculate VIG_d & d_VIG_d
  // III.1 - check for singularity cases - i.e (m==0)
  if (m == 0) {
    // obtain the readily availble VIG_d[n_min] value
    Wigner[0] = 1.0; // by definition
    // d_vig_d for the previous from current : d_VIG_d[i] = 0.*VIG_d[i-1] +
    // 0.*VIG_d[i] + ...*VIG_d[i+1]    - eqn(B.26)
    dWigner[0] = 0.0; // == zero by definition
    // obtain first term in recursive relationship eq(B.22) - from special case
    Wigner[1] = vig_x;
    // obtain all other values in VIG_d[i] from recursive relationship
    for (int i = 2; i <= nMax + 1; i++) {
      // evaluate from above
      // vig_d        : VIG_d[i] = ...*VIG_d[i-1] + ...*VIG_d[i-2] - eqn(B.22)
      const double d_temp = 1.0 / ((double)(i * i * (i - 1))) *
                                (double)((2 * i - 1) * i * (i - 1)) * vig_x *
                                Wigner[i - 1] -
                            (double)(i * (i - 1) * (i - 1)) * Wigner[i - 2];
      if (i <= nMax)
        Wigner[i] = d_temp;
      // d_vig_d      : d_VIG_d[i-1] = ...*VIG_d[i-2] + ...*VIG_d[i-1] +
      // ...*VIG_d[i]         - eqn(B.26)
      dWigner[i - 1] =
          ((-(double)(i * (i - 1) * (i - 1)) / (double)((i - 1)) *
            (2 * i - 1)) *
               Wigner[i - 2] +
           ((double)(i * i * (i - 1)) / (double)(2 * i * i - i)) * d_temp) /
          std::sin(vig_the);
    }
  }

  // III.2 - else if (m!=0), no special case
  else {
    // obtain the readily availble VIG_d[n_min] value
    const double vig_d_n_min =
        std::pow(2.0, -m) *
        (std::sqrt((double)gsl_sf_fact(2 * m)) / (double)gsl_sf_fact(m)) *
        std::pow(1.0 - vig_x, (double)m / 2.0) *
        std::pow(1.0 + vig_x, (double)m / 2.0);
    Wigner[m] = vig_d_n_min; // by definition

    // d_vig_d for the previous from current : d_VIG_d[i] = 0.*VIG_d[i-1] +
    // 0.*VIG_d[i] + ...*VIG_d[i+1]    - eqn(B.26)
    dWigner[m - 1] = 0.0; // == zero by definition

    // obtain all other values in VIG_d[i] from recursive relationship
    //   - Based on 'n_min' value and 'n' value; total recursive steps == n -
    //   n_min
    for (int i = m + 1; i <= nMax + 1; i++) {
      // evaluate from terms above
      const double d_temp =
          -((1.0 - (double)i) / 2.0) / std::sqrt((double)(i * i - m * m)) *
          ((double)((2 * i - 1) * (i - 1) * i) * vig_x * Wigner[i - 1] -
           (double)(i * (i - 1)) *
               std::sqrt((double)((i - m - 1) * (i + m - 1))) * Wigner[i - 2]);
      // populate VIG_d[i] array
      if (i <= nMax)
        Wigner[i] = d_temp;
      // d_vig_d for the previous from current : d_VIG_d[i-1] = ...*VIG_d[i-2] +
      // ...*VIG_d[i-1] + ...*VIG_d[i]        - eqn(B.26)
      dWigner[i - 1] =
          (-(double)(i * (i - 1)) *
               std::sqrt((double)((i - m - 1) * (i + m - 1))) /
               (double)((i - 1) * (2 * i - 1)) * Wigner[i - 2] +
           (double)(i * (i - 1)) * std::sqrt((double)((i * i) - (m * m))) /
               (double)(i * (2 * i - 1)) * d_temp) /
          std::sin(vig_the);
    }
  }

  // IV - if (m<0) : apply symmetry property eq(B.7) to eqs(B.22-B.24)
  if (check_m_negative) {
    for (int i = 0; i <= nMax; i++) {
      const double c = 1.0 / std::pow(-1.0, i);
      Wigner[i] *= c;   // obtain final VIG_d
      dWigner[i] *= -c; // obtain final VIG_d
    }
  }

  return std::make_tuple(Wigner, dWigner);
}

AuxCoefficients::AuxCoefficients(const Spherical<double> &R,
                                 std::complex<double> waveK, bool regular,
                                 int nMax)
    : _M(Tools::iteratorMax(nMax)), _N(Tools::iteratorMax(nMax)),
      _B(Tools::iteratorMax(nMax)), _C(Tools::iteratorMax(nMax)),
      _dn(compute_dn(nMax)) {

  const BESSEL_TYPE besselType = (regular) ? Bessel : Hankel1;

  for (CompoundIterator q(nMax, nMax); q < q.max(nMax); ++q) {

    // Wigner d function test
    std::vector<double> Wigner, dWigner;
    std::tie(Wigner, dWigner) = VIGdVIG(nMax, q.second, R);

    // call vector spherical functions

    const std::vector<SphericalP<std::complex<double>>> Pn =
        compute_Pn(nMax, Wigner);
    const std::vector<SphericalP<std::complex<double>>> Cn =
        compute_Cn(nMax, q.second, R, Wigner, dWigner);
    const std::vector<SphericalP<std::complex<double>>> Bn =
        compute_Bn(nMax, q.second, R, Wigner, dWigner);

    // call vector spherical waves
    const std::vector<SphericalP<std::complex<double>>> Mn =
        compute_Mn(nMax, q.second, R, waveK, _dn, Cn, besselType);
    const std::vector<SphericalP<std::complex<double>>> Nn =
        compute_Nn(nMax, q.second, R, waveK, _dn, Pn, Bn, besselType);

    for (int n = std::abs(q.second); n <= nMax; n++) {
      if (n != 0) {
        CompoundIterator p(n, q.second);

        _B[static_cast<long>(p)] = Tools::toProjection(R, Bn[n]);
        _C[static_cast<long>(p)] = Tools::toProjection(R, Cn[n]);
        _M[static_cast<long>(p)] = Tools::toProjection(R, Mn[n]);
        _N[static_cast<long>(p)] = Tools::toProjection(R, Nn[n]);
      }
    }
  }
}

} // namespace optimet
