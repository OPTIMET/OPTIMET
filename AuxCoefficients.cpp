#include "AuxCoefficients.h"

#include "constants.h"
#include "Tools.h"
#include "Legendre.h"
#include "CompoundIterator.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

#include <boost/math/special_functions/factorials.hpp>

namespace optimet {

std::vector<t_real> AuxCoefficients::compute_dn(t_uint nMax) {
  std::vector<t_real> dn(nMax + 1);

  dn[0] = -1000; // this should return infinity!
  for (t_uint i = 1; i <= nMax; ++i)
    dn[i] = std::sqrt((2.0 * i + 1.0) / (4.0 * consPi * (i * (i + 1))));

  return dn;
}

std::vector<SphericalP<t_complex>>
AuxCoefficients::compute_Pn(t_uint nMax, const std::vector<t_real> &Wigner) {
  std::vector<SphericalP<t_complex>> Pn(nMax + 1);

  for (t_uint i = 0; i <= nMax; ++i) {
    Pn[i].rrr = Wigner[i];
    Pn[i].the = 0.;
    Pn[i].phi = 0.;
  }

  return Pn;
}

std::vector<SphericalP<t_complex>>
AuxCoefficients::compute_Cn(t_uint nMax, t_int m, const Spherical<t_real> &R,
                            const std::vector<t_real> &Wigner,
                            const std::vector<t_real> &dWigner) {
  std::vector<SphericalP<t_complex>> Cn(nMax + 1);

  for (t_uint i = 0; i <= nMax; ++i) {
    t_real A;
    if (m == 0)
      A = 0.0;
    else if (std::abs(R.the) < 1e-10 ||
             (std::abs(R.the) - consPi + 1e-10) > 0.0)
      A = m / std::cos(R.the) * dWigner[i];
    else
      A = m / std::sin(R.the) * Wigner[i];

    Cn[i].rrr = t_complex(0.0, 0.0);
    Cn[i].the = t_complex(0.0, A);
    Cn[i].phi = t_complex(-dWigner[i], 0.0);
  }

  return Cn;
}

std::vector<SphericalP<t_complex>>
AuxCoefficients::compute_Bn(t_uint nMax, t_int m, const Spherical<t_real> &R,
                            const std::vector<t_real> &Wigner,
                            const std::vector<t_real> &dWigner) {
  std::vector<SphericalP<t_complex>> Bn(nMax + 1);

  for (t_uint i = 0; i <= nMax; ++i) {
    t_real A;
    if (m == 0)
      A = 0.0;
    else if (std::abs(R.the) < 1e-10 ||
             (std::abs(R.the) - consPi + 1e-10) > 0.0)
      A = m / std::cos(R.the) * dWigner[i];
    else
      A = m / std::sin(R.the) * Wigner[i];

    Bn[i].rrr = t_complex(0.0, 0.0);
    Bn[i].the = t_complex(dWigner[i], 0.0);
    Bn[i].phi = t_complex(0.0, A);
  }

  return Bn;
}

std::vector<SphericalP<t_complex>>
AuxCoefficients::compute_Mn(t_uint nMax, t_int m, const Spherical<t_real> &R,
                            t_complex waveK, const std::vector<t_real> &dn,
                            const std::vector<SphericalP<t_complex>> &Cn,
                            BESSEL_TYPE besselType) {
  std::vector<SphericalP<t_complex>> Mn(nMax + 1);

  std::vector<t_complex> data, ddata;
  std::tie(data, ddata) = optimet::bessel(R.rrr * waveK, besselType, 0, nMax);

  const t_real dm = std::pow(-1.0, m); // Legendre to Wigner function
  const t_complex exp_imphi(std::cos(m * R.phi), std::sin(m * R.phi));

  for (t_uint i = 0; i <= nMax; ++i) {
    const t_complex c_temp = dm * dn[i] * exp_imphi;

    Mn[i].rrr = c_temp * data[i] * Cn[i].rrr;
    Mn[i].the = c_temp * data[i] * Cn[i].the;
    Mn[i].phi = c_temp * data[i] * Cn[i].phi;
  }

  return Mn;
}

std::vector<SphericalP<t_complex>> AuxCoefficients::compute_Nn(
    t_uint nMax, t_int m, const Spherical<t_real> &R, t_complex waveK,
    const std::vector<t_real> &dn, const std::vector<SphericalP<t_complex>> &Pn,
    const std::vector<SphericalP<t_complex>> &Bn, BESSEL_TYPE besselType) {
  std::vector<SphericalP<t_complex>> Nn(nMax + 1);

  const t_complex Kr = waveK * R.rrr;

  std::vector<t_complex> data, ddata;
  std::tie(data, ddata) = optimet::bessel(R.rrr * waveK, besselType, 0, nMax);

  const t_real dm = std::pow(-1.0, m); // Legendre to Wigner function
  const t_complex exp_imphi(std::cos(m * R.phi), std::sin(m * R.phi));

  Nn[0].rrr = t_complex(0.0, 0.0);
  Nn[0].the = t_complex(0.0, 0.0);
  Nn[0].phi = t_complex(0.0, 0.0);

  for (t_uint i = 1; i <= nMax; ++i) {
    Nn[i].rrr = (1.0 / Kr) * dm * dn[i] *
                ((static_cast<t_real>(i * (i + 1)) * data[i] * Pn[i].rrr) +
                 ((Kr * ddata[i] + data[i]) * Bn[i].rrr)) *
                exp_imphi;
    Nn[i].the = (1.0 / Kr) * dm * dn[i] *
                ((static_cast<t_real>(i * (i + 1)) * data[i] * Pn[i].the) +
                 ((Kr * ddata[i] + data[i]) * Bn[i].the)) *
                exp_imphi;
    Nn[i].phi = (1.0 / Kr) * dm * dn[i] *
                ((static_cast<t_real>(i * (i + 1)) * data[i] * Pn[i].phi) +
                 ((Kr * ddata[i] + data[i]) * Bn[i].phi)) *
                exp_imphi;
  }

  return Nn;
}

std::tuple<std::vector<t_real>, std::vector<t_real>>
AuxCoefficients::VIGdVIG(t_uint nMax, t_int m, const Spherical<t_real> &R) {

  std::vector<t_real> Wigner(nMax + 1, 0.0), dWigner(nMax + 1, 0.0);

  // I and II - determine n_min and obtain vig_d_n_min
  const bool check_m_negative = (m < 0);
  // vector spherical function index : n_min=max(|l|,|m|), where l == 0
  const t_uint n_min = static_cast<t_uint>(m = std::abs(m));
  // wigner function argument : [0<=the<=PI]
  const t_real vig_the = (check_m_negative) ? consPi - R.the : R.the;
  // wigner function auxiliary variable   : x=cos(the)
  const t_real vig_x = std::cos(vig_the); // calculate in radians

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
    for (t_uint i = 2; i <= nMax + 1; ++i) {
      // evaluate from above
      // vig_d        : VIG_d[i] = ...*VIG_d[i-1] + ...*VIG_d[i-2] - eqn(B.22)
      const t_real d_temp =
          ((2 * i - 1) * vig_x * Wigner[i - 1] - (i - 1) * Wigner[i - 2]) / i;

      if (i <= nMax)
        Wigner[i] = d_temp;
      // d_vig_d      : d_VIG_d[i-1] = ...*VIG_d[i-2] + ...*VIG_d[i-1] +
      // ...*VIG_d[i]         - eqn(B.26)
      dWigner[i - 1] =
          ((i * i - i) * d_temp / (2 * i - 1) -
           i * (i * i - 2 * i + 1) * Wigner[i - 2] / (2 * i * i - 3 * i + 1)) /
          std::sin(vig_the);
    }
  }

  // III.2 - else if (m!=0), no special case
  else {
    // obtain the readily availble VIG_d[n_min] value
    using boost::math::factorial;
    Wigner[n_min] =
        std::pow(2.0, -m) *
        (std::sqrt(factorial<t_real>(2 * static_cast<unsigned int>(m))) /
         factorial<t_real>(static_cast<unsigned int>(m))) *
        std::pow(1.0 - vig_x, m / 2.0) *
        std::pow(1.0 + vig_x, m / 2.0); // by definition

    // d_vig_d for the previous from current : d_VIG_d[i] = 0.*VIG_d[i-1] +
    // 0.*VIG_d[i] + ...*VIG_d[i+1]    - eqn(B.26)
    dWigner[n_min - 1] = 0.0; // == zero by definition

    // obtain all other values in VIG_d[i] from recursive relationship
    //   - Based on 'n_min' value and 'n' value; total recursive steps == n -
    //   n_min
    for (t_uint i = n_min + 1; i <= nMax + 1; ++i) {
      // evaluate from terms above
      const t_real d_temp = ((2 * i - 1) * vig_x * Wigner[i - 1] -
                             std::sqrt((i - static_cast<t_uint>(m) - 1) *
                                       (i + static_cast<t_uint>(m) - 1)) *
                                 Wigner[i - 2]) /
                            std::sqrt(i * i - static_cast<t_uint>(m * m));
      // populate VIG_d[i] array
      if (i <= nMax)
        Wigner[i] = d_temp;
      // d_vig_d for the previous from current : d_VIG_d[i-1] = ...*VIG_d[i-2] +
      // ...*VIG_d[i-1] + ...*VIG_d[i]        - eqn(B.26)
      dWigner[i - 1] =
          ((i - 1) * std::sqrt(i * i - static_cast<t_uint>(m * m)) * d_temp /
               (2 * i - 1) -
           i * std::sqrt((i - static_cast<t_uint>(m) - 1) *
                         (i + static_cast<t_uint>(m) - 1)) *
               Wigner[i - 2] / (2 * i - 1)) /
          std::sin(vig_the);
    }
  }

  // IV - if (m<0) : apply symmetry property eq(B.7) to eqs(B.22-B.24)
  if (check_m_negative) {
    for (t_uint i = 0; i <= nMax; ++i) {
      const t_real c = 1.0 / std::pow(-1.0, i);
      Wigner[i] *= c;   // obtain final VIG_d
      dWigner[i] *= -c; // obtain final VIG_d
    }
  }

  return std::make_tuple(Wigner, dWigner);
}

AuxCoefficients::AuxCoefficients(const Spherical<t_real> &R, t_complex waveK,
                                 bool regular, t_uint nMax)
    : _M(Tools::iteratorMax(nMax)), _N(Tools::iteratorMax(nMax)),
      _B(Tools::iteratorMax(nMax)), _C(Tools::iteratorMax(nMax)),
      _dn(compute_dn(nMax)) {

  const BESSEL_TYPE besselType = (regular) ? Bessel : Hankel1;

  for (CompoundIterator q(nMax, nMax); q < q.max(nMax); ++q) {

    // Wigner d function test
    std::vector<t_real> Wigner, dWigner;
    std::tie(Wigner, dWigner) = VIGdVIG(nMax, q.second, R);

    // call vector spherical functions

    const std::vector<SphericalP<t_complex>> Pn = compute_Pn(nMax, Wigner);
    const std::vector<SphericalP<t_complex>> Cn =
        compute_Cn(nMax, q.second, R, Wigner, dWigner);
    const std::vector<SphericalP<t_complex>> Bn =
        compute_Bn(nMax, q.second, R, Wigner, dWigner);

    // call vector spherical waves
    const std::vector<SphericalP<t_complex>> Mn =
        compute_Mn(nMax, q.second, R, waveK, _dn, Cn, besselType);
    const std::vector<SphericalP<t_complex>> Nn =
        compute_Nn(nMax, q.second, R, waveK, _dn, Pn, Bn, besselType);

    for (t_uint n = static_cast<t_uint>(std::abs(q.second)); n <= nMax; ++n) {
      if (n != 0) {
        CompoundIterator p(n, q.second);

        _B[static_cast<t_uint>(p)] = Tools::toProjection(R, Bn[n]);
        _C[static_cast<t_uint>(p)] = Tools::toProjection(R, Cn[n]);
        _M[static_cast<t_uint>(p)] = Tools::toProjection(R, Mn[n]);
        _N[static_cast<t_uint>(p)] = Tools::toProjection(R, Nn[n]);
      }
    }
  }
}

} // namespace optimet
