#include "Scatterer.h"
#include "Bessel.h"
#include "HarmonicsIterator.h"
#include "Tools.h"

Scatterer::Scatterer(Spherical<double> vR_, ElectroMagnetic elmag_, double radius_, int nMax_)
    : vR(vR_), elmag(elmag_), radius(radius_), nMax(nMax_),
      sourceCoef(2 * Tools::iteratorMax(nMax)) {}

Scatterer::Scatterer(int nMax) : Scatterer({0, 0, 0}, {0, 0}, 0e0, nMax) {}

Scatterer::~Scatterer() {}

optimet::Matrix<optimet::t_complex>
Scatterer::getTLocal(optimet::t_real omega_, ElectroMagnetic const &bground) const {
  using namespace optimet;
  auto const k_s = omega_ * std::sqrt(elmag.epsilon * elmag.mu);
  auto const k_b = omega_ * std::sqrt(bground.epsilon * bground.mu);

  auto const rho = k_s / k_b;
  auto const r_0 = k_b * radius;
  auto const mu_sob = elmag.mu / bground.mu;

  auto const Jn = bessel<Bessel>(r_0, nMax);
  auto const Jrho = bessel<Bessel>(rho * r_0, nMax);
  auto const Hn = bessel<Hankel1>(r_0, nMax);

  auto const N = HarmonicsIterator::max_flat(nMax) - 1;
  Matrix<t_complex> result = Matrix<t_complex>::Zero(2 * N, 2 * N);
  for(t_uint n(1), current(0); n <= nMax; current += 2 * n + 1, ++n) {
    auto const psi = r_0 * std::get<0>(Jn)[n];
    auto const dpsi = r_0 * std::get<1>(Jn)[n] + std::get<0>(Jn)[n];

    auto const ksi = r_0 * std::get<0>(Hn)[n];
    auto const dksi = r_0 * std::get<1>(Hn)[n] + std::get<0>(Hn)[n];

    auto const psirho = r_0 * rho * std::get<0>(Jrho)[n];
    auto const dpsirho = r_0 * rho * std::get<1>(Jrho)[n] + std::get<0>(Jrho)[n];

    // TE Part
    auto const TE = (psi / ksi) * (mu_sob * dpsi / psi - rho * dpsirho / psirho) /
                    (rho * dpsirho / psirho - mu_sob * dksi / ksi);
    result.diagonal().segment(current, 2 * n + 1).fill(TE);

    // TM part
    auto const TM = (psi / ksi) * (mu_sob * dpsirho / psirho - rho * dpsi / psi) /
                    (rho * dksi / ksi - mu_sob * dpsirho / psirho);
    result.diagonal().segment(current + N, 2 * n + 1).fill(TM);
  }

  return result;
}
