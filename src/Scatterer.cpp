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

#include "Bessel.h"
#include "HarmonicsIterator.h"
#include "Scatterer.h"
#include "Tools.h"

Scatterer::Scatterer(Spherical<double> vR_, ElectroMagnetic elmag_, double radius_, int nMax_)
    : vR(vR_), elmag(elmag_), radius(radius_), nMax(nMax_),
      sourceCoef(2 * Tools::iteratorMax(nMax)) {}

Scatterer::Scatterer(int nMax) : Scatterer({0, 0, 0}, {0, 0}, 0e0, nMax) {}

Scatterer::~Scatterer() {}

optimet::Vector<optimet::t_complex>
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
  Vector<t_complex> result = Vector<t_complex>::Zero(2 * N);
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
    result.segment(current, 2 * n + 1).fill(TE);

    // TM part
    auto const TM = (psi / ksi) * (mu_sob * dpsirho / psirho - rho * dpsi / psi) /
                    (rho * dksi / ksi - mu_sob * dpsirho / psirho);
    result.segment(current + N, 2 * n + 1).fill(TM);
  }

  return result;
}

optimet::Vector<optimet::t_complex>
Scatterer::getIaux(optimet::t_real omega_, ElectroMagnetic const &bground) const {
  auto const k_s = omega_ * std::sqrt(elmag.epsilon * elmag.mu);
  auto const k_b = omega_ * std::sqrt(bground.epsilon * bground.mu);
  auto const rho = k_s / k_b;
  auto const r_0 = k_b * radius;
  auto const mu_j = elmag.mu;
  auto const mu_0 = bground.mu;

  auto Jdata = optimet::bessel<optimet::Bessel>(r_0, nMax);
  auto const Jrho = optimet::bessel<optimet::Bessel>(rho * r_0, nMax);

  optimet::Vector<optimet::t_complex> result(2 * nMax * (nMax + 2));
  auto TE = result.head(nMax * (nMax + 2));
  auto TM = result.tail(nMax * (nMax + 2));
  for(auto n = 1, i = 0; n <= nMax; ++n) {
    // obtain Riccati-Bessel functions
    auto const psi = r_0 * std::get<0>(Jdata)[n];
    auto const dpsi = r_0 * std::get<1>(Jdata)[n] + std::get<0>(Jdata)[n];
    auto const psirho = r_0 * rho * std::get<0>(Jrho)[n];
    auto const dpsirho = r_0 * rho * std::get<1>(Jrho)[n] + std::get<0>(Jrho)[n];

    for(auto m = -n; m <= n; ++m, ++i) {
      TE(i) = (mu_j * rho) / (mu_0 * rho * dpsirho * psi - mu_j * psirho * dpsi) *
              std::complex<double>(0., 1.);
      TM(i) = (mu_j * rho) / (mu_j * psi * dpsirho - mu_0 * rho * psirho * dpsi) *
              std::complex<double>(0., 1.);
    }
  }
  return result;
}
