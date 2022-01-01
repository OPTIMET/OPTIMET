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

#include "Coupling.h"

#include "CompoundIterator.h"
#include "constants.h"
#include "Types.h"
#include "TranslationAdditionCoefficients.h"

#include <cmath>
#include <tuple>

namespace optimet {

namespace {
t_complex coefficients_A(t_int n, t_int m, t_int l, t_int k, TranslationAdditionCoefficients &ta) {
  if(std::abs(k) > l)
    return 0e0;
  auto const factor = 0.5 / std::sqrt(l * (l + 1) * n * (n + 1));
  t_real const c0 = 2 * k * m;
  auto const c1 = std::sqrt((n - m) * (n + m + 1) * (l - k) * (l + k + 1));
  auto const c2 = std::sqrt((n + m) * (n - m + 1) * (l + k) * (l - k + 1));
  return factor * (c0 * ta(n, m, l, k) + c1 * ta(n, m + 1, l, k + 1) + c2 * ta(n, m - 1, l, k - 1));
}

t_complex coefficients_B(t_int n, t_int m, t_int l, t_int k, TranslationAdditionCoefficients &ta) {
  if(std::abs(k) > l)
    return 0e0;
  t_real const a0 = 2 * l + 1;
  t_real const a1 = (2 * l - 1) * l * (l + 1) * n * (n + 1);
  t_complex const factor(0, -0.5 * std::sqrt(a0 / a1));
  auto const c0 = static_cast<t_real>(2 * m) * std::sqrt((l - k) * (l + k));
  auto const c1 = std::sqrt((n - m) * (n + m + 1) * (l - k) * (l - k - 1));
  auto const c2 = std::sqrt((n + m) * (n - m + 1) * (l + k) * (l + k - 1));
  return factor * (c0 * ta(n, m, l - 1, k) + c1 * ta(n, m + 1, l - 1, k + 1) -
                   c2 * ta(n, m - 1, l - 1, k - 1));
}

std::tuple<Matrix<t_complex>, Matrix<t_complex>>
transfer_coefficients(Spherical<double> R, std::complex<double> waveK, bool regular, int n_max) {
  auto const N = Tools::iteratorMax(n_max);
  Matrix<t_complex> diagonal = Matrix<t_complex>::Zero(N, N);
  Matrix<t_complex> offdiagonal = Matrix<t_complex>::Zero(N, N);

  TranslationAdditionCoefficients ta(R, waveK, regular);

  // start at harmonic n = 1. (because n=0 spherical and hence symmetrically incompatible with
  // propagating wave)
  for(t_int n(1); n <= n_max; ++n)
    for(t_int m(-n); m <= n; ++m) {

      auto const p = flatten_indices(n, m);
      for(t_int l(1); l <= n_max; ++l)
        for(t_int k(-l); k <= l; ++k) {

          auto const q = flatten_indices(l, k);
          diagonal(p, q) = coefficients_A(n, m, l, k, ta);
          offdiagonal(p, q) = coefficients_B(n, m, l, k, ta);
        }
    }
  return std::make_tuple(diagonal, offdiagonal);
}

} // anonymous namespace

Coupling::Coupling(Spherical<t_real> relR, t_complex waveK, t_uint nMax, bool regular) {
  auto const n = Tools::iteratorMax(nMax);
  if(std::abs(relR.rrr) < errEpsilon) { // Check for NO translation case
    offdiagonal = Matrix<t_complex>::Zero(n, n);
    diagonal = Matrix<t_complex>::Identity(n, n);
  } else
    std::tie(diagonal, offdiagonal) = transfer_coefficients(relR, waveK, not regular, nMax);
}
} // namespace optimet
