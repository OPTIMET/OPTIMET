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

#ifndef TRANSLATION_ADDITION_COEFFICIENTS_H

#include "Types.h"
#include <array>
#include <map>

#include "Spherical.h"

namespace optimet {
//! Equation 1 of Appendix A in Stout (2002)
t_complex Ynm(Spherical<t_real> const &R, t_int n, t_int m);

namespace details {

//! \brief Computes translation addition coefficients for m > 0
//! \details Equations come from Stout (2002), appendix C. Negative m coefficients should be
//! obtained using the symmetry relationship. This class caches previously computed coefficients.
class CachedRecurrence {
public:
  //! Indices tuple
  typedef std::array<t_int, 4> t_indices;

  CachedRecurrence(Spherical<t_real> R, t_complex waveK, bool regular = true)
      : direction(R), waveK(waveK), regular(regular) {}

  //! \brief Returns translation addition coefficients
  //! \details n and m correspond to the same variables in Stout (2004), l and k correspond to ν and
  //! μ, respectively.
  t_complex operator()(t_int n, t_int m, t_int l, t_int k);

  bool is_regular() const { return regular; }

protected:
  //! Direction of the incident wave
  Spherical<t_real> const direction;
  //! Wavenumber of the incident wave
  t_complex const waveK;
  //! Whether this is for regular or irregular coeffs
  bool const regular;
  //! Caches known coefficients
  std::map<t_indices, t_complex> cache;

  //! Switches between recurrence relationships
  t_complex recurrence(t_int n, t_int m, t_int l, t_int k);
  t_complex initial(t_int l, t_int k);
  t_complex diagonal_recurrence(t_int n, t_int l, t_int k);
  t_complex offdiagonal_recurrence(t_int n, t_int m, t_int l, t_int k);
};

} // end of details namespace

//! \brief Computes and caches translation-addition coefficients
//! \details The coefficients are obtained for given wave. It is not possible to change it once set.
class TranslationAdditionCoefficients {
public:
  TranslationAdditionCoefficients(Spherical<t_real> R, t_complex waveK, bool regular = true)
      : positive(R, waveK, regular),
        negative(R, regular ? std::conj(waveK) : -std::conj(waveK), regular) {}

  //! \brief Computes the coefficients as per Stout (2002)
  //! \details n, m, l, k correspond to n, m, ν, μ in Stout (2002), respectively.
  t_complex operator()(t_int n, t_int m, t_int l, t_int k);

protected:
  //! Recurrence for positive m
  details::CachedRecurrence positive;
  //! Recurrence for negative m
  details::CachedRecurrence negative;
};
}

#endif
