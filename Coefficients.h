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

#ifndef OPTIMET_COEFFICIENTS_H
#define OPTIMET_COEFFICIENTS_H

#include "Types.h"
#include <cmath>

namespace optimet {
namespace coefficient {
template <class TYPE = t_real> TYPE a(t_uint n, t_int m) {
  t_uint const absm(std::abs(m));
  if(n < absm)
    return static_cast<TYPE>(0);
  return std::sqrt(static_cast<TYPE>((n + 1 + absm) * (n + 1 - absm)) /
                   static_cast<TYPE>((2 * n + 1) * (2 * n + 3)));
}

template <class TYPE = t_real> TYPE b(t_uint n, t_int m) {
  if(static_cast<t_uint>(std::abs(m)) > n)
    return static_cast<TYPE>(0);
  return (m >= 0 ? 1 : -1) * std::sqrt(static_cast<TYPE>((n - m - 1) * (n - m)) /
                                       static_cast<TYPE>((2 * n - 1) * (2 * n + 1)));
}

template <class TYPE = t_real> TYPE c(t_uint n, t_int m) {
  if(std::abs(m) > n)
    return 0;
  return (m < 0 ? -1 : 1) *
         std::sqrt((static_cast<t_int>(n) - m) * (static_cast<t_int>(n) + m + 1));
}
}
}
#endif
