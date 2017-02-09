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

#include "Types.h"
#include "HarmonicsIterator.h"

namespace optimet {
void HarmonicsIterator::throw_on_invalid_nm(t_uint n, t_int m) {
  if(not valid_nm(n, m))
    throw std::runtime_error("abs(m) > n");
}

HarmonicsIterator& HarmonicsIterator::principal(t_uint n) {
  principal_ = n;
  if(not valid_nm(n, m()))
    secondary_ = -static_cast<t_int>(n);
  return *this;
}

HarmonicsIterator& HarmonicsIterator::secondary(t_int m) {
  secondary_ = m;
  if(not valid_nm(n(), m))
    principal_ = static_cast<t_uint>(std::abs(m));
  return *this;
}

HarmonicsIterator HarmonicsIterator::operator++(int) {
  HarmonicsIterator const c(*this);
  operator++();
  return c;
}
HarmonicsIterator& HarmonicsIterator::operator++() {
  --secondary_;
  if(not valid_nm(principal_, secondary_)) {
    ++principal_;
    secondary_ = static_cast<t_int>(principal_);
  }
  return *this;
}
HarmonicsIterator HarmonicsIterator::operator--(int) {
  HarmonicsIterator const c(*this);
  operator--();
  return c;
}
HarmonicsIterator& HarmonicsIterator::operator--() {
  // Absolute limit when decrementing
  if(principal_ == std::numeric_limits<t_uint>::max())
    return *this;

  ++secondary_;
  if(not valid_nm(principal_, secondary_)) {
    --principal_;
    secondary_ = -static_cast<t_int>(principal_);
  }
  return *this;
}

} /* optimet */
