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

#ifndef OPTIMET_HARMONICS_ITERATOR_H
#define OPTIMET_HARMONICS_ITERATOR_H

#include "Types.h"
#include <iterator>
#include <iostream>

namespace optimet {
//! \brief Iterator over harmonics parameters (n, m)
//! \note In order to be more compatible with the original CompoundIterator, the m are increment
//! backwards: (0, 0) -> (1, 1) -> (1, 0) -> (1, -1) -> (2, 2) ...
class HarmonicsIterator : public std::iterator<std::bidirectional_iterator_tag, t_uint const> {
public:
  using std::iterator<std::bidirectional_iterator_tag, t_uint const>::value_type;
  using std::iterator<std::bidirectional_iterator_tag, t_uint const>::pointer;
  using std::iterator<std::bidirectional_iterator_tag, t_uint const>::reference;
  using std::iterator<std::bidirectional_iterator_tag, t_uint const>::difference_type;
  using std::iterator<std::bidirectional_iterator_tag, t_uint const>::iterator_category;

  //! Sets the principal and secondary number
  HarmonicsIterator(t_uint n, t_int m) : principal_(n), secondary_(m) { throw_on_invalid_nm(n, m); }
  //! Sets the principal number and secondary number to -n
  HarmonicsIterator(t_uint n = 0) : HarmonicsIterator(n, static_cast<t_int>(n)) {}

  //! \details Principal number
  //! \details Rather then let the secondary number become invalid, this function will modify (if
  //! necessary) and set it to the -n;
  HarmonicsIterator &principal(t_uint n);
  //! Alias for the principal number
  HarmonicsIterator &n(t_uint n) { return principal(n); }
  //! Principal number
  t_uint principal() const { return principal_; }
  //! Alias for the principal number
  t_uint n() const { return principal_; }
  //! \brief Secondary number -n ≤ m ≤ n
  //! \details Rather then let the principal number become invalid, this function will modify (if
  //! necessary) and set it to the absolute value of m.
  HarmonicsIterator &secondary(t_int m);
  //! Secondary number -n ≤ m ≤ n
  t_int secondary() const { return secondary_; }
  //! Alias for the secondary number
  t_int m() const { return secondary(); }
  //! Alias for the secondary number
  HarmonicsIterator &m(t_int m) { return secondary(m); }

  //! The flat index corresponding to the current n and m;
  value_type operator*() const { return static_cast<t_int>(n() * (n() + 1)) - m(); }
  //! The flat index corresponding to the current n and m;
  pointer operator->() const {
    flat = operator*();
    return &flat;
  }

  bool operator==(HarmonicsIterator const &c) const {
    return (n() == std::numeric_limits<t_uint>::max() and
            c.n() == std::numeric_limits<t_uint>::max()) or
           (n() == c.n() and m() == c.m());
  }
  bool operator!=(HarmonicsIterator const &c) const { return not operator==(c); }
  difference_type operator-(HarmonicsIterator const &c) const { return operator*() - *c; }
  HarmonicsIterator &operator++();
  HarmonicsIterator &operator--();
  HarmonicsIterator operator--(int);
  HarmonicsIterator operator++(int);

  //! One before given principal number
  static HarmonicsIterator least(t_uint n) {
    return n == 0 ? least() : HarmonicsIterator(n - 1, static_cast<t_int>(n - 1));
  }
  //! One before (0, 0). Usefull when decrementing to zero.
  static HarmonicsIterator least() { return HarmonicsIterator(std::numeric_limits<t_uint>::max()); }
  //! One past the given principal number
  static HarmonicsIterator end(t_uint n) { return HarmonicsIterator(n + 1); }
  //! Maximum flat index
  static t_uint max_flat(t_uint n) { return *HarmonicsIterator::end(n); }

private:
  //! Principal number n
  t_uint principal_;
  //! Secondary number -n ≤ m ≤ n
  t_int secondary_;

  //! Used only to provide something to point to
  mutable t_uint flat;

  static bool valid_nm(t_uint n, t_int m) { return static_cast<t_uint>(std::abs(m)) <= n; }
  static void throw_on_invalid_nm(t_uint n, t_int m);
};

inline std::ostream &operator<<(std::ostream &stream, HarmonicsIterator const &c) {
  return stream << "{" << c.principal() << ", " << c.secondary() << "}";
}
} /* optimet */

#endif
