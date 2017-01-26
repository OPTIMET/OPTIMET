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

#ifndef COMPOUNDITERATOR_H_
#define COMPOUNDITERATOR_H_

#include "Types.h"

namespace optimet {
//! Harmonic index n, m but flat
inline constexpr t_int flatten_indices(t_int first, t_int second) {
  return first * (first + 1) - second - 1;
}

//! Maximum flat harmonic index
inline constexpr t_int max_flat_index(t_int first) {
  return first * (first + 2);
}
}

/**
 * The CompoundIterator class implements compound iterators.
 * Iterators map p -> (n, m) -> p using a specified mapping function.
 */
class CompoundIterator {
private:
  /**
   * Private function for (n,m) -> p.
   */
  void forwardMap();
  /**
   * Private function for p -> (n,m).
   */
  void backwardMap();

public:
  int compound; /**< The compound value of the iterator. */
  int first;    /**< The first member of the iterator. */
  int second;   /**< The second member of the iterator. */

  /**
   * Default constructor for the CompoundIterator class.
   * Creates an iterator of the type 0 -> (1, 1).
   */
  CompoundIterator();

  /**
   * Alternative constructor for the CompoundIterator class.
   * Creates an iterator with compound_ and re-maps p->(n,m).
   * @param compound_ the value of the compound iterator.
   */
  CompoundIterator(int compound_);

  /**
   * Alternative constructor for the CompoundIterator class.
   * Creates and iterator with (first_, second_) and re-maps (n,m) -> p.
   * @param first_ the first value of the iterator.
   * @param second_ the second value of the iterator.
   */
  CompoundIterator(int first_, int second_);

  /**
   * Initialize a CompoundIterator.
   * @param compound_ the compound iterator value.
   */
  void init(int compound_);

  /**
   * Initialize a CompoundIterator.
   * @param first_ the first iterator value.
   * @param second_ the second iterator value.
   */
  void init(int first_, int second_);

  /**
   * Returns the maximum value needed for compound to cover an iterator up to
   * first_.
   * Assumes that second will go from -first_ to first_.
   * @param first_ the maximum value for first.
   * @return the maximum value for compound.
   */
  static int max(int first_);

  /**
   * Default destructor for the CompoundIterator class.
   */
  virtual ~CompoundIterator();

  /**
   * Operator overload for ++.
   * Increases the compound iterator and re-maps p -> (n,m).
   */
  void operator++();

  /**
   * Operator overload for ++.
   * Increases the compound iterator and re-maps p -> (n,m).
   */
  void operator++(int compound_);

  /**
   * Operator overload for --.
   * Decreases the compound iterator and re-maps p -> (n,m).
   */
  void operator--();

  /**
   * Operator overload for --.
   * Decreases the compound iterator and re-maps p -> (n,m).
   */
  void operator--(int compound_);

  /**
   * Operator overload for +.
   * Increases the compound iterator by increase_ and re-maps p -> (n,m).
   * @param increase_ the value to increase the iterator by.
   * @return the new iterator.
   */
  CompoundIterator operator+(int increase_);

  /**
   * Operator overload for +.
   * Decreases the compound iterator by decrease_ and re-maps p -> (n,m).
   * @param increase_ the value to decrease the iterator by.
   * @return the new iterator.
   */
  CompoundIterator operator-(int decrease_);

  /**
   * Operator overload for =.
   * Assigns a new compound iterator value and re-maps p -> (n,m).
   * @param compound_ the value for the compound iterator.
   */
  void operator=(int compound_);

  /**
   * Operator overload for <.
   * Compares current compound iterator with a compound value.
   * @param compound_ the compound value.
   * @return true or false.
   */
  bool operator<(int compound_);

  /**
   * Operator overload for >.
   * Compares current compound iterator with a compound value.
   * @param compound_ the compound value.
   * @return true or false.
   */
  bool operator>(int compound_);

  /**
   * Operator overload for <=.
   * Compares current compound iterator with a compound value.
   * @param compound_ the compound value.
   * @return true or false.
   */
  bool operator<=(int compound_);

  /**
   * Operator overload for >=.
   * Compares current compound iterator with a compound value.
   * @param compound_ the compound value.
   * @return true or false.
   */
  bool operator>=(int compound_);

  /**
   * Operator overload for ==.
   * Compares current compound iterator with a compound value.
   * @param compound_ the compound value.
   * @return true or false.
   */
  bool operator==(int compound_);

  /**
   * Operator overload for int typecast.
   * Returns the compound value.
   */
  operator optimet::t_uint() const;
};

#endif /* COMPOUNDITERATOR_H_ */
