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

#include <iostream>
#include "catch.hpp"

#include "Types.h"
#include "HarmonicsIterator.h"
#include "CompoundIterator.h"

using namespace optimet;

TEST_CASE("Harmonic iterator") {
  SECTION("Creation") {
    CHECK(HarmonicsIterator().n() == 0);
    CHECK(HarmonicsIterator().m() == 0);
    CHECK(HarmonicsIterator(1).n() == 1);
    CHECK(HarmonicsIterator(1).m() == 1);
    CHECK(HarmonicsIterator(2, -1).n() == 2);
    CHECK(HarmonicsIterator(2, -1).m() == -1);

    CHECK_THROWS(HarmonicsIterator(1, -2));
    CHECK_THROWS(HarmonicsIterator(3, 4));
    CHECK_THROWS(HarmonicsIterator(0, 1));
  }
  SECTION("Comparison") {
    CHECK(HarmonicsIterator(2, 2) == HarmonicsIterator(2, 2));
    CHECK(HarmonicsIterator(1, -1) == HarmonicsIterator(1, -1));
  }
  SECTION("Increment") {
    CHECK(++HarmonicsIterator(2, 2) == HarmonicsIterator(2, 1));
    CHECK(++HarmonicsIterator(2, 1) == HarmonicsIterator(2, 0));
    CHECK(++HarmonicsIterator(2, 0) == HarmonicsIterator(2, -1));
    CHECK(++HarmonicsIterator(2, -1) == HarmonicsIterator(2, -2));
    CHECK(++HarmonicsIterator(2, -2) == HarmonicsIterator(3, 3));
  }
  SECTION("Decrement") {
    CHECK(--HarmonicsIterator(2, -2) == HarmonicsIterator(2, -1));
    CHECK(--HarmonicsIterator(2, -1) == HarmonicsIterator(2, 0));
    CHECK(--HarmonicsIterator(2, 0) == HarmonicsIterator(2, 1));
    CHECK(--HarmonicsIterator(2, 1) == HarmonicsIterator(2, 2));
    CHECK(--HarmonicsIterator(2, 2) == HarmonicsIterator(1, -1));

    CHECK(++(--HarmonicsIterator(2, 2)) == HarmonicsIterator(2, 2));
    CHECK(++(--HarmonicsIterator(2, 1)) == HarmonicsIterator(2, 1));
    CHECK(++(--HarmonicsIterator(2, -2)) == HarmonicsIterator(2, -2));
    CHECK(--(++HarmonicsIterator(2, 2)) == HarmonicsIterator(2, 2));
    CHECK(--(++HarmonicsIterator(2, 1)) == HarmonicsIterator(2, 1));
    CHECK(--(++HarmonicsIterator(2, -2)) == HarmonicsIterator(2, -2));
  }
  SECTION("Dereferencing") {
    CHECK(*HarmonicsIterator() == 0);
    CHECK(*(++HarmonicsIterator()) == 1);
    HarmonicsIterator iterator;
    t_uint i(0);
    for(; iterator != HarmonicsIterator::end(3); ++i, ++iterator) {
      CHECK(*iterator == i);
    }
    CHECK(iterator.n() == 4);
    CHECK(iterator.m() == 4);
    CHECK(i == HarmonicsIterator::max_flat(3));
  }
  SECTION("Boundaries") {
    CHECK(HarmonicsIterator::end(1) == ++HarmonicsIterator(1, -1));
    CHECK(HarmonicsIterator::least(1) == --HarmonicsIterator(1, 1));
    CHECK(HarmonicsIterator::least(0) == --HarmonicsIterator(0));
    CHECK(HarmonicsIterator::least(0) == --(--HarmonicsIterator(0)));
    // Equivalent to CompoundIterator::nMax
    CHECK(HarmonicsIterator::max_flat(2) - 1 == 2 * (2 + 2));
    CHECK(HarmonicsIterator::max_flat(3) - 1 == 3 * (3 + 2));
  }
  SECTION("C++11") {
    HarmonicsIterator a = {1};
    a = {2};
    CHECK(a == HarmonicsIterator(2, 2));
    HarmonicsIterator b = {1, 1};
    b = {1, 1};
    CHECK(b == HarmonicsIterator(1, 1));
  }
}

TEST_CASE("Comparison with CompoundIterator") {

  CHECK(static_cast<t_int>(HarmonicsIterator().n()) == CompoundIterator().first);
  CHECK(HarmonicsIterator().m() == CompoundIterator().second);
  CHECK(static_cast<t_int>(*HarmonicsIterator()) == CompoundIterator().compound + 1);

  CHECK(static_cast<t_int>(HarmonicsIterator(1, 0).n()) == CompoundIterator(1).first);
  CHECK(HarmonicsIterator(1, 0).m() == CompoundIterator(1).second);
  CHECK(static_cast<t_int>(*HarmonicsIterator(1, 0)) == CompoundIterator(1).compound + 1);

  CHECK(static_cast<t_int>(HarmonicsIterator(1, -1).n()) == CompoundIterator(2).first);
  CHECK(HarmonicsIterator(1, -1).m() == CompoundIterator(2).second);
  CHECK(static_cast<t_int>(*HarmonicsIterator(1, -1)) == CompoundIterator(2).compound + 1);

  HarmonicsIterator a;
  CompoundIterator b;
  for(; a != HarmonicsIterator::end(100); ++a, ++b)
    CHECK(static_cast<int>(*a) == static_cast<int>(b) + 1);

  for(t_int i(0); i < 5; ++i)
    CHECK(HarmonicsIterator::max_flat(i) == CompoundIterator::max(i) + 1);
}
