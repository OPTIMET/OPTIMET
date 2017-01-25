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

#include "AuxCoefficients.h"

using namespace optimet;

TEST_CASE("AuxCoefficients") {

  const Spherical<double> R(0.0, 1.57079632679, 1.57079632679); // π/2

  std::vector<double> Wigner, dWigner;

  SECTION("nMax == 0") {
    std::tie(Wigner, dWigner) = AuxCoefficients::VIGdVIG(0, 0, R);

    CHECK(Wigner.size() == 1);
    CHECK(dWigner.size() == 1);
  }

  SECTION("nMax > 0, m = nMax") {
    std::tie(Wigner, dWigner) = AuxCoefficients::VIGdVIG(7, 7, R);

    CHECK(Wigner.size() == 8);
    CHECK(dWigner.size() == 8);

    CHECK(Wigner[0] == 0.0);
    CHECK(Wigner[1] == 0.0);
    CHECK(Wigner[2] == 0.0);
    CHECK(Wigner[3] == 0.0);
    CHECK(Wigner[4] == 0.0);
    CHECK(Wigner[5] == 0.0);
    CHECK(Wigner[6] == 0.0);
    CHECK(Wigner[7] == Approx(0.457681828621));

    CHECK(dWigner[0] == 0.0);
    CHECK(dWigner[1] == 0.0);
    CHECK(dWigner[2] == 0.0);
    CHECK(dWigner[3] == 0.0);
    CHECK(dWigner[4] == 0.0);
    CHECK(dWigner[5] == 0.0);
    CHECK(dWigner[6] == 0.0);
    CHECK(dWigner[7] == Approx(1.56875582046e-11));
  }

  SECTION("nMax > 0, m = -nMax") {
    std::tie(Wigner, dWigner) = AuxCoefficients::VIGdVIG(7, -7, R);

    CHECK(Wigner.size() == 8);
    CHECK(dWigner.size() == 8);

    CHECK(Wigner[0] == 0.0);
    CHECK(Wigner[1] == 0.0);
    CHECK(Wigner[2] == 0.0);
    CHECK(Wigner[3] == 0.0);
    CHECK(Wigner[4] == 0.0);
    CHECK(Wigner[5] == 0.0);
    CHECK(Wigner[6] == 0.0);
    CHECK(Wigner[7] == Approx(-0.457681828621));

    CHECK(dWigner[0] == 0.0);
    CHECK(dWigner[1] == 0.0);
    CHECK(dWigner[2] == 0.0);
    CHECK(dWigner[3] == 0.0);
    CHECK(dWigner[4] == 0.0);
    CHECK(dWigner[5] == 0.0);
    CHECK(dWigner[6] == 0.0);
    CHECK(dWigner[7] == Approx(-1.56875582046e-11));
  }

  SECTION("nMax > 0, m = 0") {
    std::tie(Wigner, dWigner) = AuxCoefficients::VIGdVIG(7, 0, R);

    CHECK(Wigner.size() == 8);
    CHECK(dWigner.size() == 8);

    CHECK(Wigner[0] == Approx(1));
    CHECK(Wigner[1] == Approx(4.89658886015e-12));
    CHECK(Wigner[2] == Approx(-0.5));
    CHECK(Wigner[3] == Approx(-7.34488329022e-12));
    CHECK(Wigner[4] == Approx(0.375));
    CHECK(Wigner[5] == Approx(9.18110411278e-12));
    CHECK(Wigner[6] == Approx(-0.3125));
    CHECK(Wigner[7] == Approx(-1.07112881316e-11));

    CHECK(dWigner[0] == Approx(0));
    CHECK(dWigner[1] == Approx(-1));
    CHECK(dWigner[2] == Approx(-1.46897665804e-11));
    CHECK(dWigner[3] == Approx(1.5));
    CHECK(dWigner[4] == Approx(3.67244164511e-11));
    CHECK(dWigner[5] == Approx(-1.875));
    CHECK(dWigner[6] == Approx(-6.42677287894e-11));
    CHECK(dWigner[7] == Approx(2.1875));
  }
}
