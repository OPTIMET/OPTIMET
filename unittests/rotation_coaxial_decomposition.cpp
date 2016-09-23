#include "catch.hpp"

#include "RotationCoaxialDecomposition.h"
#include "Coefficients.h"
#include "Types.h"
#include "constants.h"

using namespace optimet;

TEST_CASE("Brute check") {
  using optimet::coefficient::a;
  auto const nmax = 4;
  Matrix<t_complex> input((nmax + 1) * (nmax + 1), 2);
  input.fill(0);
  SECTION("Nothing in n±1") {
    SECTION("n=1") {
      input(1, 0) = 2;
      input(1, 1) = t_complex(0, 1);

      for(auto const k: {1.0, 1.5}) {
        auto out = rotation_coaxial_decomposition(k, 2, input);
        CHECK(out(1, 0).real() == Approx(2 + k * k));
        CHECK(out(1, 0).imag() == Approx(0));
        CHECK(out(1, 1).real() == Approx(0));
        CHECK(out(1, 1).imag() == Approx(-1));
      }
    }
    SECTION("n=2") {
      input(4, 0) = 2;
      input(4, 1) = t_complex(0, 1);

      auto out = rotation_coaxial_decomposition(1, 6, input);
      CHECK(out(4, 0).real() == Approx(4));
      CHECK(out(4, 0).imag() == Approx(0));
      CHECK(out(4, 1).real() == Approx(0));
      CHECK(out(4, 1).imag() == Approx(-3));
    }
  }

  SECTION("Something in n - 1") {
    input(1, 0) = 2;
    input(1, 1) = t_complex(0, 1);

    auto out = rotation_coaxial_decomposition(1, 6, input);
    CHECK(out(5, 0).real() == Approx(3 * a<>(1, -1) * 2));
    CHECK(out(5, 0).imag() == Approx(0));
    CHECK(out(5, 1).real() == Approx(0));
    CHECK(out(5, 1).imag() == Approx(3 * a<>(1, -1)));
  }

  SECTION("Something in n + 1") {
    input(5, 0) = 2;
    input(5, 1) = t_complex(0, 1);

    for(auto const k: {1.0, 1.5}) {
      auto out = rotation_coaxial_decomposition(k, 2, input);
      CHECK(out(1, 0).real() == Approx(a<>(1, -1) * 2 * k));
      CHECK(out(1, 0).imag() == Approx(0));
      CHECK(out(1, 1).real() == Approx(0));
      CHECK(out(1, 1).imag() == Approx(a<>(1, -1) * k));
    }
  }
}
