#include <iostream>
#include "catch.hpp"

#include "FastMatrixMultiply.h"


class FastMatrixMultiply: public optimet::FastMatrixMultiply {
  public:
    using optimet::FastMatrixMultiply::FastMatrixMultiply;
    const decltype(global_indices_) & indices() const { return global_indices_; }
    const decltype(rotations_) & rotations() const { return rotations_; }
    const decltype(mie_coefficients_) & mie_coefficients() const { return mie_coefficients_; }
};

TEST_CASE("Single object") {
  using namespace optimet;
  ElectroMagnetic const elmag{13.1, 1.0};
  auto const wavenumber = 3;
  auto const radius = 500e-9;
  auto const nHarmonics = 5;
  std::vector<Scatterer> const scatterers{{{0, 0, 0}, elmag, radius, nHarmonics}};


  ::FastMatrixMultiply fmm{wavenumber, scatterers};

  SECTION("Internal constructed object sanity") {
    CHECK(fmm.indices().size() == 2);
    CHECK(fmm.rotations().size() == 1);
    CHECK(fmm.mie_coefficients().size() == 2 * nHarmonics * (nHarmonics + 2));
  }
}
