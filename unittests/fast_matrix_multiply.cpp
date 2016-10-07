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
  auto const wavenumber = 2 * constant::pi / (1200 * 1e-9);
  auto const radius = 500e-9;
  auto const nHarmonics = 5;
  std::vector<Scatterer> const scatterers{{{0, 0, 0}, elmag, radius, nHarmonics}};


  ::FastMatrixMultiply fmm{wavenumber, scatterers};

  SECTION("Internal constructed object sanity") {
    CHECK(fmm.indices().size() == 2);
    CHECK(fmm.indices().front() == 0);
    CHECK(fmm.indices().back() == nHarmonics * (nHarmonics + 2));

    CHECK(fmm.mie_coefficients().size() == 2 * nHarmonics * (nHarmonics + 2));
    CHECK(fmm.mie_coefficients().isApprox(
          scatterers.front().getTLocal(wavenumber * constant::c, ElectroMagnetic())));

    CHECK(fmm.rotations().size() == 0);
  }

  SECTION("Matrix is identity") {
    Vector<t_complex> const input = Vector<t_complex>::Random(2 * nHarmonics * (nHarmonics + 2));
    auto const result = fmm(input);
    CHECK(result.isApprox(input));
  }
}
