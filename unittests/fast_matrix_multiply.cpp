#include "Aliases.h"
#include "FastMatrixMultiply.h"
#include "Geometry.h"
#include "HarmonicsIterator.h"
#include "Result.h"
#include "Solver.h"
#include "Tools.h"
#include "catch.hpp"
#include <iostream>

class FastMatrixMultiply : public optimet::FastMatrixMultiply {
public:
  using optimet::FastMatrixMultiply::FastMatrixMultiply;
  const decltype(global_indices_) &indices() const { return global_indices_; }
  const decltype(rotations_) &rotations() const { return rotations_; }
  const decltype(mie_coefficients_) &mie_coefficients() const { return mie_coefficients_; }
};

ElectroMagnetic const elmag{13.1, 1.0};
auto const wavenumber = 2 * optimet::constant::pi / (1200 * 1e-9);
auto const nHarmonics = 2;
auto const radius = 500e-9;

TEST_CASE("Single object") {
  using namespace optimet;
  std::vector<Scatterer> const scatterers{{{0, 0, 0}, elmag, radius, nHarmonics}};

  ::FastMatrixMultiply fmm{wavenumber, scatterers};

  SECTION("Internal constructed object sanity") {
    CHECK(fmm.indices().size() == 2);
    CHECK(fmm.indices().front() == 0);
    CHECK(fmm.indices().back() == 2 * nHarmonics * (nHarmonics + 2));

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

TEST_CASE("Two objects") {
  using namespace optimet;
  auto const radius = 500e-9;
  Eigen::Matrix<t_real, 3, 1> const direction = Vector<t_real>::Random(3).normalized();
  std::vector<Scatterer> const scatterers{
      {{0, 0, 0}, elmag, radius, nHarmonics},
      {direction * 3 * radius, elmag, 2 * radius, nHarmonics + 1}};

  ::FastMatrixMultiply fmm{wavenumber, scatterers};

  SECTION("Internal constructed object sanity") {
    CHECK(fmm.indices().size() == 3);
    CHECK(fmm.indices()[0] == 0);
    CHECK(fmm.indices()[1] == 2 * nHarmonics * (nHarmonics + 2));
    CHECK(fmm.indices()[2] ==
          2 * nHarmonics * (nHarmonics + 2) + 2 * (nHarmonics + 1) * (nHarmonics + 3));

    CHECK(fmm.mie_coefficients().size() == fmm.indices()[2]);
    auto const n0 = fmm.indices()[1];
    auto const T0 = scatterers.front().getTLocal(wavenumber * constant::c, ElectroMagnetic());
    CHECK(fmm.mie_coefficients().head(n0).isApprox(T0));
    auto const n1 = fmm.indices()[2] - fmm.indices()[1];
    auto const T1 = scatterers.back().getTLocal(wavenumber * constant::c, ElectroMagnetic());
    CHECK(fmm.mie_coefficients().segment(n0, n1).isApprox(T1));

    CHECK(fmm.rotations().size() == 2);
    auto const r0 = fmm.rotations().front().basis_rotation();
    auto const r1 = fmm.rotations().back().basis_rotation();
    Eigen::Matrix<t_real, 3, 1> const x(0, 0, 1);
    auto const theta = std::acos(-direction(2));
    auto const phi = std::atan2(-direction(1), -direction(0));
    CHECK(fmm.rotations().front().theta() == Approx(std::acos(direction(2))));
    CHECK(fmm.rotations().back().theta() == Approx(std::acos(-direction(2))));
    CHECK(fmm.rotations().front().phi() == Approx(std::atan2(direction(1), direction(0))));
    CHECK(fmm.rotations().back().phi() == Approx(std::atan2(-direction(1), -direction(0))));
    CHECK((r0 * x).isApprox(-(r1 * x)));
    CHECK((r0.transpose() * direction).isApprox(x));
    CHECK((r0 * r1.transpose() * r0 * r1.transpose()).isApprox(Matrix<t_real>::Identity(3, 3)));
  }
}

TEST_CASE("Transparent objects") {
  using namespace optimet;
  CHECK(HarmonicsIterator::max_flat(nHarmonics) - 1 == nHarmonics * (nHarmonics + 2));
  Eigen::Matrix<t_real, 3, 1> const direction = Vector<t_real>::Random(3).normalized();
  Geometry geometry;
  geometry.pushObject({{0, 0, 0}, geometry.bground, radius, nHarmonics});
  geometry.pushObject({direction * 3 * radius * 1.00001, geometry.bground, 2 * radius, nHarmonics});

  auto const wavelength = 1490.0e-9;
  Spherical<t_real> const vKinc{2 * consPi / wavelength, 90 * consPi / 180.0, 90 * consPi / 180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  auto excitation =
      std::make_shared<Excitation>(0, Tools::toProjection(vKinc, Eaux), vKinc, nHarmonics);
  excitation->populate();
  geometry.update(excitation);

  Solver solver(&geometry, excitation, O3DSolverIndirect, nHarmonics);
  optimet::FastMatrixMultiply fmm(geometry.bground, excitation->omega() / constant::c,
                                  geometry.objects);

  Vector<t_complex> const input = Vector<t_complex>::Random(solver.Q.size());
  Vector<t_complex> const expected = solver.S * input;
  auto const actual = fmm * input;
  CHECK(expected.isApprox(actual));
  CHECK(expected.isApprox(input));
}

TEST_CASE("Standard vs Fast matrix multiply") {
  using namespace optimet;
  CHECK(HarmonicsIterator::max_flat(nHarmonics) - 1 == nHarmonics * (nHarmonics + 2));
  auto const radius = 500.0e-9;
  Eigen::Matrix<t_real, 3, 1> const direction(0, 0, 1); //= Vector<t_real>::Random(3).normalized();
  Geometry geometry;
  geometry.pushObject({{0, 0, 0}, elmag, radius, nHarmonics});
  geometry.pushObject({direction * 3 * radius * 1.500001, elmag, 2 * radius, nHarmonics});

  auto const wavelength = 1490.0e-9;
  Spherical<t_real> const vKinc{2 * consPi / wavelength, 90 * consPi / 180.0, 90 * consPi / 180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  auto excitation =
      std::make_shared<Excitation>(0, Tools::toProjection(vKinc, Eaux), vKinc, nHarmonics);
  excitation->populate();
  geometry.update(excitation);

  Solver solver(&geometry, excitation, O3DSolverIndirect, nHarmonics);
  optimet::FastMatrixMultiply fmm(geometry.bground, excitation->omega() / constant::c,
                                  geometry.objects);

  Vector<t_complex> input = Vector<t_complex>::Zero(solver.Q.size());
  input(5) = 1;
  Vector<t_complex> const expected = solver.S * input;
  auto const actual = fmm * input;
  CAPTURE(actual.transpose());
  CAPTURE(expected.transpose());
  CAPTURE((actual - expected).transpose());
  CHECK(expected.isApprox(actual));
}
