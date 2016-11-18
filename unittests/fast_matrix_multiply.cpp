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

ElectroMagnetic const silicon{13.1, 1.0};
auto const wavenumber = 2 * optimet::constant::pi / (1200 * 1e-9);
auto const nHarmonics = 5;
auto const radius = 500e-9;

TEST_CASE("Single object") {
  using namespace optimet;
  std::vector<Scatterer> const scatterers{{{0, 0, 0}, silicon, radius, nHarmonics}};

  ::FastMatrixMultiply fmm{wavenumber, scatterers};

  SECTION("Internal constructed object sanity") {
    CHECK(fmm.indices().size() == 2);
    CHECK(fmm.indices().front() == 0);
    CHECK(fmm.indices().back() == 2 * nHarmonics * (nHarmonics + 2));

    CHECK(fmm.mie_coefficients().size() == 2 * nHarmonics * (nHarmonics + 2));
    CHECK(fmm.mie_coefficients().isApprox(
        scatterers.front().getTLocal(wavenumber * constant::c, ElectroMagnetic())));

    CHECK(fmm.rotations().size() == 1);
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
      {{0, 0, 0}, silicon, radius, nHarmonics},
      {direction * 3 * radius, silicon, 2 * radius, nHarmonics + 1}};

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

    CHECK(fmm.rotations().size() == 4);
    auto const r0 = fmm.rotations()[1].basis_rotation();
    auto const r1 = fmm.rotations()[2].basis_rotation();
    Eigen::Matrix<t_real, 3, 1> const x(0, 0, 1);
    auto const theta = std::acos(-direction(2));
    auto const phi = std::atan2(-direction(1), -direction(0));
    CHECK(fmm.rotations()[1].theta() == Approx(std::acos(direction(2))));
    CHECK(fmm.rotations()[2].theta() == Approx(std::acos(-direction(2))));
    CHECK(fmm.rotations()[1].phi() == Approx(std::atan2(direction(1), direction(0))));
    CHECK(fmm.rotations()[2].phi() == Approx(std::atan2(-direction(1), -direction(0))));
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

TEST_CASE("Transpose/conjugate/adjoint of the fast matrix multiply") {
  using namespace optimet;
  auto const radius = 500.0e-9;
  Eigen::Matrix<t_real, 3, 1> const direction(2, 1, 1); // = Vector<t_real>::Random(3).normalized();
  Geometry geometry;
  geometry.pushObject({{0, 0, 0}, silicon, radius, nHarmonics});
  geometry.pushObject({direction * 3 * radius * 1.500001, silicon, 2 * radius, nHarmonics});

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

  auto const size = geometry.objects.size() * nHarmonics * (nHarmonics + 2) * 2;
  Matrix<t_complex> expected(size, size), transpose(size, size), conjugate(size, size),
      adjoint(size, size);
  for(t_int i(0); i < size; ++i) {
    expected.col(i) = fmm(Vector<t_complex>::Unit(size, i));
    transpose.col(i) = fmm.transpose(Vector<t_complex>::Unit(size, i));
    conjugate.col(i) = fmm.conjugate(Vector<t_complex>::Unit(size, i));
    adjoint.col(i) = fmm.adjoint(Vector<t_complex>::Unit(size, i));
  }
  CAPTURE(transpose.row(0));
  CAPTURE(expected.col(0).transpose());
  CHECK(transpose.isApprox(expected.transpose()));
  CHECK(conjugate.isApprox(expected.conjugate()));
  CHECK(adjoint.isApprox(expected.adjoint()));
}

TEST_CASE("Standard vs Fast matrix multiply") {
  using namespace optimet;
  auto const radius = 500.0e-9;
  Eigen::Matrix<t_real, 3, 1> const direction = Vector<t_real>::Random(3).normalized();
  Geometry geometry;
  geometry.pushObject({{0, 0, 0}, silicon, radius, nHarmonics});
  geometry.pushObject({direction * 3 * radius * 1.500001, silicon, 2 * radius, nHarmonics});

  auto const wavelength = 1490.0e-9;
  Spherical<t_real> const vKinc{2 * consPi / wavelength, 90 * consPi / 180.0, 90 * consPi / 180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  auto excitation =
      std::make_shared<Excitation>(0, Tools::toProjection(vKinc, Eaux), vKinc, nHarmonics);
  excitation->populate();
  geometry.update(excitation);

  SECTION("Two particles") {
    Solver solver(&geometry, excitation, O3DSolverIndirect, nHarmonics);
    optimet::FastMatrixMultiply fmm(geometry.bground, excitation->omega() / constant::c,
                                    geometry.objects);

    auto const N = nHarmonics * (nHarmonics + 2);
    for(t_int i(0); i < solver.Q.size(); ++i) {
      auto const input = Vector<t_complex>::Unit(solver.Q.size(), i);
      Vector<t_complex> const expected = solver.S * input;
      Vector<t_complex> const actual = fmm * input;
      CHECK(expected.isApprox(actual));
    }
  }

  SECTION("Three particles") {
    ElectroMagnetic const other{9, 2.0};
    auto const x = Eigen::Matrix<t_real, 3, 1>::Unit(0).eval();
    geometry.pushObject(
        {direction * 1.5 * radius * 1.500001 + x * radius * 8, other, 0.5 * radius, nHarmonics});
    Solver solver(&geometry, excitation, O3DSolverIndirect, nHarmonics);
    optimet::FastMatrixMultiply fmm(geometry.bground, excitation->omega() / constant::c,
                                    geometry.objects);

    auto const N = nHarmonics * (nHarmonics + 2);
    for(t_int i(0); i < solver.Q.size(); ++i) {
      auto const input = Vector<t_complex>::Unit(solver.Q.size(), i);
      Vector<t_complex> const expected = solver.S * input;
      Vector<t_complex> const actual = fmm * input;
      CHECK(expected.isApprox(actual));
    }
  }
}

TEST_CASE("Ranged matrices") {
  using namespace optimet;
  auto const wavelength = 1490.0e-9;
  auto const radius = 500.0e-9;
  Eigen::Matrix<t_real, 3, 1> const direction = Vector<t_real>::Random(3).normalized();
  ElectroMagnetic const other{9, 2.0};

  std::vector<Scatterer> scatterers;
  scatterers.emplace_back(Vector<t_real>::Zero(3), silicon, radius, nHarmonics);
  scatterers.emplace_back(direction * 3 * radius * 1.500001, silicon, 2 * radius, nHarmonics);
  auto const x = Eigen::Matrix<t_real, 3, 1>::Unit(0).eval();
  scatterers.emplace_back(direction * 1.5 * radius * 1.500001 + x * radius * 8, other, 0.5 * radius,
                          nHarmonics);

  auto const omega = 2 * constant::pi / wavelength;
  auto const wavenumber = omega / constant::c;
  ElectroMagnetic const bground;

  auto const N = nHarmonics * (nHarmonics + 2);
  optimet::FastMatrixMultiply whole(wavenumber, scatterers);
  optimet::FastMatrixMultiply partA(wavenumber, scatterers, {0, 1}, {0, 1});
  optimet::FastMatrixMultiply partB(wavenumber, scatterers, {0, 1}, {1, 3});
  optimet::FastMatrixMultiply partC(wavenumber, scatterers, {1, 3}, {0, 3});

  CHECK(whole.rows() == whole.cols());
  CHECK(whole.rows() == 2 * N * scatterers.size());
  CHECK(partA.rows() == partA.cols());
  CHECK(partA.rows() == 2 * N);
  CHECK(partB.cols() == 2 * N);
  CHECK(partB.rows() == 2 * N * 2);
  CHECK(partC.cols() == 2 * N * 2);
  CHECK(partC.rows() == 2 * N * 3);

  SECTION("Transpose, conjugate, adjoint") {
    for(auto const &part : {partA, partB, partC}) {
      Matrix<t_complex> direct(part.rows(), part.cols()), transpose(part.cols(), part.rows()),
          adjoint(part.cols(), part.rows()), conjugate(part.rows(), part.cols());
      for(t_int i(0); i < part.cols(); ++i) {
        direct.col(i) = part(Vector<t_complex>::Unit(part.cols(), i));
        conjugate.col(i) = part.conjugate(Vector<t_complex>::Unit(part.cols(), i));
      }
      for(t_int i(0); i < part.rows(); ++i) {
        transpose.col(i) = part.transpose(Vector<t_complex>::Unit(part.rows(), i));
        adjoint.col(i) = part.adjoint(Vector<t_complex>::Unit(part.rows(), i));
      }
      CHECK(transpose.isApprox(direct.transpose()));
      CHECK(conjugate.isApprox(direct.conjugate()));
      CHECK(adjoint.isApprox(direct.adjoint()));
    }
  }

  SECTION("Whole against parts") {
    Matrix<t_complex> matrix_whole(whole.rows(), whole.cols()),
        matrix_A(partA.rows(), partA.cols()), matrix_B(partB.rows(), partB.cols()),
        matrix_C(partC.rows(), partC.cols());
    for(t_uint i(0); i < whole.rows(); ++i) {
      matrix_whole.col(i) = whole(Vector<t_complex>::Unit(whole.cols(), i));
      if(i < matrix_A.cols())
        matrix_A.col(i) = partA(Vector<t_complex>::Unit(partA.cols(), i));
      if(i < matrix_B.cols())
        matrix_B.col(i) = partB(Vector<t_complex>::Unit(partB.cols(), i));
      if(i < matrix_C.cols())
        matrix_C.col(i) = partC(Vector<t_complex>::Unit(partC.cols(), i));
    }

    Matrix<t_complex> reconstructed = Matrix<t_complex>::Zero(whole.rows(), whole.size());
    reconstructed.block(0, 0, partA.rows(), partA.cols()) = matrix_A;
    reconstructed.block(partA.rows(), 0, partB.rows(), partB.cols()) = matrix_B;
    reconstructed.block(0, partA.cols(), partC.rows(), partC.cols()) = matrix_C;
    CHECK(reconstructed.isApprox(matrix_whole));
  }
}
