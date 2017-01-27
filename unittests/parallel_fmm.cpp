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

#include "Excitation.h"
#include "FMMBelosSolver.h"
#include "Geometry.h"
#include "PreconditionedMatrixSolver.h"
#include "Reader.h"
#include "catch.hpp"
#include "mpi/FastMatrixMultiply.h"
#include <BelosTypes.hpp>
#include <iostream>

TEST_CASE("ReduceComputation") {
  using namespace optimet;
  mpi::Communicator const world;
  if(world.size() < 3)
    return;

  ElectroMagnetic const silicon{13.1, 1.0};
  auto const nprocs = 3;
  auto const nHarmonics = 5;
  auto const radius = 500e-9;
  std::vector<Scatterer> const scatterers = {
      {{0, 0, 0}, silicon, radius, nHarmonics},     {{0, 0, 0}, silicon, radius, nHarmonics + 1},
      {{0, 0, 0}, silicon, radius, nHarmonics + 2}, {{0, 0, 0}, silicon, radius, nHarmonics + 3},
      {{0, 0, 0}, silicon, radius, nHarmonics + 4}, {{0, 0, 0}, silicon, radius, nHarmonics + 5},
  };

  auto const distribution = mpi::details::vector_distribution(scatterers.size(), nprocs);
  Matrix<bool> locals = Matrix<bool>::Zero(scatterers.size(), scatterers.size());
  locals.topLeftCorner(3, 3).fill(true);
  locals.bottomRightCorner(3, 3).fill(true);
  mpi::GraphCommunicator const graph_comm(
      world, mpi::details::graph_edges(locals.array() == false, distribution), false);

  auto const nfunctions = [&scatterers](t_int i) {
    return 2 * scatterers[i].nMax * (scatterers[i].nMax + 2);
  };
  int const sizes[] = {nfunctions(0) + nfunctions(1) + nfunctions(2),
                       nfunctions(3) + nfunctions(4) + nfunctions(5)};
  Vector<int> const messages[] = {
      Vector<int>::LinSpaced(sizes[1], sizes[0], sizes[1]).eval(),
      10 * Vector<int>::LinSpaced(sizes[0] + sizes[1], 0, sizes[0] + sizes[1]).eval(),
      100 * Vector<int>::LinSpaced(sizes[0], 0, sizes[0]).eval(), Vector<int>::Zero(0).eval()};

  mpi::FastMatrixMultiply::ReduceComputation reduction(graph_comm, locals.array() == false,
                                                       distribution, scatterers);
  std::vector<int> const buffer_sizes{2 * (nfunctions(0) + nfunctions(1)),
                                      2 * (nfunctions(2) + nfunctions(3)),
                                      2 * (nfunctions(4) + nfunctions(5)), 0};
  SECTION("Send message") {
    Vector<int> in_buffer, out_buffer;
    if(auto const request =
           reduction.send(messages[std::min<int>(world.rank(), nprocs)], in_buffer, out_buffer))
      REQUIRE(request);
    else
      REQUIRE(false);

    switch(world.rank()) {
    case 0:
      REQUIRE(out_buffer.size() == 2 * (nfunctions(0) + nfunctions(1)));
      CHECK(out_buffer.head(nfunctions(0) + nfunctions(1)) ==
            messages[1].head(nfunctions(0) + nfunctions(1)));
      CHECK(out_buffer.tail(nfunctions(0) + nfunctions(1)) ==
            messages[2].head(nfunctions(0) + nfunctions(1)));
      break;
    case 1:
      REQUIRE(out_buffer.size() == 2 * (nfunctions(2) + nfunctions(3)));
      CHECK(out_buffer.head(nfunctions(3)) == messages[0].head(nfunctions(3)));
      CHECK(out_buffer.segment(nfunctions(3), nfunctions(2) + nfunctions(3)) ==
            messages[1].segment(nfunctions(0) + nfunctions(1), nfunctions(2) + nfunctions(3)));
      CHECK(out_buffer.tail(nfunctions(2)) == messages[2].tail(nfunctions(2)));
      break;
    case 2:
      REQUIRE(out_buffer.size() == 2 * (nfunctions(4) + nfunctions(5)));
      CHECK(out_buffer.head(nfunctions(4) + nfunctions(5)) ==
            messages[0].tail(nfunctions(4) + nfunctions(5)));
      CHECK(out_buffer.tail(nfunctions(4) + nfunctions(5)) ==
            messages[1].tail(nfunctions(4) + nfunctions(5)));
      break;

    default:
      REQUIRE(out_buffer.size() == 0);
    }
  }

  SECTION("Do Reduction") {
    std::vector<int> const result_sizes{nfunctions(0) + nfunctions(1),
                                        nfunctions(2) + nfunctions(3),
                                        nfunctions(4) + nfunctions(5), 0};
    auto const rank = std::min<int>(world.rank(), nprocs);
    Vector<int> const buffer =
        (world.rank() + 2) * Vector<int>::LinSpaced(buffer_sizes[rank], 0, buffer_sizes[rank]);
    Vector<int> result(result_sizes[rank]);
    result.fill(world.rank() + 1);
    reduction.reduce(result, buffer);

    switch(world.rank()) {
    case 0:
      CHECK(result ==
            result.Constant(result.size(), world.rank() + 1) +
                buffer.head(nfunctions(0) + nfunctions(1)) +
                buffer.tail(nfunctions(0) + nfunctions(1)));
      break;
    case 1:
      CHECK(result.head(nfunctions(2)) ==
            result.Constant(nfunctions(2), world.rank() + 1) +
                buffer.segment(nfunctions(3), nfunctions(2)) + buffer.tail(nfunctions(2)));
      CHECK(result.tail(nfunctions(3)) ==
            result.Constant(nfunctions(3), world.rank() + 1) +
                buffer.segment(nfunctions(3) + nfunctions(2), nfunctions(3)) +
                buffer.head(nfunctions(3)));
      break;
    case 2:
      CHECK(result ==
            result.Constant(result.size(), world.rank() + 1) +
                buffer.head(nfunctions(4) + nfunctions(5)) +
                buffer.tail(nfunctions(4) + nfunctions(5)));
      break;
    default:
      REQUIRE(result.size() == 0);
    }
  }
}

TEST_CASE("DistributeInput") {
  using namespace optimet;
  mpi::Communicator const world;
  if(world.size() < 3)
    return;

  ElectroMagnetic const silicon{13.1, 1.0};
  auto const nprocs = 3;
  auto const nHarmonics = 5;
  auto const radius = 500e-9;
  std::vector<Scatterer> const scatterers = {
      {{0, 0, 0}, silicon, radius, nHarmonics},     {{0, 0, 0}, silicon, radius, nHarmonics + 1},
      {{0, 0, 0}, silicon, radius, nHarmonics + 2}, {{0, 0, 0}, silicon, radius, nHarmonics + 3},
      {{0, 0, 0}, silicon, radius, nHarmonics + 4}, {{0, 0, 0}, silicon, radius, nHarmonics + 5},
  };

  auto const vector_distribution = mpi::details::vector_distribution(scatterers.size(), nprocs);
  Matrix<bool> locals = Matrix<bool>::Zero(scatterers.size(), scatterers.size());
  locals.topLeftCorner(3, 3).fill(true);
  locals.bottomRightCorner(3, 3).fill(true);
  mpi::GraphCommunicator const graph_comm(
      world, mpi::details::graph_edges(locals.array() == false, vector_distribution), false);

  auto const nfunctions = [&scatterers](t_int i) {
    return 2 * scatterers[i].nMax * (scatterers[i].nMax + 2);
  };
  int const sizes[] = {nfunctions(0) + nfunctions(1), nfunctions(2) + nfunctions(3),
                       nfunctions(4) + nfunctions(5), 0};
  Vector<int> const messages[] = {Vector<int>::LinSpaced(sizes[0], 0, sizes[0]).eval(),
                                  10 * Vector<int>::LinSpaced(sizes[1], 0, sizes[1]).eval(),
                                  100 * Vector<int>::LinSpaced(sizes[2], 0, sizes[2]).eval(),
                                  Vector<int>::Zero(0).eval()};

  mpi::FastMatrixMultiply::DistributeInput distribution(graph_comm, locals.array() == false,
                                                        vector_distribution, scatterers);
  std::vector<int> const buffer_sizes{sizes[1] + sizes[2], sizes[0] + sizes[1] + sizes[2],
                                      sizes[0] + sizes[1], 0};
  auto const rank = std::min<int>(world.rank(), nprocs);
  SECTION("Send message") {
    Vector<int> buffer;
    if(auto const request = distribution.send(messages[rank], buffer))
      REQUIRE(request);
    else
      REQUIRE(false);

    REQUIRE(buffer.size() == buffer_sizes[rank]);
    switch(world.rank()) {
    case 0:
      CHECK(buffer.head(sizes[1]) == messages[1]);
      CHECK(buffer.tail(sizes[2]) == messages[2]);
      break;
    case 1:
      CHECK(buffer.head(sizes[0]) == messages[0]);
      CHECK(buffer.tail(sizes[2]) == messages[2]);
      break;
    case 2:
      CHECK(buffer.head(sizes[0]) == messages[0]);
      CHECK(buffer.tail(sizes[1]) == messages[1]);
      break;
    default:
      break;
    }
  }

  SECTION("Synthesize output") {
    Vector<int> const buffer =
        Vector<int>::LinSpaced(buffer_sizes[rank], 1, buffer_sizes[rank] + 1);
    Vector<int> const result;
    distribution.synthesize(buffer, result);
    switch(world.rank()) {
    case 0:
      REQUIRE(result.size() == nfunctions(3) + nfunctions(4) + nfunctions(5));
      CHECK(result == buffer.tail(result.size()));
      break;
    case 1:
      REQUIRE(result.size() ==
              nfunctions(0) + nfunctions(1) + nfunctions(2) + nfunctions(3) + nfunctions(4) +
                  nfunctions(5));
      CHECK(result == buffer);
      break;
    case 2:
      REQUIRE(result.size() == nfunctions(0) + nfunctions(1) + nfunctions(2));
      CHECK(result == buffer.head(nfunctions(0) + nfunctions(1) + nfunctions(2)));
      break;
    default:
      REQUIRE(result.size() == 0);
      break;
    }
  }
}

template <class T0, class T1>
optimet::Vector<typename T1::Scalar>
split(std::vector<Scatterer> const &scatterers, Eigen::DenseBase<T0> const &distribution,
      Eigen::DenseBase<T1> const &input) {
  assert(distribution.size() == scatterers.size());
  int N(0);
  for(int i(0); i < distribution.size(); ++i)
    if(distribution(i))
      N += 2 * scatterers[i].nMax * (scatterers[i].nMax + 2);

  optimet::Vector<typename T1::Scalar> result(N);
  for(int i(0), j(0), k(0); i < distribution.size(); ++i) {
    auto const size = 2 * scatterers[i].nMax * (scatterers[i].nMax + 2);
    if(distribution(i)) {
      assert(input.size() >= k + size);
      assert(result.size() >= j + size);
      result.segment(j, size) = input.segment(k, size);
      j += size;
    }
    k += size;
  }
  return result;
}

TEST_CASE("MPI vs serial FMM") {
  using namespace optimet;
  mpi::Communicator const world;
  // spherical coords, ε, μ, radius, nmax
  int const nHarmonics = 3;
  std::vector<Scatterer> scatterers;
  for(int i(0); i < 10; ++i)
    scatterers.push_back({{t_real(i) * 1.5 * 2e-6, 0, 0},
                          {0.45e0 + 0.1 * t_real(i), 1.1e0},
                          (0.5 + 0.01 * t_real(i)) * 2e-6,
                          std::min(10, nHarmonics + i)});
  auto const nscatt = scatterers.size();

  // Create excitation
  auto const wavelength = 14960e-9;
  auto const wavenumber = (2 * constant::pi / wavelength) / constant::c;
  auto const distribution = mpi::details::vector_distribution(nscatt, world.size());
  FastMatrixMultiply serial(wavenumber, scatterers);
  auto const serial_input =
      world.broadcast<Vector<t_complex>>(Vector<t_complex>::Random(serial.cols()));
  CHECK(serial_input.size() == serial.cols());
  auto const serial_output =
      split(scatterers, distribution.array() == world.rank(), serial(serial_input));
  auto const transpose_serial_output =
      split(scatterers, distribution.array() == world.rank(), serial.transpose(serial_input));

  if(world.size() < 2)
    return;

  auto const parallel_input = split(scatterers, distribution.array() == world.rank(), serial_input);
  // -1 corresponds to computations using only local data
  // nscatt - 1 corresponds to computations using all data gathered from all procs
  for(int diag(-1); diag < static_cast<int>(nscatt); ++diag) {
    SECTION("Overlay communication and data with diag = " + std::to_string(diag)) {
      mpi::FastMatrixMultiply parallel(wavenumber, scatterers, diag, distribution, world);
      auto const parallel_out = parallel(parallel_input);
      REQUIRE(parallel_out.size() == serial_output.size());
      CHECK(parallel_out.isApprox(serial_output));

      auto const transpose_parallel_out = parallel.transpose(parallel_input);
      REQUIRE(transpose_parallel_out.size() == transpose_serial_output.size());
      CHECK(transpose_parallel_out.isApprox(transpose_serial_output));
    }
  }

  SECTION("Randomly patterned communications") {
    Matrix<bool> locals = Matrix<bool>::Random(nscatt, nscatt);

    mpi::FastMatrixMultiply parallel(wavenumber, scatterers, locals, distribution, world);
    auto const parallel_out = parallel(parallel_input);
    REQUIRE(parallel_out.size() == serial_output.size());
    CHECK(parallel_out.isApprox(serial_output));

    auto const transpose_parallel_out = parallel.transpose(parallel_input);
    REQUIRE(transpose_parallel_out.size() == transpose_serial_output.size());
    CHECK(transpose_parallel_out.isApprox(transpose_serial_output));
  }

  SECTION("Conjugate operation") {
    mpi::FastMatrixMultiply parallel(wavenumber, scatterers, world);

    auto const expected =
        split(scatterers, distribution.array() == world.rank(), serial.conjugate(serial_input));
    auto const actual = parallel.conjugate(parallel_input);
    CHECK(actual.isApprox(expected));
  }

  SECTION("Adjoint operation") {
    mpi::FastMatrixMultiply parallel(wavenumber, scatterers, world);

    auto const expected =
        split(scatterers, distribution.array() == world.rank(), serial.adjoint(serial_input));
    auto const actual = parallel.adjoint(parallel_input);
    CHECK(actual.isApprox(expected));
  }
}

TEST_CASE("FMM solver vs serial solver") {
  using namespace optimet;
  auto const nHarmonics = 5;
  auto const nSpheres = 5;
  auto geometry = std::make_shared<Geometry>();
  // spherical coords, ε, μ, radius, nmax
  for(t_uint i(0); i < nSpheres; ++i)
    geometry->pushObject(
        {{static_cast<t_real>(i) * 1.5 * 2e-6, 0, 0}, {5e0, 1.1e0}, 0.5 * 2e-6, nHarmonics});

  // Create excitation
  auto const wavelength = 14960e-9;
  Spherical<t_real> const vKinc{2 * consPi / wavelength, 90 * consPi / 180.0, 90 * consPi / 180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  auto const excitation =
      std::make_shared<Excitation>(0, Tools::toProjection(vKinc, Eaux), vKinc, nHarmonics);
  excitation->populate();
  geometry->update(excitation);

  optimet::mpi::Communicator world;
  optimet::Result parallel(geometry, excitation);
  optimet::solver::FMMBelos solver(geometry, excitation, world);
  solver.belos_parameters()->set("Solver", "GMRES");
  solver.belos_parameters()->set<int>("Num Blocks", 500);
  solver.belos_parameters()->set("Maximum Iterations", 4000);
  solver.belos_parameters()->set("Convergence Tolerance", 1.0e-10);
  solver.solve(parallel.scatter_coef, parallel.internal_coef);

  optimet::Result serial(geometry, excitation);
  optimet::solver::PreconditionedMatrix const serial_solver(geometry, excitation, world);
  serial_solver.solve(serial.scatter_coef, serial.internal_coef);

  REQUIRE(parallel.scatter_coef.rows() == serial.scatter_coef.rows());
  REQUIRE(parallel.scatter_coef.cols() == serial.scatter_coef.cols());
  auto const scatter_tol = 1e-6 * std::max(1., serial.scatter_coef.array().abs().maxCoeff());
  CHECK(parallel.scatter_coef.isApprox(serial.scatter_coef, scatter_tol));
  REQUIRE(parallel.internal_coef.rows() == serial.internal_coef.rows());
  REQUIRE(parallel.internal_coef.cols() == serial.internal_coef.cols());
  auto const internal_tol = 1e-6 * std::max(1., serial.internal_coef.array().abs().maxCoeff());
  CHECK(parallel.internal_coef.isApprox(serial.internal_coef, internal_tol));
}

TEST_CASE("Parallel matrix vs serial matrix") {
  using namespace optimet;
  mpi::Communicator const world;
  auto const nHarmonics = 5;
  auto const nSpheres = 8;
  auto geometry = std::make_shared<Geometry>();
  // spherical coords, ε, μ, radius, nmax
  for(t_uint i(0); i < nSpheres; ++i)
    geometry->pushObject(
        {{static_cast<t_real>(i) * 1.5 * 2e-6, 0, 0}, {0.45e0, 1.1e0}, 0.5 * 2e-6, nHarmonics});

  // Create excitation
  auto const wavelength = 14960e-9;
  Spherical<t_real> const vKinc{2 * consPi / wavelength, 90 * consPi / 180.0, 90 * consPi / 180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  auto const excitation =
      std::make_shared<Excitation>(0, Tools::toProjection(vKinc, Eaux), vKinc, nHarmonics);
  excitation->populate();
  geometry->update(excitation);

  auto const serial = preconditioned_scattering_matrix(*geometry, excitation);

  auto const diags = std::max<int>(1, geometry->objects.size() / 2 - 2);
  mpi::FastMatrixMultiply parallel(geometry->bground, excitation->wavenumber(), geometry->objects,
                                   diags, world);
  auto const distribution =
      mpi::details::vector_distribution(geometry->objects.size(), world.size());
  auto const size = serial.cols();

  for(t_int i(0); i < size; ++i) {
    Vector<t_complex> const serial_input = Vector<t_complex>::Unit(size, i);
    Vector<t_complex> const expected = serial * serial_input;
    Vector<t_complex> const parallel_input =
        split(geometry->objects, distribution.array() == world.rank(), serial_input);
    Vector<t_complex> const actual = parallel * parallel_input;
    auto const gathered = world.all_gather(actual);
    CHECK(expected.isApprox(gathered));
  }
}

TEST_CASE("Read XML") {
  std::istringstream buffer(
      "<simulation>\n"
      "  <harmonics nmax=\"6\" />\n"
      "</simulation>\n"
      "<source type=\"planewave\">\n"
      "  <wavelength value=\"1460\" />\n"
      "  <propagation theta=\"90\" phi=\"90\" />\n"
      "  <polarization Etheta.real=\"1.0\" Etheta.imag=\"0.0\" Ephi.real=\"0.0\" "
      "Ephi.imag=\"0.0\" />\n"
      "</source>\n"
      "<geometry>\n"
      "  <object type=\"sphere\">\n"
      "    <cartesian x=\"0.0\" y=\"0.0\" z=\"0.0\" />\n"
      "    <properties radius=\"500.0\" />\n"
      "    <epsilon type=\"relative\" value.real=\"13.0\" value.imag=\"0.0\" />\n"
      "    <mu type=\"relative\" value.real=\"1.0\" value.imag=\"0.0\" />\n"
      "  </object>\n"
      "</geometry>\n"
      "<output type=\"field\">\n"
      "  <grid type=\"cartesian\">\n"
      "    <x min=\"-1000\" max=\"1000\" steps=\"21\" />\n"
      "    <y min=\"-0.1\" max=\"0.1\" steps=\"2\" />\n"
      "    <z min=\"-1000\" max=\"1000\" steps=\"21\" />\n"
      "  </grid>\n"
      "</output>\n"
      "<FMM subdiagonals=\"2\"/>\n"
      "<ParameterList name=\"Belos\">\n"
      "  <Parameter name=\"Solver\" type=\"string\" value=\"GMRES\"/>\n"
      "  <Parameter name=\"Maximum Iterations\" type=\"int\" value=\"4000\"/>\n"
      "  <Parameter name=\"Output Frequency\" type=\"int\" value=\"20\"/>\n"
      "</ParameterList>\n");

  auto const run = optimet::simulation_input(buffer);
  CHECK(run.do_fmm);
  CHECK(run.fmm_subdiagonals == 2);
  auto const solver = optimet::solver::factory(run);
  CHECK_NOTHROW(std::dynamic_pointer_cast<optimet::solver::FMMBelos>(solver));
}
