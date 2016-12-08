#include "catch.hpp"
#include "mpi/FastMatrixMultiply.h"
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
      world, mpi::GraphCommunicator::symmetrize(
                 mpi::details::non_local_graph_edges(locals.array() == false, distribution)),
      false);

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
      world, mpi::GraphCommunicator::symmetrize(
                 mpi::details::non_local_graph_edges(locals.array() == false, vector_distribution)),
      false);

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

  if(world.size() < 2)
    return;

  auto const parallel_input = split(scatterers, distribution.array() == world.rank(), serial_input);
  SECTION("Only reduction -- all computation use local data") {
    mpi::FastMatrixMultiply parallel(wavenumber, scatterers, Matrix<bool>::Ones(nscatt, nscatt),
                                     distribution, world);
    auto const parallel_out = parallel(parallel_input);
    REQUIRE(parallel_out.size() == serial_output.size());
    CHECK(parallel_out.isApprox(serial_output));
  }

  SECTION("Only distribution -- all computation use non-local (and local) data") {
    mpi::FastMatrixMultiply parallel(wavenumber, scatterers, Matrix<bool>::Zero(nscatt, nscatt),
                                     distribution, world);
    auto const parallel_out = parallel(parallel_input);
    REQUIRE(parallel_out.size() == serial_output.size());
    CHECK(parallel_out.isApprox(serial_output));
  }

  for(int diag(1); diag < nscatt - 1; ++diag) {
    SECTION("Overlay communication and data with diag = " + std::to_string(diag)) {
      Matrix<bool> locals = Matrix<bool>::Zero(nscatt, nscatt);
      for(t_int d(0); d < nscatt; ++d) {
        auto const start = std::max(0, d - diag);
        auto const end = std::min<t_int>(nscatt, d + diag + 1);
        locals.row(d).segment(start, end - start).fill(true);
      }

      mpi::FastMatrixMultiply parallel(wavenumber, scatterers, locals, distribution, world);
      auto const parallel_out = parallel(parallel_input);
      REQUIRE(parallel_out.size() == serial_output.size());
      CHECK(parallel_out.isApprox(serial_output));
    }
  }

  SECTION("Randomly patterned communications") {
    Matrix<bool> locals = Matrix<bool>::Random(nscatt, nscatt);
    for(t_int d(0); d < nscatt; ++d)
      locals.row(d).tail(nscatt - d) = locals.col(d).tail(nscatt - d).eval();

    mpi::FastMatrixMultiply parallel(wavenumber, scatterers, locals, distribution, world);
    auto const parallel_out = parallel(parallel_input);
    REQUIRE(parallel_out.size() == serial_output.size());
    CHECK(parallel_out.isApprox(serial_output));
  }
}
