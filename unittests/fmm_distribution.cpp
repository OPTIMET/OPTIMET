#include "catch.hpp"
#include "mpi/FastMatrixMultiply.h"
#include <iostream>

TEST_CASE("Check distributions") {
  using namespace optimet;
  SECTION("Column distributions") {
    using optimet::mpi::details::vector_distribution;
    CHECK(vector_distribution(0, 4).size() == 0);

    CHECK(vector_distribution(1, 4) == Vector<t_int>::Zero(1));

    CHECK(vector_distribution(2, 4)(0) == 0);
    CHECK(vector_distribution(2, 4)(1) == 1);

    CHECK(vector_distribution(4, 4)(0) == 0);
    CHECK(vector_distribution(4, 4)(1) == 1);
    CHECK(vector_distribution(4, 4)(2) == 2);
    CHECK(vector_distribution(4, 4)(3) == 3);

    CHECK(vector_distribution(5, 4).head(2) == Vector<t_int>::Zero(2));
    CHECK(vector_distribution(5, 4)(2) == 1);
    CHECK(vector_distribution(5, 4)(3) == 2);
    CHECK(vector_distribution(5, 4)(4) == 3);

    CHECK(vector_distribution(5, 3).head(2) == Vector<t_int>::Zero(2));
    CHECK(vector_distribution(5, 3).segment(2, 2) == Vector<t_int>::Ones(2));
    CHECK(vector_distribution(5, 3)(4) == 2);

    CHECK(vector_distribution(5, 2).head(3) == Vector<t_int>::Zero(3));
    CHECK(vector_distribution(5, 2).segment(3, 2) == Vector<t_int>::Ones(2));
  }

  SECTION("Local vs non-local") {
    using optimet::mpi::details::local_interactions;
    CHECK(local_interactions(0).rows() == 0);
    CHECK(local_interactions(0).cols() == 0);

    CHECK(local_interactions(1).all());

    CHECK(local_interactions(2).all());

    CHECK(local_interactions(3).block(0, 0, 2, 2).all());
    CHECK(local_interactions(3).block(1, 1, 2, 2).all());
    CHECK(not local_interactions(3)(2, 0));
    CHECK(local_interactions(3) == local_interactions(3).transpose());
    CHECK(local_interactions(3) == local_interactions(3).reverse());
    CHECK(local_interactions(3) == local_interactions(3).reverse().transpose());

    CHECK(local_interactions(4) == local_interactions(4).transpose());
    CHECK(local_interactions(4) == local_interactions(4).reverse());
    CHECK(local_interactions(4) == local_interactions(4).reverse().transpose());
    CHECK(local_interactions(4).block(0, 0, 2, 2).all());
    CHECK(local_interactions(4).block(1, 1, 2, 2).all());
    CHECK(not local_interactions(4).row(0).tail(2).any());
    CHECK(not local_interactions(4).col(3).head(2).any());
  }
}

TEST_CASE("Graph communicators creation") {
  using namespace optimet;
  using mpi::details::vector_distribution;
  using mpi::details::local_graph_edges;
  using mpi::details::non_local_graph_edges;

  auto const nprocs = 3;
  auto const nscatt = 6;
  auto const vecdist = vector_distribution(nscatt, nprocs);
  SECTION("Nothing allowed") {
    CHECK(local_graph_edges(Matrix<bool>::Zero(nscatt, nscatt), vecdist).size() == nprocs);
    for(auto const &connections : local_graph_edges(Matrix<bool>::Zero(nscatt, nscatt), vecdist))
      CHECK(connections.size() == 0);
    CHECK(non_local_graph_edges(Matrix<bool>::Zero(nscatt, nscatt), vecdist).size() == nprocs);
    for(auto const &connections :
        non_local_graph_edges(Matrix<bool>::Zero(nscatt, nscatt), vecdist))
      CHECK(connections.size() == 0);
  }

  SECTION("All allowed") {
    auto const locals = local_graph_edges(Matrix<bool>::Ones(nscatt, nscatt), vecdist);
    CHECK(locals.size() == nprocs);
    for(t_int i(0); i < locals.size(); ++i) {
      CHECK(locals[i].size() == nprocs - 1);
      for(t_uint n(0); n < nprocs; ++n)
        CHECK((locals[i].count(n) == 1) == (n != i));
    }

    auto const nonlocals = non_local_graph_edges(Matrix<bool>::Ones(nscatt, nscatt), vecdist);
    CHECK(nonlocals.size() == nprocs);
    for(t_int i(0); i < locals.size(); ++i) {
      CHECK(nonlocals[i].size() == nprocs - 1);
      for(t_uint n(0); n < nprocs; ++n)
        CHECK((nonlocals[i].count(n) == 1) == (n != i));
    }
  }

  SECTION("Half-half") {
    Matrix<bool> allowed = Matrix<bool>::Zero(nscatt, nscatt);
    allowed.rightCols(nscatt / 2).fill(true);
    auto const locals = local_graph_edges(allowed, vecdist);
    CHECK(locals.size() == nprocs);
    CHECK(locals[0].size() == 0);
    CHECK(locals[1].size() == 2);
    CHECK(locals[1].count(0) == 1);
    CHECK(locals[1].count(2) == 1);
    CHECK(locals[2].size() == 2);
    CHECK(locals[2].count(0) == 1);
    CHECK(locals[2].count(1) == 1);

    auto const nonlocals = non_local_graph_edges(allowed, vecdist);
    CHECK(nonlocals.size() == nprocs);
    CHECK(nonlocals[0].size() == 2);
    CHECK(nonlocals[0].count(1) == 1);
    CHECK(nonlocals[0].count(2) == 1);
    CHECK(nonlocals[1].size() == 1);
    CHECK(nonlocals[1].count(1) == 0);
    CHECK(nonlocals[2].size() == 1);
    CHECK(nonlocals[2].count(0) == 0);
  }
}

TEST_CASE("Input neighbor count") {
  using namespace optimet;
  auto const nprocs = 3;
  auto const nHarmonics = 5;
  std::vector<Scatterer> const scatterers = {
      {{0, 0, 0}, {0, 0}, 0, nHarmonics},     {{0, 0, 0}, {0, 0}, 0, nHarmonics + 1},
      {{0, 0, 0}, {0, 0}, 0, nHarmonics + 2}, {{0, 0, 0}, {0, 0}, 0, nHarmonics + 3},
      {{0, 0, 0}, {0, 0}, 0, nHarmonics + 4}, {{0, 0, 0}, {0, 0}, 0, nHarmonics + 5},
  };
  auto const nfunctions = [&scatterers](int i) {
    return 2 * scatterers[i].nMax * (scatterers[i].nMax + 2);
  };
  auto const nscatt = scatterers.size();
  auto const vecdist = mpi::details::vector_distribution(nscatt, nprocs);

  SECTION("No nonlocal") {
    for(auto rank = 0; rank < nprocs; ++rank) {
      auto const empty = mpi::details::neighborhood_input_counts(Matrix<bool>::Zero(nscatt, nscatt),
                                                                 vecdist, scatterers, rank);
      CHECK(empty.size() == 0);
    }
  }
  SECTION("All local") {
    SECTION("Rank 0") {
      auto const all_but_one = mpi::details::neighborhood_input_counts(
          Matrix<bool>::Ones(nscatt, nscatt), vecdist, scatterers, 0);

      REQUIRE(all_but_one.size() == nprocs - 1);
      CHECK(all_but_one[0] == nfunctions(2) + nfunctions(3));
      CHECK(all_but_one[1] == nfunctions(4) + nfunctions(5));
    }

    SECTION("Rank 1") {
      auto const all_but_one = mpi::details::neighborhood_input_counts(
          Matrix<bool>::Ones(nscatt, nscatt), vecdist, scatterers, 1);

      REQUIRE(all_but_one.size() == nprocs - 1);
      CHECK(all_but_one[0] == nfunctions(0) + nfunctions(1));
      CHECK(all_but_one[1] == nfunctions(4) + nfunctions(5));
    }

    SECTION("Rank 2") {
      auto const all_but_one = mpi::details::neighborhood_input_counts(
          Matrix<bool>::Ones(nscatt, nscatt), vecdist, scatterers, 2);

      REQUIRE(all_but_one.size() == nprocs - 1);
      CHECK(all_but_one[0] == nfunctions(0) + nfunctions(1));
      CHECK(all_but_one[1] == nfunctions(2) + nfunctions(3));
    }
  }

  SECTION("Particular peculiarities") {
    Matrix<bool> non_locals = Matrix<bool>::Zero(nscatt, nscatt);
    non_locals.col(3).fill(true);
    auto counts = mpi::details::neighborhood_input_counts(non_locals, vecdist, scatterers, 0);

    REQUIRE(counts.size() == 1);
    CHECK(counts[0] == nfunctions(3));

    non_locals.col(5).fill(true);
    counts = mpi::details::neighborhood_input_counts(non_locals, vecdist, scatterers, 0);
    REQUIRE(counts.size() == 2);
    CHECK(counts[0] == nfunctions(3));
    CHECK(counts[1] == nfunctions(5));
  }
}
