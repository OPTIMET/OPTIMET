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

#include "catch.hpp"
#include <iostream>

#include "Aliases.h"
#include "Geometry.h"
#include "Scatterer.h"
#include "ScalapackSolver.h"
#include "MatrixBelosSolver.h"
#include "Tools.h"
#include "Types.h"
#include "constants.h"
#include <BelosTypes.hpp>

using namespace optimet;

TEST_CASE("Scalapack vs Belos") {
  auto geometry = std::make_shared<Geometry>();
  // spherical coords, ε, μ, radius, nmax
  auto const nSpheres = 10;
  auto const nHarmonics = 10;
  for(t_uint i(0); i < 10; ++i) {
    t_real const ii(i);
    geometry->pushObject({{ii * 1.5 * 2e-6, 0, 0},
                         {0.45e0 + 0.1 * ii, 1.1e0},
                         (0.5 + 0.01 * ii) * 2e-6,
                         nHarmonics});
  }

  // Create excitation
  auto const wavelength = 14960e-9;
  Spherical<t_real> const vKinc{2 * consPi / wavelength, 90 * consPi / 180.0, 90 * consPi / 180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  auto excitation =
      std::make_shared<Excitation>(0, Tools::toProjection(vKinc, Eaux), vKinc, nHarmonics);
  excitation->populate();
  geometry->update(excitation);

  optimet::solver::Scalapack const scalapack_solver(geometry, excitation);
  optimet::Result scalapack(geometry, excitation);
  scalapack_solver.solve(scalapack.scatter_coef, scalapack.internal_coef);

  optimet::solver::MatrixBelos solver(geometry, excitation);
  solver.belos_parameters()->set<int>("Num Blocks", solver.scattering_size());
  solver.belos_parameters()->set("Maximum Iterations", 4000);
  solver.belos_parameters()->set("Convergence Tolerance", 1.0e-14);
  // solver.belos_parameters()->set(
  //     "Verbosity", Belos::MsgType::Warnings + Belos::MsgType::Debug + Belos::MsgType::Errors +
  //                      Belos::MsgType::StatusTestDetails + Belos::MsgType::IterationDetails +
  //                      Belos::MsgType::FinalSummary);
  // solver.belos_parameters()->set("Output Frequency", 1);

  // known to fail: "CGPOLY", "FLEXIBLE GMRES", "RECYCLING CG", "RCG", "PCPG", "MINRES", "LSQR",
  // "SEED CG",
  auto const names = {
      "BICGSTAB",
      "BLOCK GMRES",
      "CG",                 // "PSEUDO BLOCK CG",
      "GMRES",              // "PSEUDO BLOCK GMRES",
      "GMRESPOLY",          // "HYBRID BLOCK GMRES", "SEED GMRES"
      "PSEUDO BLOCK TFQMR", // "PSEUDO BLOCK TRANSPOSE-FREE QMR",
      "GCRODR",             // "RECYCLING GMRES",
      "STOCHASTIC CG",      // "PSEUDO BLOCK STOCHASTIC CG"
      "TFQMR",              // "TRANSPOSE-FREE QMR",
      "BLOCK CG",
      "BLOCK GMRES",
      "FIXED POINT",
  };
  for(auto const name : {"GMRES", "CG"}) {
    SECTION(name) {
      solver.belos_parameters()->set("Solver", name);
      optimet::Result belos(geometry, excitation);
      solver.solve(belos.scatter_coef, belos.internal_coef);

      auto const scatter_tol = 1e-6 * std::max(1., scalapack.scatter_coef.array().abs().maxCoeff());
      CHECK(belos.scatter_coef.isApprox(scalapack.scatter_coef, scatter_tol));
      auto const internal_tol =
          1e-6 * std::max(1., scalapack.internal_coef.array().abs().maxCoeff());
      CHECK(belos.internal_coef.isApprox(belos.internal_coef, internal_tol));
    }
  }
}
