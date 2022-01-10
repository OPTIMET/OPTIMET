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

#ifndef OPTIMET_RUN_H_
#define OPTIMET_RUN_H_

#include "CompoundIterator.h"
#include "Excitation.h"
#include "Geometry.h"
#include "Types.h"
#include "mpi/Communicator.h"
#include "scalapack/Context.h"
#include "scalapack/Parameters.h"
#include <array>
#include <memory>

#ifdef OPTIMET_BELOS
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterListExceptions.hpp>
#endif
namespace optimet {
/**
 * The Run class implements a single instance of a run.
 * Run has several components:
 *  1. The Geometry - needed on all nodes so must be copied.
 *  2. The Excitation - specific excitation data; iteration of wavelength,
 *            etc. will be done using this.
 *  3. The Solver - solves the problem.
 *  4. The Result - the final data once request is processed.
 */
class Run {
public:
  //! The Geometry of the case
  std::shared_ptr<Geometry> geometry;
  //! The Excitation of the case
  std::shared_ptr<Excitation> excitation;
  //! Parameters needed to setup parallel computations
  scalapack::Parameters parallel_params;
#ifdef OPTIMET_BELOS
  Teuchos::RCP<Teuchos::ParameterList> belos_params;
#endif

  //!  The maximum value of the n iterator
  t_int nMax;
  t_int nMaxS; // second harmonic number of spherical harmonics

  //! This bit will be moved to the case or where it is appropiate
  t_int projection;
  std::array<t_real, 9> params;
  //! Output type required: 0 -> Field, 1 -> Cross Sections, 2 -> Scattering Coefficients
  t_int outputType;
  //! Output only one mode (harmonic) in the field profile
  bool singleMode;
  //! Index of single mode to output in the field profile
  CompoundIterator singleModeIndex;
  //! Get dominant mode automatically
  bool dominantAuto;
  //! Get only one or both components: 0 -> Both, 1 - > TE, 2 - > TM
  t_int singleComponent;
  scalapack::Context context;
  mpi::Communicator communicator;

  //! Wether to run with FMM or concrete matrix
  bool do_fmm;
  //! Number of subdiagonals when setting up fmm local vs non-local mpi distribution
  t_int fmm_subdiagonals;

  /**
   * Params:
   *    -> for Field see OutputGrid
   *    -> For cross section only 3 used: params[0] - initial lambda, params[1]
   * - final lambda, params[2] - number of steps
   */

  /**
   * Default constructor for the Case class.
   * Does NOT initialize the instance.
   */
  Run() : geometry(new Geometry), context(scalapack::Context::Squarest()){};

  /**
   * Default destructor for the Case class.
   */
  virtual ~Run(){};
};
}
#endif /* RUN_H_ */
