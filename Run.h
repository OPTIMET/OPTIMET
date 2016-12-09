#ifndef RUN_H_
#define RUN_H_

#include "Geometry.h"
#include "CompoundIterator.h"
#include "Excitation.h"
#include "scalapack/Parameters.h"
#include <memory>

#ifdef OPTIMET_BELOS
#include <Teuchos_ParameterListExceptions.hpp>
#include <Teuchos_ParameterList.hpp>
#endif
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
  std::shared_ptr<Geometry> geometry;     /**< The Geometry of the case. */
  std::shared_ptr<optimet::Excitation> excitation; /**< The Excitation of the case. */
  //! Parameters needed to setup parallel computations
  optimet::scalapack::Parameters parallel_params;
#ifdef OPTIMET_BELOS
  Teuchos::RCP<Teuchos::ParameterList> belos_params;
#endif

  int nMax; /**< The maximum value of the n iterator. */

  // This bit will be moved to the case or where it is appropiate
  int projection;
  double params[9];
  int outputType;  /**< Output type required: 0 -> Field, 1 -> Cross Sections, 2
                      -> Scattering Coefficients. */
  bool singleMode; /**< Output only one mode (harmonic) in the field profile. */
  CompoundIterator singleModeIndex; /**< Index of single mode to output in the
                                       field profile. */
  bool dominantAuto;                /**< Get dominant mode automatically. */
  int singleComponent; /**< Get only one or both components: 0 -> Both, 1 - >
                          TE, 2 - > TM. */

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
  Run() : geometry(new Geometry) {};

  /**
   * Default destructor for the Case class.
   */
  virtual ~Run() {};
};

#endif /* RUN_H_ */
