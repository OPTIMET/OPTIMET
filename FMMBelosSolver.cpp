#include "FMMBelosSolver.h"
#include "PreconditionedMatrix.h"
#include "scalapack/LinearSystemSolver.h"
#include <Kokkos_View.hpp>
#include <Teuchos_RCP.hpp>

namespace optimet {
namespace solver {
namespace {
//! Type for belos to figure out how to apply FMM
typedef std::reference_wrapper<mpi::FastMatrixMultiply const> FMMOperator;
//! Type for belos to figure out how to apply FMM
typedef Tpetra::MultiVector<t_complex> TpetraVector;
}
}
}

namespace Belos {
template <>
//! Partial specialization of OperatorTraits for Tpetra objects.
class OperatorTraits<optimet::t_complex, optimet::solver::TpetraVector,
                     optimet::solver::FMMOperator> {
public:
  static void Apply(optimet::solver::FMMOperator const &Op, optimet::solver::TpetraVector const &X,
                    optimet::solver::TpetraVector &Y, const ETrans trans = NOTRANS) {
    using namespace optimet;
    if(X.getLocalLength() != Y.getLocalLength())
      throw std::runtime_error("Local lengths of X and Y are different");
    if(X.getLocalLength() == 0)
      return;
    auto const input =
        Eigen::Map<Vector<t_complex> const>(X.getData(0).getRawPtr(), X.getLocalLength());
    auto output = Vector<t_complex>::Map(Y.getDataNonConst(0).getRawPtr(), Y.getLocalLength());
    if(trans == Belos::TRANS)
      output = Op.get().transpose(input);
    else if(trans == Belos::CONJTRANS)
      output = Op.get().adjoint(input);
    else
      output = Op.get() * input;
  }

  static bool HasApplyTranspose(optimet::solver::FMMOperator const &) { return true; }
};
}

namespace optimet {
namespace solver {
namespace {
Teuchos::RCP<const Teuchos::Comm<int>> teuchos_communicator(mpi::Communicator const &comm);
Teuchos::RCP<Tpetra::MultiVector<t_complex>>
tpetra_vector(t_uint nglobals, Vector<t_complex> &x,
              Teuchos::RCP<const Teuchos::Comm<int>> const &comm);
Teuchos::RCP<Tpetra::MultiVector<t_complex>>
tpetra_vector(t_uint nglobals, Vector<t_complex> const &x,
              Teuchos::RCP<const Teuchos::Comm<int>> const &comm);
}

void FMMBelos::update() {
  if(geometry and incWave and communicator_.is_valid()) {
    auto const diags = subdiagonals == std::numeric_limits<t_int>::max() ?
                           std::max<int>(1, geometry->objects.size() / 2 - 2) :
                           subdiagonals;
    fmm_ = std::make_shared<mpi::FastMatrixMultiply>(geometry->bground, incWave->wavenumber(),
                                                     geometry->objects, diags, communicator_);
    auto const distribution =
        mpi::details::vector_distribution(geometry->objects.size(), communicator_.size());
    auto const first = std::find(distribution.data(), distribution.data() + distribution.size(),
                                 communicator_.rank()) -
                       distribution.data();
    auto const last =
        std::find_if(distribution.data() + first, distribution.data() + distribution.size(),
                     [this](t_int value) { return value != communicator_.rank(); }) -
        distribution.data();
    Q = source_vector(geometry->objects.begin() + first, geometry->objects.begin() + last, incWave);
  } else {
    fmm_ = nullptr;
    Q = Vector<t_complex>::Zero(0);
  }
}

void FMMBelos::solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_,
                     mpi::Communicator const &comm) const {
  if(*comm != *communicator_)
    throw std::runtime_error(
        "Solver must be used with the same communicator it is constructed with");
  return solve(X_sca_, X_int_);
}

void FMMBelos::solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_) const {
  auto const distribution =
      mpi::details::vector_distribution(geometry->objects.size(), communicator_.size());

  // used to create vector
  auto const nglobals = geometry->scatterer_size();
  auto const nlocals = Q.size();
  X_sca_.resize(nlocals);
  X_sca_.fill(0);

  auto const tcom = teuchos_communicator(communicator_);
  auto const x = tpetra_vector(nglobals, X_sca_, tcom);
  auto const b = tpetra_vector(nglobals, Q, tcom);
  auto Aptr = Teuchos::rcp(new FMMOperator(*fmm_));

  typedef Belos::LinearProblem<t_complex, TpetraVector, FMMOperator> BelosLinearProblem;
  auto const problem = rcp(new BelosLinearProblem(Aptr, x, b));
  // Tell the solver what problem you want to solve.
  if(not problem->setProblem())
    throw std::runtime_error("Could not setup up Belos problem");

  typedef Belos::SolverFactory<t_complex, TpetraVector, FMMOperator> BelosSolverFactory;
  auto solver = BelosSolverFactory().create(belos_params_->get("Solver", "GMRES"), belos_params_);
  solver->setProblem(problem);

  // Print out parameters for given verbosity
  if(belos_params_->get<int>("Verbosity", 0) & Belos::MsgType::FinalSummary) {
    auto const out = belos_params_->get<Teuchos::RCP<std::ostream>>(
        "Output Stream", Teuchos::rcp(&std::cout, false));
    if(communicator_.rank() == 0)
      solver->getCurrentParameters()->print(*out);
  }

  if(solver->solve() != Belos::Converged)
    throw std::runtime_error("Belos optimizer did not converge");

  X_sca_ = Eigen::Map<Vector<t_complex> const>(x->getData(0).getRawPtr(), x->getLocalLength());
  X_sca_ = communicator_.all_gather(X_sca_);
  X_sca_ = AbstractSolver::convertIndirect(X_sca_);
  X_int_ = AbstractSolver::solveInternal(X_sca_);
}

namespace {
// Construct teuchos mpi wrapper, making sure it owns it.
Teuchos::RCP<const Teuchos::Comm<int>> teuchos_communicator(mpi::Communicator const &comm) {
  MPI_Comm raw_comm;
  MPI_Comm_dup(*comm, &raw_comm);
  Teuchos::RCP<Teuchos::OpaqueWrapper<MPI_Comm>> opaque_comm =
      Teuchos::opaqueWrapper(raw_comm, MPI_Comm_free);
  return Teuchos::rcp(new Teuchos::MpiComm<int>(opaque_comm));
}

Teuchos::RCP<Tpetra::MultiVector<t_complex>>
tpetra_vector(t_uint nglobals, Vector<t_complex> &x,
              Teuchos::RCP<const Teuchos::Comm<int>> const &comm) {
  Teuchos::ArrayView<t_complex> const array_view(x.data(), x.size());
  auto const array_map = Teuchos::rcp(new Tpetra::Map<>(nglobals, x.size(), 0, comm));
  return Teuchos::rcp(new Tpetra::MultiVector<t_complex>(array_map, array_view, x.size(), 1));
}

Teuchos::RCP<Tpetra::MultiVector<t_complex>>
tpetra_vector(t_uint nglobals, Vector<t_complex> const &x,
              Teuchos::RCP<const Teuchos::Comm<int>> const &comm) {
  Teuchos::ArrayView<t_complex const> const array_view(x.data(), x.size());
  auto const array_map = Teuchos::rcp(new Tpetra::Map<>(nglobals, x.size(), 0, comm));
  return Teuchos::rcp(new Tpetra::MultiVector<t_complex>(array_map, array_view, x.size(), 1));
}
}
}
}
