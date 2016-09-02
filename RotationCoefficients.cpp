#include "HarmonicsIterator.h"
#include "RotationCoefficients.h"
#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#ifdef OPTIMET_BELOS
#include <Teuchos_ArrayViewDecl.hpp>
#include <Teuchos_DefaultComm.hpp>
#endif

namespace optimet {

RotationCoefficients::Real RotationCoefficients::a(t_uint n, t_int m) {
  t_uint const absm(std::abs(m));
  if(n < absm)
    return static_cast<Real>(0);
  return std::sqrt(static_cast<Real>((n + 1 + absm) * (n + 1 - absm)) /
                   static_cast<Real>((2 * n + 1) * (2 * n + 3)));
}

RotationCoefficients::Real RotationCoefficients::b(t_uint n, t_int m) {
  if(static_cast<t_uint>(std::abs(m)) > n)
    return static_cast<Real>(0);
  return (m >= 0 ? 1 : -1) * std::sqrt(static_cast<Real>((n - m - 1) * (n - m)) /
                                       static_cast<Real>((2 * n - 1) * (2 * n + 1)));
}

RotationCoefficients::Coefficients
RotationCoefficients::factors(t_uint n, t_int m, t_int mu) const {
  auto const factor = std::exp(Complex(0, chi_)) / b(n + 1, m - 1);
  auto const half_factor = static_cast<Real>(0.5) * factor;
  auto const c0 =
      half_factor * b(n + 1, -mu - 1) * std::exp(Complex(0, phi_)) * (1 - std::cos(theta_));
  auto const c1 =
      -half_factor * b(n + 1, mu - 1) * std::exp(Complex(0, -phi_)) * (1 + std::cos(theta_));
  auto const c2 = -factor * a(n, mu) * std::sin(theta_);
  return Coefficients(c0, c1, c2);
}

RotationCoefficients::Complex RotationCoefficients::with_caching(t_uint n, t_int m, t_int mu) {
  if(static_cast<t_uint>(std::abs(m)) > n or static_cast<t_uint>(std::abs(mu)) > n)
    return static_cast<Real>(0);
  if(m < 0)
    return std::conj(with_caching(n, -m, -mu));

  auto const prior = cache.find(std::make_tuple(n, m, mu));
  if(prior != cache.end())
    return prior->second;

  auto const value = recursion(n, m, mu);
  cache[std::make_tuple(n, m, mu)] = value;
  return value;
}

RotationCoefficients::Complex RotationCoefficients::recursion(t_uint n, t_int m, t_int mu) {
  if(static_cast<t_uint>(std::abs(m)) > n or static_cast<t_uint>(std::abs(mu)) > n)
    return static_cast<Real>(0);
  if(m < 0)
    return std::conj(with_caching(n, -m, -mu));

  if(m == 0)
    return initial(n, mu);

  auto const factors = this->factors(n, m, mu);
  return factors(0) * with_caching(n + 1, m - 1, mu + 1) +
         factors(1) * with_caching(n + 1, m - 1, mu - 1) +
         factors(2) * with_caching(n + 1, m - 1, mu);
}

#ifdef OPTIMET_BELOS
namespace belos {
Teuchos::RCP<Map const> map(t_uint nmax, t_uint nobjects, mpi::Communicator const &comm) {
  typedef int ordinal;
  auto const n = nobjects / comm.size() + (nobjects % comm.size() > comm.rank() ? 1 : 0);
  auto const per_object = 2 * (HarmonicsIterator::max_flat(nmax) - 1);
  auto const global = nobjects * per_object;
  auto const local = n * per_object;

#ifdef OPTIMET_MPI
  // Construct teuchos mpi wrapper, making sure it owns it.
  MPI_Comm raw_comm;
  MPI_Comm_dup(*comm, &raw_comm);
  Teuchos::RCP<Teuchos::OpaqueWrapper<MPI_Comm>> opaque_comm =
      Teuchos::opaqueWrapper(raw_comm, MPI_Comm_free);
  Teuchos::RCP<const Teuchos::Comm<ordinal>> tcomm =
      Teuchos::rcp(new Teuchos::MpiComm<ordinal>(opaque_comm));
#else
  Teuchos::RCP<const Teuchos::Comm<ordinal>> tcomm = Teuchos::DefaultComm<ordinal>::getComm();
#endif

  return Teuchos::rcp(new Tpetra::Map<>(global, local, 0, tcomm));
}

CrsGraph crs_rotation_graph(t_uint nmax, t_uint nobjects, mpi::Communicator const &comm) {
  typedef int ordinal;
  auto const row_map = map(nmax, nobjects, comm);
  auto const col_map = map(nmax, nobjects, comm);
  // See equation 5.15 in Nail A. Gumerov, Ramani Duraiswami, 2003.
  std::vector<ordinal> columns;
  auto const max_entry_per_row = 2 * nmax + 1;
  auto const offset = nmax * (nmax + 2);
  CrsGraph result(row_map, col_map, max_entry_per_row);
  for(t_uint object(0), i(0); object < nobjects; ++object)
    for(t_uint n(1); n < nmax; ++n) {
      columns.resize(2 * n + 1);
      // entries to rotate harmonics in Φ
      std::iota(columns.begin(), columns.end(), i);
      for(t_uint c(0), j(i); c < 2 * n + 1; ++c)
        result.insertLocalIndices(j, columns);
      // entries to rotate harmonics in Ψ
      std::iota(columns.begin(), columns.end(), i + offset);
      for(t_uint c(0), j(i + offset); c < 2 * n + 1; ++c, ++i, ++j)
        result.insertLocalIndices(j, columns);
    }
  result.fillComplete();
  return result;
}
}
#endif
}
