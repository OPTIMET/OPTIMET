#ifndef OPTIMET_TYPES_H
#define OPTIMET_TYPES_H

#include <Eigen/Core>
#include <complex>
#include <functional>

#cmakedefine OPTIMET_MPI
#cmakedefine OPTIMET_BELOS
#cmakedefine OPTIMET_CHAR_ARCH
#cmakedefine OPTIMET_LONG_ARCH
#cmakedefine OPTIMET_ULONG_ARCH

#ifdef OPTIMET_BELOS
#include <Tpetra_CrsGraph_decl.hpp>
#include <Tpetra_CrsGraph_decl.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_Map_decl.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#endif

namespace optimet {
//! Root of the type hierarchy for signed integers
typedef int t_int;
//! Root of the type hierarchy for unsigned integers
typedef std::size_t t_uint;
//! Root of the type hierarchy for real numbers
typedef double t_real;
//! Root of the type hierarchy for (real) complex numbers
typedef std::complex<t_real> t_complex;

//! \brief A vector of a given type
//! \details Operates as mathematical vector.
//! \note Scalapack probably expects column-major. Best not to offend it.
template <class T = t_complex> using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>;
//! \brief A matrix of a given type
//! \details Operates as mathematical matrix.
//! \note Scalapack probably expects column-major. Best not to offend it.
template <class T = t_complex>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

#ifdef OPTIMET_BELOS
namespace belos {
//! Map type across OPTIMET's use of Belos
typedef Tpetra::Map<> Map;
//! Special map type to define Tpetra sparse matrices
typedef Tpetra::CrsGraph<> CrsGraph;
//! Simplifies access to a tpetra vector type
template <class SCALAR> using Vector = Tpetra::MultiVector<SCALAR>;
//! Simplifies access to a tpetra vector type
template <class SCALAR> using SparseMatrix = Tpetra::CrsMatrix<SCALAR>;
}
#endif
}
#endif
