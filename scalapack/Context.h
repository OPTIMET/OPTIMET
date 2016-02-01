#ifndef OPTIMET_SCALAPACK_CONTEXT_H
#define OPTIMET_SCALAPACK_CONTEXT_H
#include "Types.h"

#include <memory>
#include "scalapack/Parameters.h"
#include "scalapack/InitExit.h"

namespace optimet {
namespace scalapack {

#ifdef OPTIMET_MPI
//! A context for a distributed array
class Context {
  //! Holds actual data associated with the context
  struct Impl {
    //! The blacs context
    int context;
    //! The number of rows
    t_int rows;
    //! The number of cols
    t_int cols;
    //! Index of this row
    t_int row;
    //! Index of this row
    t_int col;
  };

public:
  //! Constructs a context
  Context(t_uint rows, t_uint cols);
  //! Constructs using the default system context
  Context() : Context(1, global_size()) {};
  //! Constructs a context
  Context(scalapack::Sizes const &c) : Context(c.rows, c.cols) {}
  //! Constructs a context from a gridmap
  Context(Context const &system, Matrix<t_uint> const & gridmap);
  //! Constructs a context from a gridmap
  Context(Matrix<t_uint> const & gridmap) : Context(Context(), gridmap) {};

  virtual ~Context(){};

  //! Whether this is a valid context for this process
  bool is_valid() const { return impl->row >= 0 and impl->col >= 0; }
  //! The number of rows
  decltype(Impl::rows) rows() const { return is_valid() ? impl->rows : 0; }
  //! The number of cols
  decltype(Impl::cols) cols() const { return is_valid() ? impl->cols : 0; }
  //! Total number of processes
  t_uint size() const { return is_valid() ? static_cast<t_uint>(cols() * rows()) : 0; }
  //! Index of this row
  decltype(Impl::row) row() const { return is_valid() ? impl->row : 0; }
  //! Index of this row
  decltype(Impl::col) col() const { return is_valid() ? impl->col : 0; }
  //! Returns the Blacs context in a way blacs undersands
  decltype(Impl::context) operator*() const { return impl->context; }

  //! Check contexts are the same
  bool operator==(Context const &c) const { return **this == *c; }
  //! Check contexts are the same
  bool operator!=(Context const &c) const { return not operator==(c); }

  //! Map from grid to rank
  Matrix<t_uint> rank_map() const;
  //! \brief Split into nrows * ncols smaller grids
  //! \details Each axis of N points is split into n chunks of N / n + (rank < N % n ? 1: 0) points.
  //! Doing this to both axis results in grids. If the  n > N for both axis, then the output grid is
  //! equivalent to the input grid.
  Context split(t_uint nrows, t_uint ncols) const;
  //! Creates a subcontext
  Context subcontext(Matrix<t_uint> map) const { return Context(*this, map); }

  //! \brief Creates the largest squarest context
  //! \see optimet::scalapack::squarest_largest_grid()
  static Context Squarest(t_uint n, t_real skew) { return Context(squarest_largest_grid(n, skew)); }
  //! \brief Creates the largest squarest context
  //! \see optimet::scalapack::squarest_largest_grid()
  static Context Squarest(t_uint n) { return Context(squarest_largest_grid(n)); }
  //! \brief Creates the largest squarest context
  //! \see optimet::scalapack::squarest_largest_grid()
  static Context Squarest() { return Context(squarest_largest_grid(global_size())); }

private:
  //! Holds data associated with the context
  std::shared_ptr<Impl const> impl;

  //! Deletes a blacs context
  static void delete_with_finalize(Impl *impl);
  //! \brief Deletes the default system blacs context
  //! \details Does not gridexit the context
  static void delete_without_finalize(Impl *impl);
};
#else
//! A fake scalapack context when compiled without scalapack
class Context {
public:
  //! Constructs a context
  Context(t_uint, t_uint) {}
  Context(t_uint) {}
  Context() {}

  //! Whether this is a valid context for this process
  bool is_valid() const { return true; }
  //! The number of rows
  t_int rows() const { return 0; }
  //! The number of cols
  t_int cols() const { return 0; }
  //! Total number of processes
  t_uint size() const { return 1; }
  //! Index of this row
  t_int row() const { return 0; }
  //! Index of this row
  t_int col() const { return 0; }
  //! Returns the Blacs context in a way blacs undersands
  t_int operator*() const {
    throw std::runtime_error("No context since compiled without scalapack");
  }

  //! \brief Creates the largest squarest context
  //! \see optimet::scalapack::squarest_largest_grid()
  static Context Squarest(t_uint, t_real) { return Context(); }
  //! \brief Creates the largest squarest context
  //! \see optimet::scalapack::squarest_largest_grid()
  static Context Squarest(t_uint) { return Context(); }
  //! \brief Creates the largest squarest context
  //! \see optimet::scalapack::squarest_largest_grid()
  static Context Squarest() { return Context(); }
};
#endif

} /* scalapack */
} /* optimet */
#endif /* ifdef OPTIMET_MPI */
#endif /* ifndef OPTIMET_BLACS_CONTEXT */
