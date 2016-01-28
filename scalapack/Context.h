#ifndef OPTIMET_SCALAPACK_CONTEXT_H
#define OPTIMET_SCALAPACK_CONTEXT_H
#include "Types.h"

#ifdef OPTIMET_MPI
#include <memory>
#include "scalapack/InitExit.h"

namespace optimet {
namespace scalapack {

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
  Context(Sizes const &c) : Context(c.rows, c.cols) {}
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

private:
  //! Holds data associated with the context
  std::shared_ptr<Impl const> impl;

  //! Deletes a blacs context
  static void delete_context(Impl *impl);
  //! \brief Deletes the default system blacs context
  //! \details Does not gridexit the context
  static void delete_default_system_context(Impl *impl);
};

} /* scalapack */
} /* optimet */
#endif /* ifdef OPTIMET_MPI */
#endif /* ifndef OPTIMET_BLACS_CONTEXT */
