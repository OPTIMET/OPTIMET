#ifndef OPTIMET_SCALAPACK_CONTEXT_H
#define OPTIMET_SCALAPACK_CONTEXT_H

#include "Types.h"
#include <memory>

namespace optimet {
namespace scalapack {

//! A context for a distributed array
class Context {
  //! Holds actual data associated with the context
  struct Impl {
    //! The blacs context
    int context;
    //! The number of rows
    t_uint rows;
    //! The number of cols
    t_uint cols;
    //! Index of this row
    t_uint row;
    //! Index of this row
    t_uint col;
  };

public:
  //! Constructs a context
  Context(t_uint rows, t_uint cols);

  virtual ~Context(){};

  //! Whether this is a valid context for this process
  bool is_valid() const { return static_cast<bool>(impl); }
  //! The number of rows
  decltype(Impl::rows) rows() const { return is_valid() ? impl->rows : 0; }
  //! The number of cols
  decltype(Impl::cols) cols() const { return is_valid() ? impl->cols : 0; }
  //! Total number of processes
  decltype(Impl::cols) size() const { return cols() * rows(); }
  //! Index of this row
  decltype(Impl::row) row() const { return is_valid() ? impl->row : 0; }
  //! Index of this row
  decltype(Impl::col) col() const { return is_valid() ? impl->col : 0; }
  //! Returns the Blacs context in a way blacs undersands
  decltype(Impl::context) operator*() const {
    assert(is_valid());
    return impl->context;
  }

private:
  //! Holds data associated with the context
  std::shared_ptr<Impl const> impl;

  //! Deletes a blacs context
  static void delete_context(Impl *impl);
};

} /* scalapack */
} /* optimet */
#endif /* ifndef OPTIMET_BLACS_CONTEXT */
