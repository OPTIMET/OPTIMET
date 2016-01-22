#ifndef OPTIMET_BLACS_CONTEXT
#define OPTIMET_BLACS_CONTEXT

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

  virtual ~Context() {};

  //! Whether this is a valid context for this process
  bool is_valid() const { return static_cast<bool>(impl); }
  //! The number of rows
  decltype(Impl::rows) rows() const { return impl->rows; }
  //! The number of cols
  decltype(Impl::cols) cols() const { return impl->cols; }
  //! Index of this row
  decltype(Impl::row) row() const { return impl->row; }
  //! Index of this row
  decltype(Impl::col) col() const { return impl->col; }
  //! Returns the Blacs context in a way blacs undersands
  decltype(Impl::context) operator*() const { return impl->context; }

private:
  //! Holds data associated with the context
  std::shared_ptr<Impl const> impl;

  //! Deletes a blacs context
  static void delete_context(Impl *impl);
};

} /* scalapack */
} /* optimet */
#endif /* ifndef OPTIMET_BLACS_CONTEXT */
