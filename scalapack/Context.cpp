#include <exception>
#include "scalapack/Context.h"
#include "scalapack/InitExit.h"
#include "scalapack/Blacs.h"

namespace optimet {
namespace scalapack {

void Context::delete_context(Context::Impl *const impl) {
  if(not finalized())
    OPTIMET_FC_GLOBAL(blacs_gridexit, BLACS_GRIDEXIT)(&impl->context);
  decrement_ref();
  delete impl;
}

Context::Context(t_uint rows, t_uint cols) : impl(nullptr) {
  if(rows * cols == 0)
    throw std::runtime_error("Context of size 0");
  int nrows = static_cast<t_uint>(rows);
  int ncols = static_cast<t_uint>(cols);
  int this_row = -1, this_col = -1;
  int context;
  int rank, size = 0;
  OPTIMET_FC_GLOBAL(blacs_pinfo, BLACS_PINFO)(&rank, &size);
  if(size == 0) {
    size = static_cast<int>(rows * cols);
    OPTIMET_FC_GLOBAL(blacs_setup, BLACS_SETUP)(&rank, &size);
  }

  this_row = -1;
  this_col = 0;
  OPTIMET_FC_GLOBAL(blacs_get, BLACS_GET)(&this_row, &this_col, &context);
  char order = Matrix<>::IsRowMajor ? 'R' : 'C';
  OPTIMET_FC_GLOBAL(blacs_gridinit, BLACS_GRIDINIT)(&context, &order, &nrows, &ncols);
  OPTIMET_FC_GLOBAL(blacs_gridinfo, BLACS_GRIDINFO)(&context, &nrows, &ncols, &this_row, &this_col);
  if(nrows >= 0 and ncols >= 0 and this_row >= 0 and this_col >= 0) {
    Impl const data{context, static_cast<t_uint>(nrows), static_cast<t_uint>(ncols),
                    static_cast<t_uint>(this_row), static_cast<t_uint>(this_col)};
    impl = std::shared_ptr<Impl const>(new Impl(data), &Context::delete_context);
    if(impl)
      increment_ref();
  }
}

} /* scalapack  */
} /* optimet  */
