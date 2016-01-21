#include <exception>
#include "BlacsContext.h"
#include "BlacsExit.h"
#include "Blacs.h"

namespace optimet {

void BlacsContext::delete_context(BlacsContext::Impl *const impl) {
  OPTIMET_FC_GLOBAL(blacs_gridexit, BLACS_GRIDEXIT)(&impl->context);
  decrement_blacs_ref();
  delete impl;
}

BlacsContext::BlacsContext(t_uint rows, t_uint cols) : impl(nullptr) {

  int nrows = static_cast<t_uint>(rows);
  int ncols = static_cast<t_uint>(cols);
  int this_row = -1, this_col = -1;
  int context;
  OPTIMET_FC_GLOBAL(sl_init, SL_INIT)(&context, &nrows, &ncols);
  if(nrows >= 0 and ncols >= 0 and this_row >= 0 and this_col >= 0) {
    Impl const data{context, static_cast<t_uint>(nrows), static_cast<t_uint>(ncols),
                    static_cast<t_uint>(this_row), static_cast<t_uint>(this_col)};
    impl = std::shared_ptr<Impl const>(new Impl(data), &BlacsContext::delete_context);
    if(impl)
      increment_blacs_ref();
  }
}

} /* optimet  */
