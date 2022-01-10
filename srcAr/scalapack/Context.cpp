// (C) University College London 2017
// This file is part of Optimet, licensed under the terms of the GNU Public License
//
// Optimet is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Optimet is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Optimet. If not, see <http://www.gnu.org/licenses/>.

#include "scalapack/Blacs.h"
#include "scalapack/Context.h"
#include "scalapack/InitExit.h"
#include <exception>
#include <iostream>

namespace optimet {
namespace scalapack {

void Context::delete_with_finalize(Context::Impl *const impl) {
  if(impl->row >= 0 and impl->col >= 0) {
    if(not finalized())
      OPTIMET_FC_GLOBAL_(blacs_gridexit, BLACS_GRIDEXIT)(&impl->context);
    decrement_ref();
  }
  delete impl;
}
void Context::delete_without_finalize(Context::Impl *const impl) {
  if(impl->row >= 0 and impl->col >= 0)
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
  OPTIMET_FC_GLOBAL_(blacs_pinfo, BLACS_PINFO)(&rank, &size);
  if(size == 0) {
    size = static_cast<int>(rows * cols);
    OPTIMET_FC_GLOBAL_(blacs_setup, BLACS_SETUP)(&rank, &size);
  }

  this_row = -1;
  this_col = 0;
  OPTIMET_FC_GLOBAL_(blacs_get, BLACS_GET)(&this_row, &this_col, &context);
  char order = 'R';
  OPTIMET_FC_GLOBAL_(blacs_gridinit, BLACS_GRIDINIT)(&context, &order, &nrows, &ncols);
  OPTIMET_FC_GLOBAL_(blacs_gridinfo, BLACS_GRIDINFO)
  (&context, &nrows, &ncols, &this_row, &this_col);
  Impl const data{context, static_cast<t_int>(nrows), static_cast<t_int>(ncols),
                  static_cast<t_int>(this_row), static_cast<t_int>(this_col)};
  impl = std::shared_ptr<Impl const>(new Impl(data), &Context::delete_with_finalize);
  if(impl and this_row >= 0 and this_col >= 0)
    increment_ref();
}

Context::Context(Context const &system, Matrix<t_uint> const &gridmap) : impl(nullptr) {
  if(gridmap.size() == 0 or not system.is_valid()) {
    impl =
        std::shared_ptr<Impl const>(new Impl{-1, 0, 0, -1, -1}, &Context::delete_without_finalize);
    return;
  }
  // convert to a type fortran can deal with
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> procs = gridmap.cast<t_int>();
  // gets system context associated with input context
  int context = *system, ldau = procs.rows(), nrows = procs.rows(), ncols = procs.cols(), ten = 10;
  int sys;
  OPTIMET_FC_GLOBAL_(blacs_get, BLACS_GET)(&context, &ten, &sys);
  // Creates the context
  OPTIMET_FC_GLOBAL(blacs_gridmap, BLACS_GRIDMAP)(&sys, procs.data(), &ldau, &nrows, &ncols);
  // Creates the matrix
  int this_row, this_col;
  OPTIMET_FC_GLOBAL_(blacs_gridinfo, BLACS_GRIDINFO)(&sys, &nrows, &ncols, &this_row, &this_col);
  Impl const data{sys, static_cast<t_int>(nrows), static_cast<t_int>(ncols),
                  static_cast<t_int>(this_row), static_cast<t_int>(this_col)};
  impl = std::shared_ptr<Impl const>(new Impl(data), &Context::delete_with_finalize);
  if(impl and this_row >= 0 and this_col >= 0)
    increment_ref();
}

Matrix<t_uint> Context::rank_map() const {
  int context = **this;
  Matrix<t_uint> result(rows(), cols());
  for(int i(0); i < rows(); ++i)
    for(int j(0); j < cols(); ++j)
      result(i, j) =
          static_cast<t_uint>(OPTIMET_FC_GLOBAL_(blacs_pnum, BLACS_PNUM)(&context, &i, &j));
  return result;
}

Context Context::split(t_uint nrows, t_uint ncols) const {
  if(nrows == 0 or ncols == 0)
    throw std::domain_error("nrows and ncols must be > 0");
  if(not is_valid())
    return *this;
  auto const full_map = rank_map();
  Context result(*this);
  for(t_uint i(0), current_r(0); i < nrows; ++i) {
    auto const ni = rows() / nrows + (i < (rows() % nrows) ? 1 : 0);
    for(t_uint j(0), current_c(0); j < ncols; ++j) {
      auto const nj = cols() / ncols + (j < (cols() % ncols) ? 1 : 0);
      auto const intermediate = subcontext(full_map.block(current_r, current_c, ni, nj));
      if(static_cast<t_uint>(row()) >= current_r and static_cast<t_uint>(row()) < current_r + ni and
         static_cast<t_uint>(col()) >= current_c and static_cast<t_uint>(col()) < current_c + nj) {
        assert(result == *this);
        result = intermediate;
      }
      current_c += nj;
    }
    current_r += ni;
  }
  return result;
}

t_uint Context::process_number(t_uint row, t_uint col) const {
  if(not is_valid())
    throw std::runtime_error("Invalid context when checking process number");
  int context = **this;
  int prow = static_cast<int>(row);
  int pcol = static_cast<int>(col);
  return OPTIMET_FC_GLOBAL_(blacs_pnum, BLACS_PNUM)(&context, &prow, &pcol);
}

Index Context::process_coordinates(t_uint pnum) const {
  if(not is_valid())
    throw std::runtime_error("Invalid context when checking process coordinates");
  int context = **this;
  int num = static_cast<int>(pnum);
  int prow, pcol;
  OPTIMET_FC_GLOBAL_(blacs_pcoord, BLACS_PCOORD)(&context, &num, &prow, &pcol);
  return {static_cast<t_uint>(prow), static_cast<t_uint>(pcol)};
}
} /* scalapack  */
} /* optimet  */
