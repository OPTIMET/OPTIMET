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

#include "Bessel.h"
#include "Coefficients.h"
#include "FastMatrixMultiply.h"
#include "RotationCoaxialDecomposition.h"
#include "Types.h"
#include <Eigen/Dense>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace optimet {
namespace {
void range_sanity(t_int Nscatterers, t_int rows, t_int cols) {
  if(Nscatterers < 0)
    throw std::out_of_range("Negative number of scatterers");
  if(rows != Nscatterers)
    throw std::out_of_range("Size of couplings and scatterers do not match");
  if(cols != Nscatterers)
    throw std::out_of_range("Size of couplings and scatterers do not match");
}
}

std::vector<std::pair<t_uint, t_uint>>
FastMatrixMultiply::compute_indices(Matrix<bool> const &couplings) {
  Indices result;
  for(t_int i(0); i < couplings.rows(); ++i)
    for(t_int j(0); j < couplings.rows(); ++j)
      if(couplings(i, j))
        result.emplace_back(i, j);
  return result;
}

std::vector<t_uint> FastMatrixMultiply::compute_offsets(std::vector<Scatterer> const &scatterers,
                                                        Vector<bool> const &couplings) {
  std::vector<t_uint> result(couplings.size() + 1);
  result[0] = 0;
  for(t_uint i(0); i < couplings.size(); ++i) {
    if(couplings(i))
      result[i + 1] = 2 * scatterers[i].nMax * (scatterers[i].nMax + 2) + result[i];
    else
      result[i + 1] = result[i];
  }
  return result;
}

std::vector<Rotation>
FastMatrixMultiply::compute_rotations(std::vector<Scatterer> const &scatterers,
                                      Matrix<bool> const &couplings) {
  range_sanity(scatterers.size(), couplings.rows(), couplings.cols());

  std::vector<Rotation> result;
  result.reserve(couplings.count());

  auto const chi = constant::pi;
  Eigen::Matrix<t_real, 3, 1> const z(0, 0, 1);
  for(t_int i(0); i < couplings.rows(); ++i)
    for(t_int j(0); j < couplings.rows(); ++j) {
      if(not couplings(i, j))
        continue;
      if(i == j) {
        result.emplace_back(0, 0, 0, 1);
        continue;
      }
      auto const &in_scatt = scatterers[j];
      auto const &out_scatt = scatterers[i];
      auto const a2 =
          (out_scatt.vR.toEigenCartesian() - in_scatt.vR.toEigenCartesian()).normalized().eval();
      auto const theta = std::acos(a2(2));
      auto const phi = std::atan2(a2(1), a2(0));
      result.emplace_back(theta, phi, chi, std::max(in_scatt.nMax, out_scatt.nMax));
      assert((result.back().basis_rotation().adjoint() * a2).isApprox(Vector<t_real>::Unit(3, 2)));
      assert((result.back().basis_rotation() * Vector<t_real>::Unit(3, 2)).isApprox(a2));
    }
  return result;
}

std::vector<CachedCoAxialRecurrence::Functor>
FastMatrixMultiply::compute_coaxial_translations(t_complex wavenumber,
                                                 std::vector<Scatterer> const &scatterers,
                                                 Matrix<bool> const &couplings) {
  range_sanity(scatterers.size(), couplings.rows(), couplings.cols());
  std::vector<CachedCoAxialRecurrence::Functor> result;
  result.reserve(couplings.count());
  for(t_int i(0); i < couplings.rows(); ++i)
    for(t_int j(0); j < couplings.rows(); ++j) {
      if(not couplings(i, j))
        continue;
      if(i == j) {
        result.push_back(CachedCoAxialRecurrence(0, 10, false).functor(1));
        continue;
      }
      auto const &in_scatt = scatterers[j];
      auto const &out_scatt = scatterers[i];
      auto const Orad = in_scatt.vR.toEigenCartesian();
      auto const Ononrad = out_scatt.vR.toEigenCartesian();
      CachedCoAxialRecurrence tca((Orad - Ononrad).stableNorm(), wavenumber, false);
      result.push_back(tca.functor(std::max(in_scatt.nMax, out_scatt.nMax) + nplus));
    }
  return result;
}

Vector<t_complex>
FastMatrixMultiply::compute_mie_coefficients(ElectroMagnetic const &background, t_real wavenumber,
                                             std::vector<Scatterer> const &scatterers,
                                             Matrix<bool> const &couplings) {
  auto const outs = couplings.colwise().any().eval();
  // First figure out total size
  t_uint total_size = 0;
  for(t_int i(0); i < outs.size(); ++i)
    if(outs(i))
      total_size += nfunctions(scatterers[i].nMax);
  Vector<t_complex> result(2 * total_size);

  // then actually compute vector
  for(t_uint i(0), j(0); i < outs.size(); ++i) {
    if(not outs(i))
      continue;
    auto const v = scatterers[i].getTLocal(wavenumber * constant::c, background);
    assert(result.size() >= j + v.size());
    result.segment(j, v.size()) = v;
    // incrementing here avoids problems with variable incrementation order
    j += v.size();
  }
  return result;
}

Eigen::Array<t_real, Eigen::Dynamic, 2>
FastMatrixMultiply::compute_normalization(std::vector<Scatterer> const &scatterers) {
  if(scatterers.size() == 0)
    return Eigen::Array<t_real, Eigen::Dynamic, 2>::Zero(0, 2);
  auto const cmp = [](Scatterer const &a, Scatterer const &b) { return a.nMax < b.nMax; };
  auto const nMax = std::max_element(scatterers.begin(), scatterers.end(), cmp)->nMax + nplus;
  Eigen::Array<t_real, Eigen::Dynamic, 2> result(nfunctions(nMax), 2);
  for(t_int n(1), i(0); n <= nMax; i += 2 * n + 1, ++n) {
    result.col(0).segment(i, 2 * n + 1).fill(1e0 / std::sqrt((n * (n + 1)) / 2));
    result.col(1).segment(i, 2 * n + 1).fill(1e0 / std::sqrt((n * (n + 1)) / 2));
    for(t_int m(-n); m <= n; ++m) {
      if(m >= 0 or m % 2 == 0)
        result(i + n + m, 0) *= -1;
      else
        result(i + n + m, 1) *= -1;
    }
  }
  return result;
}

void FastMatrixMultiply::operator()(Vector<t_complex> const &in, Vector<t_complex> &out) const {
  if(in.size() != cols())
    throw std::runtime_error("Incorrect incident vector size");
  out.resize(rows());
  if(out.size() == 0)
    return;
  out.fill(0);

  // Adds identity component (left-hand-side of Eq 106 in Gumerov, Duraiswami 2007)
  // It is added only in those instances which compute the diagonal couplings/self-interactions.
  for(Indices::size_type i(0); i < indices_.size(); ++i)
    if(indices_[i].first == indices_[i].second) {
      auto const n = 2 * nfunctions(incident_nmax(i));
      out.segment(translate_offset(i), n) = in.segment(incident_offset(i), n);
    }

  // Adds right-hand-side of Eq 106 in Gumerov, Duraiswami 2007
  translation(mie_coefficients_.array() * in.array(), out);
}

void FastMatrixMultiply::transpose(Vector<t_complex> const &in, Vector<t_complex> &out) const {
  if(in.size() != rows())
    throw std::runtime_error("Incorrect incident vector size");
  out.resize(cols());
  if(cols() == 0)
    return;
  out.fill(0);

  // Adds right-hand-side of Eq 106 in Gumerov, Duraiswami 2007
  translation_transpose(in, out);
  // Adds mie coefficient last when transposing
  out.array() *= mie_coefficients_.array();

  // Adds identity component (left-hand-side of Eq 106 in Gumerov, Duraiswami 2007)
  for(Indices::size_type i(0); i < indices_.size(); ++i)
    if(is_self_interaction(i)) {
      auto const n = 2 * nfunctions(incident_nmax(i));
      out.segment(incident_offset(i), n) += in.segment(translate_offset(i), n);
    }
}

t_int FastMatrixMultiply::max_nmax() const {
  t_uint nmax = 0;
  for(auto const &indices : indices_)
    nmax = std::max<t_uint>(scatterers_[indices.first].nMax,
                            std::max<t_uint>(scatterers_[indices.second].nMax, nmax));
  return nmax + nplus;
}

void FastMatrixMultiply::translation(Vector<t_complex> const &input, Vector<t_complex> &out) const {
  typedef Eigen::Matrix<t_complex, Eigen::Dynamic, 2> Matrixified;

  // create a work matrix with appropriate size
  // It should have nplus (degree) more harmonics than the maximum object + the n = 0 term (1
  // element)
  Eigen::Matrix<t_complex, Eigen::Dynamic, 4> work(nfunctions(max_nmax()) + 1, 4);

  // Adds left-hand-side of Eq 106 in Gumerov, Duraiswami 2007
  // This is done one at a time for each scatterer -> translated location pair
  // e.g. for each scatterer and particle on which the EM field impinges.
  for(Indices::size_type i(0); i < indices_.size(); ++i) {
    // no self-interaction
    if(is_self_interaction(i))
      continue;
    auto const in_rows = nfunctions(incident_nmax(i));
    Eigen::Map<const Matrixified> const incident(input.data() + incident_offset(i), in_rows, 2);
    auto const out_rows = nfunctions(translate_nmax(i));
    Eigen::Map<Matrixified> translate(out.data() + translate_offset(i), out_rows, 2);
    remove_translation(incident, translate, work, i);
  }
}

void FastMatrixMultiply::translation_transpose(Vector<t_complex> const &input,
                                               Vector<t_complex> &out) const {
  typedef Eigen::Matrix<t_complex, Eigen::Dynamic, 2> Matrixified;

  // create a work matrix with appropriate size
  Eigen::Matrix<t_complex, Eigen::Dynamic, 4> work(nfunctions(max_nmax()) + 1, 4);

  // Adds left-hand-side of Eq 106 in Gumerov, Duraiswami 2007
  // This is done one at a time for each scatterer -> translated location pair
  // e.g. for each scatterer and particle on which the EM field impinges.
  for(Indices::size_type i(0); i < indices_.size(); ++i) {
    if(is_self_interaction(i))
      continue;
    auto const in_rows = nfunctions(incident_nmax(i));
    Eigen::Map<Matrixified> const incident(out.data() + incident_offset(i), in_rows, 2);
    auto const out_rows = nfunctions(translate_nmax(i));
    Eigen::Map<const Matrixified> translate(input.data() + translate_offset(i), out_rows, 2);
    remove_translation_transpose(translate, incident, work, i);
  }
}

Vector<t_complex> FastMatrixMultiply::operator()(Vector<t_complex> const &in) const {
  Vector<t_complex> result(rows());
  operator()(in, result);
  return result;
}

Vector<t_complex> FastMatrixMultiply::transpose(Vector<t_complex> const &in) const {
  Vector<t_complex> result(cols());
  transpose(in, result);
  return result;
}
} // optimet namespace
