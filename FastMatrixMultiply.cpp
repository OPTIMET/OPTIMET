#include "Bessel.h"
#include "Coefficients.h"
#include "FastMatrixMultiply.h"
#include "RotationCoaxialDecomposition.h"
#include "Types.h"
#include <Eigen/Dense>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <iostream>

namespace optimet {
namespace {
void range_sanity(t_int Nscatterers, FastMatrixMultiply::Range incident,
                  FastMatrixMultiply::Range translate) {
  if(incident.first < 0 or translate.first < 0)
    throw std::out_of_range("Start of range must be positive");
  if(incident.first > incident.second)
    throw std::out_of_range("End of range before beginning of range");
  if(translate.first > translate.second)
    throw std::out_of_range("End of range before beginning of range");
  if(incident.second < Nscatterers)
    throw std::out_of_range("End of range past last item");
  if(translate.second < Nscatterers)
    throw std::out_of_range("End of range past last item");
}
}
std::vector<t_uint> FastMatrixMultiply::compute_indices(std::vector<Scatterer> const &scatterers) {
  std::vector<t_uint> result{0u};
  for(auto const &scatterer : scatterers)
    result.push_back(result.back() + 2 * scatterer.nMax * (scatterer.nMax + 2));
  return result;
}

std::vector<Rotation>
FastMatrixMultiply::compute_rotations(std::vector<Scatterer> const &scatterers, Range incident,
                                      Range translate) {
  range_sanity(scatterers.size(), incident, translate);
  std::vector<Rotation> result;
  result.reserve((incident.second - incident.first) * (translate.second - translate.first));
  auto in_begin = scatterers.cbegin() + incident.first;
  auto const in_end = scatterers.cbegin() + incident.second;
  auto const chi = constant::pi;
  Eigen::Matrix<t_real, 3, 1> const z(0, 0, 1);
  for(; in_begin != in_end; ++in_begin) {
    auto out_begin = scatterers.cbegin() + translate.first;
    auto const out_end = scatterers.cbegin() + translate.second;
    auto const x0 = in_begin->vR.toEigenCartesian();
    for(; out_begin != out_end; ++out_begin) {
      if(out_begin == in_begin)
        continue;
      auto const a2 = (out_begin->vR.toEigenCartesian() - x0).normalized().eval();
      auto const theta = std::acos(a2(2));
      auto const phi = std::atan2(a2(1), a2(0));
      result.emplace_back(theta, phi, chi, std::max(in_begin->nMax, out_begin->nMax));
    }
  }
  return result;
}

std::vector<CachedCoAxialRecurrence::Functor>
FastMatrixMultiply::compute_coaxial_translations(t_complex wavenumber,
                                                 std::vector<Scatterer> const &scatterers,
                                                 Range incident, Range translate) {
  range_sanity(scatterers.size(), incident, translate);
  std::vector<CachedCoAxialRecurrence::Functor> result;
  result.reserve((incident.second - incident.first) * (translate.second - translate.first));
  auto in_begin = scatterers.cbegin() + incident.first;
  auto const in_end = scatterers.cbegin() + incident.second;
  for(; in_begin != in_end; ++in_begin) {
    auto out_begin = scatterers.cbegin() + translate.first;
    auto const out_end = scatterers.cbegin() + translate.second;
    auto const Orad = in_begin->vR.toEigenCartesian();
    for(; out_begin != out_end; ++out_begin) {
      if(out_begin == in_begin)
        continue;
      auto const Ononrad = out_begin->vR.toEigenCartesian();
      CachedCoAxialRecurrence tca((Orad - Ononrad).stableNorm(), wavenumber, false);
      result.push_back(tca.functor(std::max(in_begin->nMax, out_begin->nMax) + nplus));
    }
  }
  return result;
}

Vector<t_complex>
FastMatrixMultiply::compute_mie_coefficients(ElectroMagnetic const &background, t_real wavenumber,
                                             std::vector<Scatterer> const &scatterers,
                                             Range incident) {
  auto in_begin = scatterers.cbegin() + incident.first;
  auto const in_end = scatterers.cbegin() + incident.second;
  // First figure out total size
  t_uint total_size = 0;
  for(; in_begin != in_end; ++in_begin)
    total_size += in_begin->nMax * (in_begin->nMax + 2);
  Vector<t_complex> result(2 * total_size);

  // then actually compute vector
  in_begin = scatterers.cbegin() + incident.first;
  for(t_uint i(0); in_begin != in_end; ++in_begin) {
    auto const v = in_begin->getTLocal(wavenumber * constant::c, background);
    assert(result.size() >= i + v.size());
    result.segment(i, v.size()) = v;
    // incrementing here avoids problems with variable incrementation order
    i += v.size();
  }
  return result;
}

void FastMatrixMultiply::operator()(Vector<t_complex> const &in, Vector<t_complex> &out) const {
  auto const in_offset =
      global_indices_[incident_range_.second] - global_indices_[incident_range_.first];
  auto const out_offset =
      global_indices_[translate_range_.second] - global_indices_[translate_range_.first];
  if(in.size() != in_offset)
    throw std::runtime_error("Incorrect incident vector size");
  out.resize(out_offset);
  out.fill(0);

  // Adds identity component (left-hand-side of Eq 106 in Gumerov, Duraiswami 2007)
  if(std::max(translate_range_.first, incident_range_.first) <
     std::min(translate_range_.second, incident_range_.second)) {
    auto const start = global_indices_[std::max(translate_range_.first, incident_range_.first)];
    auto const end = global_indices_[std::min(translate_range_.second, incident_range_.second)];
    auto const n = end - start;
    out.segment(start - translate_range_.first, n) = in.segment(start - incident_range_.first, n);
  }

  // Adds right-hand-side of Eq 106 in Gumerov, Duraiswami 2007
  translation(mie_coefficients_.array() * in.array(), out);
}

void FastMatrixMultiply::translation(Vector<t_complex> const &input, Vector<t_complex> &out) const {

  auto const in_size =
      global_indices_[incident_range_.second] - global_indices_[incident_range_.first];
  auto const out_size =
      global_indices_[translate_range_.second] - global_indices_[translate_range_.first];
  if(input.size() != in_size)
    throw std::runtime_error("Incorrect incident vector size");
  if(out.size() != out_size)
    throw std::runtime_error("Incorrect outgoing vector size");

  typedef Eigen::Matrix<t_complex, Eigen::Dynamic, 2> Matrixified;

  // create a work matrix with appropriate size
  auto in_begin = scatterers_.cbegin() + incident_range_.first;
  auto const in_end = scatterers_.cbegin() + incident_range_.second;
  auto const cmp_nmax = [](Scatterer const &a, Scatterer const &b) { return a.nMax < b.nMax; };
  auto const max_nMax =
      std::max(std::max_element(in_begin, in_end, cmp_nmax)->nMax,
               std::max_element(scatterers_.begin() + translate_range_.first,
                                scatterers_.begin() + translate_range_.second, cmp_nmax)
                   ->nMax) +
      nplus;
  Eigen::Matrix<t_complex, Eigen::Dynamic, 4> work(nfunctions(max_nMax), 4);

  // Adds left-hand-side of Eq 106 in Gumerov, Duraiswami 2007
  // This is done one at a time for each scatterer -> translated location pair
  // e.g. for each scatterer and particle on which the EM field impinges.
  auto i_rotation = rotations_.cbegin();
  auto i_translation = coaxial_translations_.cbegin();
  for(t_int i(0); in_begin != in_end; ++in_begin) {
    auto const in_rows = nfunctions(in_begin->nMax);
    Eigen::Map<const Matrixified> const incident(input.data() + i, in_rows, 2);
    auto out_begin = scatterers_.cbegin() + translate_range_.first;
    auto const out_end = scatterers_.cbegin() + translate_range_.second;
    for(t_int j(0); out_begin != out_end; ++out_begin) {
      auto const out_rows = nfunctions(out_begin->nMax);
      // no self-interaction
      if(in_begin == out_begin) {
        j += 2 * out_rows;
        continue;
      }
      // add field from j to i
      Eigen::Map<Matrixified> outgoing(out.data() + j, out_rows, 2);
      remove_translation(incident, *in_begin, *out_begin, *i_rotation, *i_translation, wavenumber_,
                         outgoing, work);
      j += 2 * out_rows;
      ++i_rotation;
      ++i_translation;
    }
    i += 2 * in_rows;
  }
}

Vector<t_complex> FastMatrixMultiply::operator()(Vector<t_complex> const &in) const {
  auto const offset =
      global_indices_[translate_range_.second] - global_indices_[translate_range_.first];
  Vector<t_complex> result(offset * 2);
  operator()(in, result);
  return result;
}

Matrix<t_complex> field_coefficients(t_real wavenumber, Vector<t_complex> const &coeffs,
                                     std::vector<Scatterer> const &objects) {
  Matrix<t_complex> results(coeffs.size() / 2, 6);
  t_int row(0);
  for(auto const &object : objects) {
    auto const N = object.nMax;
    auto const size = N * (N + 2);
    auto qcoeffs = coeffs.segment(2 * row, 2 * size);
    auto qout = results.middleRows(row, size);
    field_coefficients(wavenumber, qcoeffs, qout);
    row += size;
  }
  return results;
}

Fields::Result Fields::operator()(Position const &position) const {
  t_int row(0);
  for(auto const &object : objects) {
    auto const N = 2 * object.nMax * (object.nMax + 2);
    auto const r = position - object.vR.toEigenCartesian();
    if(r.stableNorm() <= object.radius)
      return single_sphere(r, object.nMax, inside.middleRows(row, N));
    row += N;
  }
  return compute_outside(position);
}

Fields::Result Fields::compute_outside(Position const &position) const {
  Result result = Result::Zero();
  t_int row(0);
  for(auto const &object : objects) {
    auto const N = 2 * object.nMax * (object.nMax + 2);
    auto const r = position - object.vR.toEigenCartesian();
    result += single_sphere(r, object.nMax, outside.middleRows(row, N));
    row += N;
  }
  return result;
};

} // optimet namespace
