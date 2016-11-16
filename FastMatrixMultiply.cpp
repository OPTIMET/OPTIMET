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
      if(out_begin == in_begin) {
        result.emplace_back(0, 0, 0, 1);
        continue;
      }
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
      if(out_begin == in_begin) {
        result.push_back(CachedCoAxialRecurrence(0, 10, false).functor(1));
        continue;
      }
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

Eigen::Array<t_real, Eigen::Dynamic, 2>
FastMatrixMultiply::compute_normalization(std::vector<Scatterer> const &scatterers) {
  auto const cmp = [](Scatterer const &a, Scatterer const &b) { return a.nMax < b.nMax; };
  auto const nMax = std::max_element(scatterers.begin(), scatterers.end(), cmp)->nMax + nplus;
  Eigen::Array<t_real, Eigen::Dynamic, 2> result(nfunctions(nMax), 2);
  for(t_int n(1), i(0); n <= nMax; ++n) {
    result.col(0).segment(i, 2 * n + 1).fill(1e0 / std::sqrt((n * (n + 1)) / 2));
    result.col(1).segment(i, 2 * n + 1).fill(-1e0 / std::sqrt((n * (n + 1)) / 2));
    for(t_int m(-n); m <= n; ++m, ++i)
      if(m > 0 and m % 2 == 1)
        result.row(i) *= -1;
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

void FastMatrixMultiply::transpose(Vector<t_complex> const &in, Vector<t_complex> &out) const {
  auto const out_offset =
      global_indices_[incident_range_.second] - global_indices_[incident_range_.first];
  auto const in_offset =
      global_indices_[translate_range_.second] - global_indices_[translate_range_.first];
  if(in.size() != in_offset)
    throw std::runtime_error("Incorrect incident vector size");
  out.resize(out_offset);
  out.fill(0);

  // Adds right-hand-side of Eq 106 in Gumerov, Duraiswami 2007
  translation_transpose(in, out);
  // Adds mie coefficient last when transposing
  out.array() *= mie_coefficients_.array();

  // Adds identity component (left-hand-side of Eq 106 in Gumerov, Duraiswami 2007)
  if(std::max(translate_range_.first, incident_range_.first) <
     std::min(translate_range_.second, incident_range_.second)) {
    auto const start = global_indices_[std::max(translate_range_.first, incident_range_.first)];
    auto const end = global_indices_[std::min(translate_range_.second, incident_range_.second)];
    auto const n = end - start;
    out.segment(start - translate_range_.first, n) += in.segment(start - incident_range_.first, n);
  }
}

t_int FastMatrixMultiply::max_nmax() const {
  auto const cmp_nmax = [](Scatterer const &a, Scatterer const &b) { return a.nMax < b.nMax; };
  auto const in_max = std::max_element(scatterers_.cbegin() + incident_range_.first,
                                       scatterers_.cbegin() + incident_range_.second, cmp_nmax)
                          ->nMax;
  auto const out_max = std::max_element(scatterers_.cbegin() + translate_range_.first,
                                        scatterers_.cbegin() + translate_range_.second, cmp_nmax)
                           ->nMax;
  return std::max(in_max, out_max) + nplus;
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
  auto const ntranslates = translate_range_.second - translate_range_.first;
  auto const nincidents = incident_range_.second - incident_range_.first;
  for(t_int i(0), iparticle(0); iparticle < nincidents; ++iparticle) {
    auto const in_rows = nfunctions(incident_nmax(iparticle));
    Eigen::Map<const Matrixified> const incident(input.data() + i, in_rows, 2);
    for(t_int j(0), jparticle(0); jparticle < ntranslates; ++jparticle) {
      auto const out_rows = nfunctions(translate_nmax(jparticle));
      // no self-interaction
      if(incident_range_.first + iparticle != translate_range_.first + jparticle) {
        // add field from j to i
        Eigen::Map<Matrixified> outgoing(out.data() + j, out_rows, 2);
        remove_translation(incident, outgoing, work, iparticle, jparticle);
      }
      j += 2 * out_rows;
    }
    i += 2 * in_rows;
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
  auto const ntranslates = translate_range_.second - translate_range_.first;
  auto const nincidents = incident_range_.second - incident_range_.first;
  for(t_int i(0), iparticle(0); iparticle < nincidents; ++iparticle) {
    auto const in_rows = nfunctions(incident_nmax(iparticle));
    Eigen::Map<Matrixified> const incident(out.data() + i, in_rows, 2);
    for(t_int j(0), jparticle(0); jparticle < ntranslates; ++jparticle) {
      auto const out_rows = nfunctions(translate_nmax(jparticle));
      // no self-interaction
      if(incident_range_.first + iparticle != translate_range_.first + jparticle) {
        // add field from j to i
        Eigen::Map<const Matrixified> translate(input.data() + j, out_rows, 2);
        remove_translation_transpose(translate, incident, work, iparticle, jparticle);
      }
      j += 2 * out_rows;
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

Vector<t_complex> FastMatrixMultiply::transpose(Vector<t_complex> const &in) const {
  auto const offset =
      global_indices_[incident_range_.second] - global_indices_[incident_range_.first];
  Vector<t_complex> result(offset * 2);
  transpose(in, result);
  return result;
}

} // optimet namespace
