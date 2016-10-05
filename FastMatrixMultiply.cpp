#include "FastMatrixMultiply.h"
#include "RotationCoaxialDecomposition.h"
#include "Types.h"
#include <Eigen/Dense>

namespace optimet {
std::vector<t_uint> FastMatrixMultiply::compute_indices(std::vector<Scatterer> const &scatterers) {
  std::vector<t_uint> result{0u};
  for(auto const &scatterer : scatterers)
    result.push_back(result.back() + 2 * scatterer.nMax + 1);
  return result;
}

std::vector<Rotation>
FastMatrixMultiply::compute_rotations(std::vector<Scatterer> const &scatterers, Range incident,
                                      Range translate) {
  if(incident.first < 0 or translate.first < 0)
    throw std::out_of_range("Start of range must be positive");
  if(incident.first > incident.second)
    throw std::out_of_range("End of range before beginning of range");
  if(translate.first > translate.second)
    throw std::out_of_range("End of range before beginning of range");
  if(incident.second < scatterers.size())
    throw std::out_of_range("End of range past last item");
  if(translate.second < scatterers.size())
    throw std::out_of_range("End of range past last item");

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
      auto const a2 = (x0 - out_begin->vR.toEigenCartesian()).eval();
      auto const a1 = z.cross(a2).normalized().eval();
      auto const a0 = (a1.stableNorm() > 1e-12 ? a1 : Eigen::Matrix<t_real, 3, 1>(0, 1, 0))
                          .cross(a2)
                          .normalized()
                          .eval();
      auto const theta = a2(2);
      auto const phi = std::atan2(a2(1), a0(2));
      result.emplace_back(theta, phi, chi, std::max(in_begin->nMax, out_begin->nMax));
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
    auto const rows = in_begin->nMax * (in_begin->nMax + 2);
    auto const v = in_begin->getTLocal(wavenumber * constant::pi, background);
    assert(result.size() >= i + rows);
    assert(v.size() >= 2 * total_size);

    result.segment(i, rows) = v.head(rows);
    result.segment(i + total_size, rows) = v.tail(rows);
    // incrementing here avoids problems with variable incrementation order
    i += rows;
  }
  return result;
}

void FastMatrixMultiply::operator()(Vector<t_complex> const &in, Vector<t_complex> &out) const {
  auto const offset =
      global_indices_[incident_range_.second] - global_indices_[incident_range_.first];
  auto const local_size = 2 * offset;
  if(in.size() != local_size)
    throw std::runtime_error("Incorrect incident vector size");
  out.resize(local_size);
  out.fill(0);
  // The matrix form makes it easier to focus on one scattering particle at a time
  // Φ in first column and Ψ in second
  Eigen::Map<Matrix<t_complex>> out_matrix(out.data(), offset, 2);

  // Adds identity component (left-hand-side of Eq 106 in Gumerov, Duraiswami 2007)
  if(translate_range_.first < incident_range_.second and
     translate_range_.second > incident_range_.first) {
    auto const start = std::max(translate_range_.first, incident_range_.first);
    auto const end = std::min(translate_range_.second, incident_range_.second);
    auto const n = end - start;
    out.segment(start - translate_range_.first, n) = in.segment(start - incident_range_.first, n);
  }

  // Adds right-hand-side of Eq 106 in Gumerov, Duraiswami 2007
  // first applies Mie coefficients to effective incident field
  // It transforms incident from R to S basis (S is radiating component)
  auto scattered = apply_mie_coefficients(in);
  // The matrix form makes it easier to focus on one scattering particle at a time
  // Φ in first column and Ψ in second
  Eigen::Map<Matrix<t_complex> const> const matrix_scattered(scattered.data(), offset, 2);

  // create a work matrix with appropriate size
  auto in_begin = scatterers_.cbegin() + incident_range_.first;
  auto const in_end = scatterers_.cbegin() + incident_range_.second;
  auto const cmp_nmax = [](Scatterer const &a, Scatterer const &b) { return a.nMax < b.nMax; };
  auto const max_nMax =
      2 * std::max(std::max_element(in_begin, in_end, cmp_nmax)->nMax,
                   std::max_element(scatterers_.begin() + translate_range_.first,
                                    scatterers_.begin() + translate_range_.second, cmp_nmax)
                       ->nMax);
  Matrix<t_complex> work0(2 * max_nMax + 1, 1), work1(2 * max_nMax + 1, 1);

  // Adds left-hand-side of Eq 106 in Gumerov, Duraiswami 2007
  // This is done one at a time for each scatterer -> translated location pair
  // e.g. for each scatterer and particle on which the EM field impinges.
  auto rotation = rotations_.cbegin();
  for(t_int i(0); in_begin != in_end; ++in_begin) {
    auto const in_rows = 2 * in_begin->nMax + 1;
    auto const incident = matrix_scattered.topRows(i + in_rows).bottomRows(in_rows);
    auto out_begin = scatterers_.cbegin() + translate_range_.first;
    auto const out_end = scatterers_.cbegin() + translate_range_.second;
    for(t_int j(0); out_begin != out_end; ++out_begin, ++rotation) {
      auto const out_rows = 2 * out_begin->nMax + 1;
      auto const rows = std::max(out_rows, in_rows);
      // Apply rotation from external basis to basis with z = q -> q'
      auto rotated = work0.topRows(in_rows);
      rotation->transpose(incident, rotated);
      // if in_rows < out_rows fill missing coefficients with zeros
      if(in_rows < rows)
        rotated.bottomRows(rows - in_rows).fill(0);

      // Apply S|R coaxial translation from q to q'
      auto translated = work1.topRows(rows);
      apply_translation(translated, work0.topRows(rows));

      // Apply field coaxial translation
      auto const tz =
          (out_begin->vR.toEigenCartesian() - in_begin->vR.toEigenCartesian()).stableNorm();
      auto cotranslated = work1.topRows(rows);
      rotation_coaxial_decomposition(wavenumber_, tz, translated, cotranslated);

      // Apply rotation from basis with z = q -> q' to external basis
      auto rotated_back = out_matrix.topRows(j + out_rows).bottomRows(out_rows);
      if(out_rows < rows) {
        // if in_rows > out_rows, then we need only some of the coefficients we have computed
        auto rotated_large = work0.topRows(rows);
        rotation->conjugate(cotranslated, rotated_large);
        rotated_back = rotated_large.topLeftCorner(out_rows, out_rows);
      } else
        rotation->conjugate(cotranslated, rotated_back);
      j += out_rows;
    }
    i += in_rows;
  }
}

Vector<t_complex> FastMatrixMultiply::apply_mie_coefficients(Vector<t_complex> const &in) const {
  Vector<t_complex> result = Vector<t_complex>::Zero(in.size());
  assert(2 * (global_indices_[incident_range_.second] - global_indices_[incident_range_.first]));
  auto const size = in.size() / 2;
  // Φ coeffs
  result.head(size) = mie_coefficients_.array() * in.head(size).array();
  // Ψ coeffs
  result.tail(size) = mie_coefficients_.array() * in.tail(size).array();
  return result;
}

} // optimet namespace
