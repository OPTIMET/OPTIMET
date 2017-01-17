#include "Coupling.h"
#include "PreconditionedMatrix.h"
#include "Types.h"
#include "scalapack/BroadcastToOutOfContext.h"

namespace optimet {
#ifdef OPTIMET_SCALAPACK
Vector<t_complex> distributed_source_vector(Vector<t_complex> const &input,
                                            scalapack::Context const &context,
                                            scalapack::Sizes const &blocks) {
  if(not context.is_valid())
    return Vector<t_complex>::Zero(0);
  auto const serial = context.serial();
  t_uint const n(input.size());
  Eigen::Map<Matrix<t_complex> const> const input_map(input.data(), serial.is_valid() ? n : 0,
                                                      serial.is_valid() ? 1 : 0);
  scalapack::Matrix<t_complex const *> const serial_vector(input_map, serial, {n, 1}, {n, 1});
  scalapack::Matrix<t_complex> result(context, {n, 1}, blocks);
  serial_vector.transfer_to(context, result);
  if(result.local().cols() == 0)
    return Vector<t_complex>::Zero(0);
  return result.local();
}

Vector<t_complex> gather_all_source_vector(t_uint n, Vector<t_complex> const &input,
                                           scalapack::Context const &context,
                                           scalapack::Sizes const &blocks) {
  if(not context.is_valid())
    return Vector<t_complex>::Zero(0);
  auto const serial = context.serial();
  auto const parallel_vector = map_cmatrix(input, context, {n, 1}, blocks);
  scalapack::Matrix<t_complex> result(context.serial(), {n, 1}, {n, 1});
  parallel_vector.transfer_to(context, result);
  return context.broadcast(result.local(), 0, 0);
}
Vector<t_complex> gather_all_source_vector(scalapack::Matrix<t_complex> const &matrix) {
  scalapack::Matrix<t_complex> result(matrix.context().serial(), {matrix.rows(), 1},
                                      {matrix.rows(), 1});
  matrix.transfer_to(matrix.context(), result);
  auto const result_vector = matrix.context().broadcast(result.local(), 0, 0);
  if(result_vector.size() == 0)
    return Vector<t_complex>::Zero(0);
  return result_vector;
}
#endif

Matrix<t_complex>
preconditioned_scattering_matrix(std::vector<Scatterer>::const_iterator const &first,
                                 std::vector<Scatterer>::const_iterator const &end_first,
                                 std::vector<Scatterer>::const_iterator const &second,
                                 std::vector<Scatterer>::const_iterator const &end_second,
                                 ElectroMagnetic const &bground,
                                 std::shared_ptr<Excitation const> incWave) {
  auto const nMax = first->nMax;
  auto const n = nMax * (nMax + 2);
  if(first == end_first or second == end_second)
    return Matrix<t_complex>::Zero(2 * n * (end_first - first), 2 * n * (end_second - second));

  Matrix<t_complex> result(2 * n * (end_first - first), 2 * n * (end_second - second));
  size_t y(0);
  for(auto iterj(second); iterj != end_second; ++iterj, y += 2 * n) {
    Vector<t_complex> const factor = -iterj->getTLocal(incWave->omega(), bground);
    size_t x(0);
    for(auto iteri(first); iteri != end_first; ++iteri, x += 2 * n) {
      if(iteri == iterj) {
        result.block(x, y, 2 * n, 2 * n) = Matrix<t_complex>::Identity(2 * n, 2 * n);
      } else {
        Coupling const AB(iteri->vR - iterj->vR, incWave->waveK, nMax);
        result.block(x, y, n, n) = AB.diagonal.transpose();
        result.block(x + n, y + n, n, n) = AB.diagonal.transpose();
        result.block(x, y + n, n, n) = AB.offdiagonal.transpose();
        result.block(x + n, y, n, n) = AB.offdiagonal.transpose();
        result.block(x, y, 2 * n, 2 * n).array().transpose().colwise() *= factor.array();
      }
    }
  }
  return result;
}

Matrix<t_complex> preconditioned_scattering_matrix(std::vector<Scatterer> const &objects,
                                                   ElectroMagnetic const &bground,
                                                   std::shared_ptr<Excitation const> incWave) {
  return preconditioned_scattering_matrix(objects.begin(), objects.end(), objects.begin(),
                                          objects.end(), bground, incWave);
}

Matrix<t_complex> preconditioned_scattering_matrix(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave) {
  if(geometry.objects.size() == 0)
    return Matrix<t_complex>(0, 0);
  // Check nMax is same accross all objects
  auto const nMax = geometry.objects.front().nMax;
  for(auto const &scatterer : geometry.objects)
    if(scatterer.nMax != nMax)
      throw std::runtime_error("All objects must have same number of harmonics");
  return preconditioned_scattering_matrix(geometry.objects, geometry.bground, incWave);
}

#ifdef OPTIMET_SCALAPACK
Matrix<t_complex> preconditioned_scattering_matrix(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave,
                                                   scalapack::Context const &context,
                                                   scalapack::Sizes const &blocks) {
  // construct an n by 1 context
  auto const nobj = geometry.objects.size();
  if(nobj == 0)
    return Matrix<t_complex>::Zero(0, 0);
  auto rank_map = context.rank_map();
  rank_map.resize(1, context.size());
  auto const linear_context = context.subcontext(rank_map.leftCols(std::min(context.size(), nobj)));

  auto const nMax = geometry.objects.front().nMax;
  auto const remainder = linear_context.is_valid() ? nobj % linear_context.size() : 0;
  auto const nloc = linear_context.is_valid() ? nobj / linear_context.size() : 0;
  t_uint const n = nMax * (nMax + 2);
  scalapack::Sizes const non_cyclic{linear_context.is_valid() ? nobj * n * 2 : 1,
                                    linear_context.is_valid() ? nloc * n * 2 : 1};
  scalapack::Matrix<t_complex> linear_matrix(linear_context, {nobj * n * 2, nobj * n * 2},
                                             non_cyclic);
  if(linear_context.is_valid()) {
    assert(geometry.objects.end() > geometry.objects.begin() + nloc * linear_context.col());
    assert(geometry.objects.end() >= geometry.objects.begin() + nloc * (1 + linear_context.col()));
    assert(nloc * 2 * n > 0);
    linear_matrix.local().leftCols(nloc * 2 * n) = preconditioned_scattering_matrix(
        geometry.objects.begin(), geometry.objects.end(),
        geometry.objects.begin() + nloc * linear_context.col(),
        geometry.objects.begin() + nloc * (1 + linear_context.col()), geometry.bground, incWave);
  }

  if(remainder > 0 and linear_context.is_valid()) {
    auto const remainder_context = linear_context.subcontext(rank_map.leftCols(remainder));
    scalapack::Matrix<t_complex> remainder_matrix(
        remainder_context, {nobj * n * 2, remainder * n * 2}, {nobj * n * 2, 2 * n});
    if(remainder_context.is_valid())
      remainder_matrix.local() = preconditioned_scattering_matrix(
          geometry.objects.begin(), geometry.objects.end(),
          geometry.objects.begin() + nloc * linear_context.cols() + remainder_context.col(),
          geometry.objects.begin() + nloc * linear_context.cols() + remainder_context.col() + 1,
          geometry.bground, incWave);

    scalapack::Matrix<t_complex> transfered(linear_context, {nobj * n * 2, remainder * n * 2},
                                            {nobj * n * 2, nloc * n * 2});
    remainder_matrix.transfer_to(linear_context, transfered);
    if(transfered.local().cols() > 0)
      linear_matrix.local().rightCols(transfered.local().cols()) = transfered.local();
  }

  scalapack::Matrix<t_complex> distributed_matrix(context, linear_matrix.sizes(), blocks);

  linear_matrix.transfer_to(context, distributed_matrix);
  return distributed_matrix.local();
}
#else
Matrix<t_complex> preconditioned_scattering_matrix(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave,
                                                   scalapack::Context const &,
                                                   scalapack::Sizes const &) {
  return preconditioned_scattering_matrix(geometry, incWave);
}
#endif

Vector<t_complex> source_vector(std::vector<Scatterer>::const_iterator first,
                                std::vector<Scatterer>::const_iterator const &last,
                                std::shared_ptr<Excitation const> incWave) {
  if(first == last)
    return Vector<t_complex>::Zero(0);
  auto const nMax = first->nMax;
  auto const flatMax = nMax * (nMax + 2);
  Vector<t_complex> result(2 * flatMax * (last - first));
  for(size_t i(0); first != last; ++first, i += 2 * flatMax)
    incWave->getIncLocal(first->vR, result.data() + i, nMax);
  return result;
}

Vector<t_complex>
source_vector(std::vector<Scatterer> const &objects, std::shared_ptr<Excitation const> incWave) {
  return source_vector(objects.begin(), objects.end(), incWave);
}

Vector<t_complex>
source_vector(Geometry const &geometry, std::shared_ptr<Excitation const> incWave) {
  if(geometry.objects.size() == 0)
    return Vector<t_complex>(0, 0);
  // Check nMax is same accross all objects
  auto const nMax = geometry.objects.front().nMax;
  for(auto const &scatterer : geometry.objects)
    if(scatterer.nMax != nMax)
      throw std::runtime_error("All objects must have same number of harmonics");
  return source_vector(geometry.objects, incWave);
}

Vector<t_complex> local_source_vector(Geometry const &geometry,
                                      std::shared_ptr<Excitation const> incWave,
                                      Vector<t_complex> const &input_coeffs) {
  if(geometry.objects.size() == 0)
    return Vector<t_complex>(0, 0);
  // Check nMax is same accross all objects
  auto const nMax = geometry.objects.front().nMax;
  for(auto const &scatterer : geometry.objects)
    if(scatterer.nMax != nMax)
      throw std::runtime_error("All objects must have same number of harmonics");

  Geometry copy_geometry(geometry);
  copy_geometry.setSourcesSingle(incWave, input_coeffs.data(), nMax);
  // Get the IncLocal matrices
  // These correspond directly to the Beta*a in Stout2002 Eq. 10 as
  // they are already translated.
  // we are in the SH case -> get the local sources from the geometry
  auto const flatMax = nMax * (nMax + 2);
  Vector<t_complex> result(flatMax * 2 * copy_geometry.objects.size());
  for(size_t i = 0; i < copy_geometry.objects.size(); i++)
    copy_geometry.getSourceLocal(i, incWave, nMax, result.data() + i * 2 * flatMax);
  return result;
}
}
