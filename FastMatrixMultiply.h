#ifndef FAST_MATRIX_MULTIPLY_H
#define FAST_MATRIX_MULTIPLY_H

#include "Bessel.h"
#include "CoAxialTranslationCoefficients.h"
#include "RotationCoaxialDecomposition.h"
#include "RotationCoefficients.h"
#include "Scatterer.h"
#include "Types.h"
#include <utility>
#include <vector>

namespace optimet {
//! \brief Fast multiplication of effective incident field by transfer matrix
//! \details The multiplication occurs for a set of scattering particles receiving the EM
//! radiation
//! to a set of locations where the fields are checked. Each incident effective field is
//! composed of coefficients for both Φ and Ψ potentials (in R basis).
class FastMatrixMultiply {
  //! \brief Type for couplings and offsets
  //! \details `first` refers to rows (out) and `second` to columns (input)
  typedef std::vector<std::pair<t_uint, t_uint>> Indices;

public:
  //! Creates the fast matrix multiply object
  //! \param[in] em_background: Electro-magnetic properties of the background material
  //! \param[in] wavenumber: Angular wave-number of the incident plane-wave
  //! \param[in] scatterers: all spherical scatterers in the problem
  //! \param[in] incident: Indices defining the scatterers this instance will take care of.
  //!     In practice, the input vector to the multiplication will contain the coefficients of the
  //!     effective incident field expanded on a spherical basis set at the location of the
  //!     scatterers in this range.
  //! \param[in] translate: Indices defining the location for which to expand the field on output.
  //!     In practice, the output vector of the multiplication will contain the coefficients of
  //!     the
  //!     spherical basis set used to expand the field at the location of the scatterers in this
  //!     range.
  FastMatrixMultiply(ElectroMagnetic const &em_background, t_real wavenumber,
                     std::vector<Scatterer> const &scatterers, Matrix<bool> const &couplings)
      : em_background_(em_background), wavenumber_(wavenumber), scatterers_(scatterers),
        indices_(compute_indices(couplings)),
        incident_offsets_(compute_offsets(scatterers, couplings.colwise().any())),
        translate_offsets_(compute_offsets(scatterers, couplings.rowwise().any())),
        rotations_(compute_rotations(scatterers, couplings)),
        mie_coefficients_(
            compute_mie_coefficients(em_background, wavenumber, scatterers, couplings)),
        coaxial_translations_(compute_coaxial_translations(wavenumber, scatterers, couplings)),
        normalization_(compute_normalization(scatterers)) {}
  FastMatrixMultiply(t_real wavenumber, std::vector<Scatterer> const &scatterers,
                     Matrix<bool> const &couplings)
      : FastMatrixMultiply(ElectroMagnetic(), wavenumber, scatterers, couplings) {}
  FastMatrixMultiply(ElectroMagnetic const &em_background, t_real wavenumber,
                     std::vector<Scatterer> const &scatterers)
      : FastMatrixMultiply(em_background, wavenumber, scatterers,
                           Matrix<bool>::Ones(scatterers.size(), scatterers.size())) {}
  FastMatrixMultiply(t_real wavenumber, std::vector<Scatterer> const &scatterers)
      : FastMatrixMultiply(ElectroMagnetic(), wavenumber, scatterers) {}

  //! Total size of the problem
  t_uint size() const { return rows() * cols(); }

  //! \brief Applies fast matrix multiplication to effective incident field
  void operator()(Vector<t_complex> const &in, Vector<t_complex> &out) const;
  //! \brief Applies fast matrix multiplication to effective incident field
  Vector<t_complex> operator()(Vector<t_complex> const &in) const;
  //! \brief Applies fast matrix multiplication to effective incident field
  Vector<t_complex> operator*(Vector<t_complex> const &in) const { return operator()(in); }

  //! \brief computes transpose operation
  void transpose(Vector<t_complex> const &in, Vector<t_complex> &out) const;
  //! \brief Applies fast matrix multiplication to effective incident field
  Vector<t_complex> transpose(Vector<t_complex> const &in) const;

  //! \brief computes conjugate operation
  void conjugate(Vector<t_complex> const &in, Vector<t_complex> &out) const {
    operator()(in.conjugate(), out);
    out = out.conjugate();
  }
  //! \brief Applies fast matrix multiplication to effective incident field
  Vector<t_complex> conjugate(Vector<t_complex> const &in) const {
    return operator()(in.conjugate()).conjugate();
  }

  //! \brief computes conjugate operation
  void adjoint(Vector<t_complex> const &in, Vector<t_complex> &out) const {
    transpose(in.conjugate(), out);
    out = out.conjugate();
  }
  //! \brief Applies fast matrix multiplication to effective incident field
  Vector<t_complex> adjoint(Vector<t_complex> const &in) const {
    return transpose(in.conjugate()).conjugate();
  }

  //! Number of columns of the fast matrix
  t_uint cols() const { return incident_offsets_.back(); }
  //! Number of columns of the fast matrix
  t_uint rows() const { return translate_offsets_.back(); }

  //! Couplings that this object will compute
  Indices const &couplings() const { return indices_; }

protected:
  static int const nplus = 1;
  //! The properties of the background
  ElectroMagnetic em_background_;
  //! Wavenumber of the incident wave
  t_real wavenumber_;
  //! Scatterers for which to compute matrix product
  std::vector<Scatterer> const scatterers_;
  //! Couplings to compute in this instance
  std::vector<std::pair<t_uint, t_uint>> const indices_;
  //! Offsets for contiguous input vectors
  std::vector<t_uint> const incident_offsets_;
  //! Offsets for contiguous output vectors
  std::vector<t_uint> const translate_offsets_;
  //! Rotations for owned objects
  std::vector<Rotation> const rotations_;
  //! Mie coefficients
  Vector<t_complex> const mie_coefficients_;
  //! Co-axial translations
  std::vector<CachedCoAxialRecurrence::Functor> coaxial_translations_;
  //! Normalization factors between Gumerov and Stout
  Eigen::Array<t_real, Eigen::Dynamic, 2> const normalization_;

  //! Computes index of each particle i in global input vector
  static std::vector<std::pair<t_uint, t_uint>> compute_indices(Matrix<bool> const &couplings);
  //! Computes offsets for output and input vectors
  static std::vector<t_uint>
  compute_offsets(std::vector<Scatterer> const &scatterers, Vector<bool> const &couplings);
  //! Computes rotations between relevant pairs of particles
  static std::vector<Rotation>
  compute_rotations(std::vector<Scatterer> const &scatterers, Matrix<bool> const &couplings);
  //! Computes co-axial translations between relevant pairs of particles
  static std::vector<CachedCoAxialRecurrence::Functor>
  compute_coaxial_translations(t_complex wavenumber_, std::vector<Scatterer> const &scatterers,
                               Matrix<bool> const &couplings);
  //! Computes mie coefficient for each particles
  static Vector<t_complex>
  compute_mie_coefficients(ElectroMagnetic const &background, t_real wavenumber,
                           std::vector<Scatterer> const &scatterers, Matrix<bool> const &couplings);
  //! Normalization factors between Gumerov and Stout
  static Eigen::Array<t_real, Eigen::Dynamic, 2>
  compute_normalization(std::vector<Scatterer> const &scatterers);

  //! Number of basis function for given nmax
  static constexpr t_int nfunctions(t_int nmax) { return nmax * (nmax + 2); }

  //! Max nmax across incident and translate range
  t_int max_nmax() const;

  //! nMax for given incident particle
  t_int incident_nmax(Indices::size_type i) const { return scatterers_[indices_[i].second].nMax; }
  //! nMax for given translate particle
  t_int translate_nmax(Indices::size_type j) const { return scatterers_[indices_[j].first].nMax; }
  //! Data offset of input vector for coupling i
  t_uint incident_offset(Indices::size_type i) const {
    return incident_offsets_[indices_[i].second];
  }
  //! Data offset of output vector for coupling i
  t_uint translate_offset(Indices::size_type j) const {
    return translate_offsets_[indices_[j].first];
  }

  t_real tz(Indices::size_type i) const {
    return (scatterers_[indices_[i].first].vR.toEigenCartesian() -
            scatterers_[indices_[i].second].vR.toEigenCartesian())
        .stableNorm();
  }

  //! True if coupling particle with itself
  bool is_self_interaction(Indices::size_type i) const {
    return indices_[i].first == indices_[i].second;
  }

  //! Adds co-axial rotation/translation for a given particle pair
  template <class T0, class T1, class T2>
  void remove_translation(Eigen::MatrixBase<T0> const &input, Eigen::MatrixBase<T1> const &out,
                          Eigen::MatrixBase<T2> const &work, Indices::size_type i) const;
  //! Adds co-axial rotation/translation for a given particle pair
  template <class T0, class T1, class T2>
  void
  remove_translation_transpose(Eigen::MatrixBase<T0> const &input, Eigen::MatrixBase<T1> const &out,
                               Eigen::MatrixBase<T2> const &work, Indices::size_type i) const;
  //! Apply translation to each particle pair
  void translation(Vector<t_complex> const &in, Vector<t_complex> &out) const;
  //! Apply translation to each particle pair
  void translation_transpose(Vector<t_complex> const &in, Vector<t_complex> &out) const;
};

template <class T0, class T1, class T2>
void FastMatrixMultiply::remove_translation(Eigen::MatrixBase<T0> const &input,
                                            Eigen::MatrixBase<T1> const &out,
                                            Eigen::MatrixBase<T2> const &work,
                                            Indices::size_type i) const {
  auto const in_rows = input.rows();
  auto const out_rows = out.rows();

  auto const max_rows =
      nfunctions(std::lround(std::sqrt(std::max(in_rows, out_rows) + 1)) - 1 + nplus) + 1;
  assert(work.rows() >= max_rows);
  assert(work.cols() == 4);
  const_cast<Eigen::MatrixBase<T2> &>(work).fill(0);

  // work matrices: we will alternatively use one then the other for input and output
  auto alpha = const_cast<Eigen::MatrixBase<T2> &>(work).leftCols(2).topRows(max_rows);
  auto beta = const_cast<Eigen::MatrixBase<T2> &>(work).rightCols(2).topRows(max_rows);

  // First, we take into account Gumerov's very special normalization and notations
  // It adds +/-1 factors, as well as normalization constants
  // But appears when comparing to original calculation
  // We add it second since: (i) the normalization is same for the same n, (ii) it avoids a copy
  // There is no n=0 term at this juncture
  alpha.middleRows(1, in_rows) = input.array() * normalization_.topRows(in_rows).array();

  // Then we apply the rotation - without n=0 term
  rotations_[i](alpha.middleRows(1, in_rows), beta.middleRows(1, in_rows));

  // Then perform co-axial translation - this may create n=0 term
  coaxial_translations_[i](beta, alpha);

  // Then apply field-coaxial-tranlation transform thing - n=0 term may be used to create n=1 term.
  // n=0 term itself becomes zero (thereby choosing a gauge, apparently)
  rotation_coaxial_decomposition(wavenumber_, tz(i), alpha, beta);

  // // Rotate back - remove n=0 term since it is zero
  rotations_[i].adjoint(beta.middleRows(1, out_rows), alpha.middleRows(1, out_rows));

  // Finally, add back into output vector with normalization
  const_cast<Eigen::MatrixBase<T1> &>(out).array() -=
      alpha.middleRows(1, out_rows).array() / normalization_.topRows(out_rows).array();
}

template <class T0, class T1, class T2>
void FastMatrixMultiply::remove_translation_transpose(Eigen::MatrixBase<T0> const &input,
                                                      Eigen::MatrixBase<T1> const &out,
                                                      Eigen::MatrixBase<T2> const &work,
                                                      t_uint i) const {
  auto const in_rows = input.rows();
  auto const out_rows = out.rows();

  auto const max_rows =
      nfunctions(std::lround(std::sqrt(std::max(in_rows, out_rows) + 1)) - 1 + nplus) + 1;
  assert(work.rows() >= max_rows);
  assert(work.cols() == 4);
  const_cast<Eigen::MatrixBase<T2> &>(work).fill(0);

  // work matrices: we will alternatively use one then the other for input and output
  auto alpha = const_cast<Eigen::MatrixBase<T2> &>(work).leftCols(2).topRows(max_rows);
  auto beta = const_cast<Eigen::MatrixBase<T2> &>(work).rightCols(2).topRows(max_rows);

  // First, we take into account Gumerov's very special normalization and notations
  // It adds +/-1 factors, as well as normalization constants
  // But appears when comparing to original calculation
  // We add it second since: (i) the normalization is same for the same n, (ii) it avoids a copy
  // There is no n=0 term at this juncture
  alpha.middleRows(1, in_rows) = input.array() / normalization_.topRows(in_rows).array();

  // Then we apply the rotation - without n=0 term
  rotations_[i].conjugate(alpha.middleRows(1, in_rows), beta.middleRows(1, in_rows));

  // Then apply field-coaxial-tranlation transform thing - n=0 term may be used to create n=1 term.
  // n=0 term itself becomes zero (thereby choosing a gauge, apparently)
  rotation_coaxial_decomposition_transpose(wavenumber_, tz(i), beta, alpha);

  // Then perform co-axial translation - this may create n=0 term
  coaxial_translations_[i].transpose(alpha, beta);

  // Rotate back - remove n=0 term since it is zero
  rotations_[i].transpose(beta.middleRows(1, out_rows), alpha.middleRows(1, out_rows));

  // Finally, add back into output vector
  const_cast<Eigen::MatrixBase<T1> &>(out).array() -=
      alpha.middleRows(1, out_rows).array() * normalization_.topRows(out_rows).array();
}
}

#endif
