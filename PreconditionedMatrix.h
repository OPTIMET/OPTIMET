#ifndef OPTIMET_PRECONDITIONNED_MATRIX_H
#define OPTIMET_PRECONDITIONNED_MATRIX_H

#include "Excitation.h"
#include "Geometry.h"
#include "Types.h"
#include "scalapack/Context.h"
#include "scalapack/Matrix.h"

namespace optimet {
//! \brief Computes source vector
Vector<t_complex>
source_vector(Geometry const &geometry, std::shared_ptr<Excitation const> incWave);
//! \brief Computes source vector from fundamental frequency
Vector<t_complex> local_source_vector(Geometry const &geometry,
                                      std::shared_ptr<Excitation const> incWave,
                                      Vector<t_complex> const &input_coeffs);

//! Computes preconditioned scattering matrix
Matrix<t_complex> preconditioned_scattering_matrix(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave);

//! Computes preconditioned scattering matrix in paralllel
Matrix<t_complex> preconditioned_scattering_matrix(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave,
                                                   scalapack::Context const &context,
                                                   scalapack::Sizes const &blocks);
//! Distributes the source vectors
Vector<t_complex> distributed_source_vector(Vector<t_complex> const &input,
                                            scalapack::Context const &context,
                                            scalapack::Sizes const &blocks);

#ifdef OPTIMET_MPI
//! Gather the distributed vector into a single vector
Vector<t_complex> gather_all_source_vector(t_uint n, Vector<t_complex> const &input,
                                           scalapack::Context const &context,
                                           scalapack::Sizes const &blocks);
//! Gather the distributed vector into a single vector
Vector<t_complex> gather_all_source_vector(scalapack::Matrix<t_complex> const &matrix);
#endif
}
#endif
