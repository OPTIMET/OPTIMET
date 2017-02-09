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
//! \brief Computes source vector from a range of scatterers
Vector<t_complex> source_vector(std::vector<Scatterer>::const_iterator first,
                                std::vector<Scatterer>::const_iterator const &last,
                                std::shared_ptr<Excitation const> incWave);
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

#ifdef OPTIMET_SCALAPACK
//! Gather the distributed vector into a single vector
Vector<t_complex> gather_all_source_vector(t_uint n, Vector<t_complex> const &input,
                                           scalapack::Context const &context,
                                           scalapack::Sizes const &blocks);
//! Gather the distributed vector into a single vector
Vector<t_complex> gather_all_source_vector(scalapack::Matrix<t_complex> const &matrix);
#endif
}
#endif
