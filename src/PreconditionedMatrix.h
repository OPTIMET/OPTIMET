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
#include "mpi/Communicator.h"
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
                                
   
 Vector<t_complex> source_vectorSH(Geometry &geometry, std::shared_ptr<Excitation const> incWave, 
                                   Vector<t_complex> &internalCoef_FF_, Vector<t_complex> &scatteredCoef_FF_,
                                    std::vector<double *> CGcoeff);  
                                
 //! \brief Computes source SH vector from a range of scatterers                               
Vector<t_complex> source_vectorSH(Geometry &geometry,std::vector<Scatterer>::const_iterator first,
                                std::vector<Scatterer>::const_iterator const &last,
                                std::shared_ptr<Excitation const> incWave, Vector<t_complex> &internalCoef_FF_, 
                                Vector<t_complex> &scatteredCoef_FF_, std::vector<double *> CGcoeff);

// source vector needed for SH arbitrary shapes
Vector<t_complex> source_vectorSHarb1(Geometry &geometry, std::shared_ptr<Excitation const> incWave, 
                                 Vector<t_complex> &internalCoef_FF_, Vector<t_complex> &scatteredCoef_FF_, std::vector<double*>CGcoeff);                                 

//! \brief Computes source SH vector from a range of scatterers adapted for parallelization                              
Vector<t_complex> source_vectorSH_parallel(Geometry &geometry, int gran1, int gran2,
                                std::shared_ptr<Excitation const> incWave, Vector<t_complex> &internalCoef_FF_, std::vector<double *> CGcoeff);                                 
// source vectors needed for SH arbitrary shapes and parallel
Vector<t_complex> source_vectorSH_parallelAR3(Geometry &geometry, int gran1, int gran2,
                                std::shared_ptr<Excitation const> incWave, Vector<t_complex> &internalCoef_FF_, Vector<t_complex> &scatteredCoef_FF_,
                                  std::vector<double *> CGcoeff);


Vector<t_complex> source_vectorSH_parallelAR1(Geometry &geometry, int gran1, int gran2,
                                std::shared_ptr<Excitation const> incWave, Vector<t_complex> &internalCoef_FF_, Vector<t_complex> &scatteredCoef_FF_,
                                  std::vector<double *> CGcoeff);

// computes the distributed source vector for SH on many nodes
Vector<t_complex> distributed_source_vector_SH_Mnode(Geometry &geometry,
                                           std::shared_ptr<Excitation const> incWave,
                                           Vector<t_complex> &X_int_, Vector<t_complex> &X_sca_, std::vector<double *> CGcoeff);

Vector<t_complex> distributed_vector_SH_AR1(Geometry &geometry,
                                           std::shared_ptr<Excitation const> incWave,
                                           Vector<t_complex> &X_int_, Vector<t_complex> &X_sca_, std::vector<double *> CGcoeff);                      

//! Computes preconditioned scattering matrix in serial for fundamental frequency
Matrix<t_complex> preconditioned_scattering_matrix(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave);
                                                                                                  
 
//! Computes preconditioned scattering matrix in serial for SH frequency, scattering coefficients
 
Matrix<t_complex> preconditioned_scattering_matrixSH(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave); 
                                                                                                                                                                                                        
                                                 
//! Computes preconditioned scattering matrix in paralllel
Matrix<t_complex> preconditioned_scattering_matrix(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave,
                                                   scalapack::Context const &context,
                                                   scalapack::Sizes const &blocks);

//! Computes preconditioned scattering matrix for SH in paralllel
Matrix<t_complex> preconditioned_scattering_matrix_SH(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave,
                                                   scalapack::Context const &context,
                                                   scalapack::Sizes const &blocks);
//! Distributes the source vectors
Vector<t_complex> distributed_source_vector(Vector<t_complex> const &input,
                                            scalapack::Context const &context,
                                            scalapack::Sizes const &blocks);


Vector<t_complex> distributed_source_vector_SH(Geometry &geometry, Vector<t_complex> &VecMnod,
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
