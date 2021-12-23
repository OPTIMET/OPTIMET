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

struct Matrix_ACA{
Matrix<t_complex> U;
Matrix<t_complex> V;
Matrix<t_complex> S_sub;
int dim;
};

//! Computes source vector at FF
Vector<t_complex>
source_vector(Geometry const &geometry, std::shared_ptr<Excitation const> incWave);
//! Computes source vector from a range of scatterers
Vector<t_complex> source_vector(std::vector<Scatterer>::const_iterator first,
                                std::vector<Scatterer>::const_iterator const &last,
                                std::shared_ptr<Excitation const> incWave, Geometry const &geometry);
                                
   
 Vector<t_complex> source_vectorSH(Geometry &geometry, std::shared_ptr<Excitation const> incWave, 
                                   Vector<t_complex> &internalCoef_FF_, Vector<t_complex> &scatteredCoef_FF_,
                                    std::vector<double *> CGcoeff);  
                                
 //! Computes source SH vector from a range of scatterers                               
Vector<t_complex> source_vectorSH(Geometry &geometry,std::vector<Scatterer>::const_iterator first,
                                std::vector<Scatterer>::const_iterator const &last,
                                std::shared_ptr<Excitation const> incWave, Vector<t_complex> &internalCoef_FF_, 
                                Vector<t_complex> &scatteredCoef_FF_, std::vector<double *> CGcoeff);

// the second part of the SH source vector
Vector<t_complex> source_vectorSH_K1ana(Geometry &geometry,std::shared_ptr<Excitation const> incWave,
Vector<t_complex> &internalCoef_FF_, Vector<t_complex> &scatteredCoef_FF_, std::vector<double *> CGcoeff);
                                 

//! Computes source SH vector from a range of scatterers adapted for parallelization
#ifdef OPTIMET_MPI                              
Vector<t_complex> source_vectorSH_parallel(Geometry &geometry, int gran1, int gran2,
                                std::shared_ptr<Excitation const> incWave, Vector<t_complex> &internalCoef_FF_, std::vector<double *> CGcoeff);         
#endif
// computes the distributed source vector for SH on many compute nodes
#ifdef OPTIMET_MPI
Vector<t_complex> distributed_source_vector_SH_Mnode(Geometry &geometry,
                                           std::shared_ptr<Excitation const> incWave,
                                           Vector<t_complex> &X_int_, Vector<t_complex> &X_sca_, std::vector<double *> CGcoeff);
#endif
// computation of the second part of SH source vector in parallel
#ifdef OPTIMET_MPI
Vector<t_complex> source_vectorSH_K1ana_parallel(Geometry &geometry,
                                           std::shared_ptr<Excitation const> incWave,
                                           Vector<t_complex> &X_int_, Vector<t_complex> &X_sca_, std::vector<double *> CGcoeff);
#endif 

//! Computes preconditioned scattering matrix in serial for fundamental frequency
Matrix<t_complex> preconditioned_scattering_matrix(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave);

// Computes the scattering matrix at FF with ACA compression in parallel
#ifdef OPTIMET_MPI
void Scattering_matrix_ACA_FF_parallel(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave, std::vector<Matrix_ACA>& S_comp);

// Computes scattering matrix at SH with ACA compression in parallel
void Scattering_matrix_ACA_SH_parallel(Geometry const &geometry,
                         std::shared_ptr<Excitation const> incWave, std::vector<Matrix_ACA>& S_comp); 
#endif
// Computes scattering matrix at FF with ACA compression serial
void Scattering_matrix_ACA_FF(Geometry const &geometry,
                                                    std::shared_ptr<Excitation const> incWave, std::vector<Matrix_ACA>& S_comp);

// Computes scattering matrix at SH with ACA compression serial
void Scattering_matrix_ACA_SH(Geometry const &geometry,
                         std::shared_ptr<Excitation const> incWave, std::vector<Matrix_ACA>& S_comp);

// ACA algorithm, adaptive cross approximation
void ACA_compression(Matrix<t_complex> &U, Matrix<t_complex> &V, Matrix<t_complex> &CoupMat);

//search for the index of the largest absolute element in row/column for ACA algorithm
int getMaxInd(Vector<t_complex> &RowCol, Vector<int> &K, int kmax);

//the gmres solver for ACA compressed matrices with restarts                                                   
Vector<t_complex> Gmres_Zcomp(std::vector<Matrix_ACA>const &S_comp, Vector<t_complex>const &Y, double tol, int maxit, int no_rest, Geometry const &geometry);  

// matrix-vector product for compressed matrices in parallel
#ifdef OPTIMET_MPI
Vector<t_complex> matvec_parallel(std::vector<Matrix_ACA>const &S_comp, Vector<t_complex> &J, Geometry const &geometry);
#endif

// matrix-vector product for compressed matrices in serial
Vector<t_complex> matvec(std::vector<Matrix_ACA>const &S_comp, Vector<t_complex> &J, Geometry const &geometry);

Matrix<t_complex> det_approx (double beta, int n, Matrix<t_complex> &H);
                                                                                                  
//! Computes preconditioned scattering matrix in serial for SH frequency 
Matrix<t_complex> preconditioned_scattering_matrixSH(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave); 
#ifdef OPTIMET_SCALAPACK                                                                                                                         //! Computes preconditioned scattering matrix at FF in paralllel
Matrix<t_complex> preconditioned_scattering_matrix(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave,
                                                   scalapack::Context const &context,
                                                   scalapack::Sizes const &blocks);

//! Computes preconditioned scattering matrix at SH in paralllel
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
#endif

#ifdef OPTIMET_SCALAPACK
//! Gather the distributed vector into a single vector
Vector<t_complex> gather_all_source_vector(t_uint n, Vector<t_complex> const &input,
                                           scalapack::Context const &context,
                                           scalapack::Sizes const &blocks);

Vector<t_complex> gather_all_source_vector(scalapack::Matrix<t_complex> const &matrix);

#endif
}
#endif
