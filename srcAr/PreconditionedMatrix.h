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
#include <Eigen/LU>
 
namespace optimet {
//Computes source vector
Vector<t_complex>
source_vector(Geometry const &geometry, std::shared_ptr<Excitation const> incWave);
//Computes source vector from a range of scatterers
Vector<t_complex> source_vector(std::vector<Scatterer>::const_iterator first,
                                std::vector<Scatterer>::const_iterator const &last,
                                std::shared_ptr<Excitation const> incWave);  
                                                                 
#ifdef OPTIMET_MPI                                 
// source vectors needed for SH arbitrary shapes and parallel
Vector<t_complex> source_vectorSH_parallelAR3(Geometry &geometry, int gran1, int gran2,
          std::shared_ptr<Excitation const> incWave, Vector<t_complex> &internalCoef_FF_, Vector<t_complex> &scatteredCoef_FF_, int objIndex);


Vector<t_complex> source_vectorSH_parallelAR1(Geometry &geometry, int gran1, int gran2,
              std::shared_ptr<Excitation const> incWave, Vector<t_complex> &internalCoef_FF_, Vector<t_complex> &scatteredCoef_FF_, int objIndex);

// computes the distributed source vector for SH on many nodes
Vector<t_complex> distributed_source_vector_SH_Mnode(Geometry &geometry,
                                           std::shared_ptr<Excitation const> incWave,
                                           Vector<t_complex> &X_int_, Vector<t_complex> &X_sca_, Matrix<t_complex> &TmatrixSH);

Vector<t_complex> distributed_vector_SH_AR1(Geometry &geometry,
                                           std::shared_ptr<Excitation const> incWave,
                                           Vector<t_complex> &X_int_, Vector<t_complex> &X_sca_);


// computes scattering matrix of many targets at FF
Matrix<t_complex> ScatteringMatrixFF(Matrix<t_complex> &TMatrixFF, Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave);

// computes scattering matrix of many targets at SH
Matrix<t_complex> ScatteringMatrixSH(Matrix<t_complex> &TMatrixSH, Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave);

// computes the Qmatrix for FF, single target
Vector<t_complex>getQmatrix_FF(Geometry const &geometry, ElectroMagnetic const &bground,
                                 std::shared_ptr<Excitation const> incWave, int gran1, int gran2, int objIndex);


// computes the RgQmatrix for FF, single target
Vector<t_complex>getRgQmatrix_FF(Geometry const &geometry, ElectroMagnetic const &bground,
                                 std::shared_ptr<Excitation const> incWave, int gran1, int gran2, int objIndex);

// Computes the Qmatrix for SH, single target
Vector<t_complex>getQmatrix_SH(Geometry const &geometry, ElectroMagnetic const &bground,
                                 std::shared_ptr<Excitation const> incWave, int gran1, int gran2, int objIndex);

// Computes the RgQmatrix for SH, single target
Vector<t_complex>getRgQmatrix_SH(Geometry const &geometry, ElectroMagnetic const &bground,
                                 std::shared_ptr<Excitation const> incWave, int gran1, int gran2, int objIndex);
                                                 
//! Computes Tmatrix and RgQmatrix in paralllel, single target, FF
Matrix<t_complex> getTRgQmatrix_FF_parr(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave
                                                   );

//! Computes Tmatrix and RgQmatrix in paralllel, single target, SH
Matrix<t_complex> getTRgQmatrix_SH_parr(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave
                                                   );
#endif
}
#endif
