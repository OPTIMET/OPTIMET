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

#include "Coupling.h"
#include "PreconditionedMatrix.h"
#include "Types.h"
#include "scalapack/BroadcastToOutOfContext.h"
#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Sparse>
using namespace std::chrono;

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



Vector<t_complex> distributed_source_vector_SH(Geometry &geometry, Vector<t_complex> &VecMnod,
                                           scalapack::Context const &context,
                                            scalapack::Sizes const &blocks) {

  auto const nobj = geometry.objects.size();
  if(nobj == 0)
     return Vector<t_complex>::Zero(0);

  auto const serial = context.serial();
  t_uint const N(VecMnod.size());
  Eigen::Map<Matrix<t_complex> const> const input_map(VecMnod.data(), serial.is_valid() ? N : 0,
                                                      serial.is_valid() ? 1 : 0);
  scalapack::Matrix<t_complex const *> const serial_vector(input_map, serial, {N, 1}, {N, 1});
  scalapack::Matrix<t_complex> dist_vector(context, {N, 1}, blocks);
  serial_vector.transfer_to(context, dist_vector);
  if(dist_vector.local().cols() == 0)
   return Vector<t_complex>::Zero(0);
  
 return dist_vector.local();
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


#ifdef OPTIMET_MPI
Vector<t_complex> source_vectorSH_K1ana_parallel(Geometry &geometry,
                                           std::shared_ptr<Excitation const> incWave,
                                           Vector<t_complex> &X_int_, Vector<t_complex> &X_sca_, std::vector<double *> CGcoeff) {
  auto const nobj = geometry.objects.size();
  if(nobj == 0)
     return Vector<t_complex>::Zero(0);

  int gran, gran1, gran2, kk, objInd;
  mpi::Communicator communicator;
  int rank = communicator.rank();
  int size = communicator.size();
  auto const nMaxS = geometry.objects.front().nMaxS;

  t_uint const pMax = nMaxS * (nMaxS + 2);

  int TMax = nobj * pMax;
  
  Vector<t_complex> result, resultK, resultKK, result3, result1;
  Vector<t_complex> X_int_proc;

  int sizeFFint;

    if (rank==0){
    sizeFFint =X_int_.size();
    X_int_proc = X_int_;
    }

    MPI_Bcast(&sizeFFint, 1, MPI_INT, 0, MPI_COMM_WORLD);
    X_int_proc.resize(sizeFFint); // broadcasting internal and scattering FF field coeff
    MPI_Bcast(&X_int_proc(0), sizeFFint, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

   if (rank < (TMax % size)) {
    gran1 = rank * (TMax/size + 1);
    gran2 = gran1 + TMax/size + 1;
    } else {
    gran1 = rank * (TMax/size) + (TMax % size);
    gran2 = gran1 + (TMax/size);
    }

   // Analytical for spheres
   if (geometry.objects[0].scatterer_type == "sphere"){

    int sizeVec = 4 * (gran2 - gran1);

     result.resize(4*nobj*pMax);
     resultKK.resize(4*nobj*pMax);
     resultKK.setZero();

     resultK.resize(2*nobj*pMax);
     resultK.setZero();

    Vector<t_complex> resultProc(sizeVec);

    Vector<int> sizesProc(size), disps(size);

    MPI_Allgather (&sizeVec, 1, MPI_INT, &sizesProc(0), 1, MPI_INT, MPI_COMM_WORLD);

   for (int kk = 0; kk < size; kk++)
   disps(kk) = (kk > 0) ? (disps(kk-1) + sizesProc(kk-1)) : 0;

    resultProc = source_vectorSH_parallel(geometry, gran1, gran2, incWave, X_int_proc, CGcoeff);

    MPI_Gatherv (&resultProc(0), sizeVec, MPI_DOUBLE_COMPLEX, &result(0), &sizesProc(0), &disps(0), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);


  if (rank == 0){

  for (int ranki = 0; ranki < size; ranki++){
  
    if (ranki < (TMax % size)) {
    gran1 = ranki * (TMax/size + 1);
    gran2 = gran1 + TMax/size + 1;
    } else {
    gran1 = ranki * (TMax/size) + (TMax % size);
    gran2 = gran1 + (TMax/size);
    }

      sizeVec = 4 * (gran2 - gran1);

  int brojac (0);
  for (int ii = gran1; ii < gran2; ii++){

  objInd = ii / pMax;

  kk = ii - objInd * pMax;

  gran = objInd * 4 * pMax + kk;

  resultKK(gran) = result(disps(ranki) + brojac);
  resultKK(gran + pMax) = result(disps(ranki) + (sizeVec/4) + brojac);
  resultKK(gran + 2*pMax) = result(disps(ranki) + 2*(sizeVec/4) + brojac);
  resultKK(gran + 3*pMax) = result(disps(ranki) + 3*(sizeVec/4) + brojac);

  brojac++;

  } // for ii

  }  // for ranki


  int i1 = 0;
  int ii1 = 0;

   for(int kk = 0; kk != nobj; kk++)  {

    resultK.segment(ii1, 2*pMax) =
         (geometry.objects[kk].getIauxSH2(incWave->omega(), geometry.bground).array() * resultKK.segment(i1+2*pMax, 2*pMax).array());

  i1 += 4 * pMax;
  ii1 += 2 * pMax;

  }
  
  } // if rank0
  MPI_Bcast(&resultK(0), 2 * nobj *pMax, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  
  } // if sphere

return resultK;

}

 
Vector<t_complex> distributed_source_vector_SH_Mnode(Geometry &geometry,
                                           std::shared_ptr<Excitation const> incWave,
                                           Vector<t_complex> &X_int_, Vector<t_complex> &X_sca_, std::vector<double *> CGcoeff) {
auto const nobj = geometry.objects.size();
  if(nobj == 0)
     return Vector<t_complex>::Zero(0);
  
  int gran, gran1, gran2, kk, objInd;
  mpi::Communicator communicator;
  int rank = communicator.rank();
  int size = communicator.size();
  auto const nMaxS = geometry.objects.front().nMaxS;
  auto const k_b_SH = 2.0 * incWave->omega() * std::sqrt(geometry.bground.epsilon * geometry.bground.mu);

  t_uint const pMax = nMaxS * (nMaxS + 2);

  int TMax = nobj * pMax;

  Vector<t_complex> result, resultK, resultKK, result3, result1, resultK3, resultK1;
  Matrix<t_complex> TmatrixSH (2*pMax, 2*pMax);
  Vector<t_complex> X_int_proc, X_sca_proc;
  
  int sizeFFint;
    
    if (rank==0){
    sizeFFint =X_int_.size();
    X_int_proc = X_int_;
    X_sca_proc = X_sca_;
    }
    
    MPI_Bcast(&sizeFFint, 1, MPI_INT, 0, MPI_COMM_WORLD);
    X_int_proc.resize(sizeFFint); // broadcasting internal and scattering FF field coeff
    X_sca_proc.resize(sizeFFint);
    MPI_Bcast(&X_int_proc(0), sizeFFint, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&X_sca_proc(0), sizeFFint, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    
    if (rank < (TMax % size)) {
    gran1 = rank * (TMax/size + 1);
    gran2 = gran1 + TMax/size + 1;
    } else {
    gran1 = rank * (TMax/size) + (TMax % size);
    gran2 = gran1 + (TMax/size);
    }

   // Analytical for spheres
    if (geometry.objects[0].scatterer_type == "sphere"){

    int sizeVec = 4 * (gran2 - gran1);

     result.resize(4*nobj*pMax);
     resultKK.resize(4*nobj*pMax);
     resultKK.setZero();

     resultK.resize(2*nobj*pMax);
     resultK.setZero();
          
    Vector<t_complex> resultProc(sizeVec);

    Vector<int> sizesProc(size), disps(size);
  
    MPI_Allgather (&sizeVec, 1, MPI_INT, &sizesProc(0), 1, MPI_INT, MPI_COMM_WORLD);
  
   for (int kk = 0; kk < size; kk++)
   disps(kk) = (kk > 0) ? (disps(kk-1) + sizesProc(kk-1)) : 0;
      
    resultProc = source_vectorSH_parallel(geometry, gran1, gran2, incWave, X_int_proc, CGcoeff);

    MPI_Gatherv (&resultProc(0), sizeVec, MPI_DOUBLE_COMPLEX, &result(0), &sizesProc(0), &disps(0), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD); 


  if (rank == 0){ 
  
  for (int ranki = 0; ranki < size; ranki++){

     if (ranki < (TMax % size)) {
    gran1 = ranki * (TMax/size + 1);
    gran2 = gran1 + TMax/size + 1;
    } else {
    gran1 = ranki * (TMax/size) + (TMax % size);
    gran2 = gran1 + (TMax/size);
    }

      sizeVec = 4 * (gran2 - gran1);

  int brojac (0);
  for (int ii = gran1; ii < gran2; ii++){

  objInd = ii / pMax;

  kk = ii - objInd * pMax; 

  gran = objInd * 4 * pMax + kk;

  resultKK(gran) = result(disps(ranki) + brojac);
  resultKK(gran + pMax) = result(disps(ranki) + (sizeVec/4) + brojac);
  resultKK(gran + 2*pMax) = result(disps(ranki) + 2*(sizeVec/4) + brojac);
  resultKK(gran + 3*pMax) = result(disps(ranki) + 3*(sizeVec/4) + brojac);

  brojac++; 
  
  } // for ii

  }  // for ranki
 

  int i1 = 0;
  int ii1 = 0;
 
   for(int kk = 0; kk != nobj; kk++)  {
 
    resultK.segment(ii1, 2*pMax) =
         (geometry.objects[kk].getTLocalSH1_outer(incWave->omega(), geometry.bground).array() * resultKK.segment(i1, 2*pMax).array()) +

        (geometry.objects[kk].getTLocalSH2_outer(incWave->omega(), geometry.bground).array() * resultKK.segment(2*pMax + i1, 2*pMax).array());
       
  i1 += 4 * pMax;
  ii1 += 2 * pMax;

  }
   } // if rank0

   MPI_Bcast(&resultK(0), 2 * nobj *pMax, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  } // if sphere

return resultK;

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
  Matrix<t_complex> Tmatrix (2 * n , 2 * n); 
   
  if(first == end_first or second == end_second)
    return Matrix<t_complex>::Zero(2 * n * (end_first - first), 2 * n * (end_second - second));

  Matrix<t_complex> result(2 * n * (end_first - first), 2 * n * (end_second - second));
  
   size_t x(0);
   for(auto iteri(first); iteri != end_first; ++iteri, x += 2 * n) {

     iteri->getTLocal(Tmatrix, incWave->omega(), bground);

     size_t y(0);
     for(auto iterj(second); iterj != end_second; ++iterj, y += 2 * n) {

    
      if(iteri == iterj) {
        
        result.block(x, y, 2 * n, 2 * n) = Matrix<t_complex>::Identity(2 * n, 2 * n);

        
      } else {
      
        Coupling const AB(iteri->vR - iterj->vR, incWave->waveK, nMax);
 
        result.block(x, y, n, n) = AB.diagonal.transpose();
        result.block(x + n, y + n, n, n) = AB.diagonal.transpose();
        result.block(x, y + n, n, n) = AB.offdiagonal.transpose();
        result.block(x + n, y, n, n) = AB.offdiagonal.transpose();
        result.block(x, y, 2 * n, 2 * n) = -(Tmatrix)*result.block(x, y, 2 * n, 2 * n);
        
  
    }
  
      }
  
        }
  
  return result;
}
#ifdef OPTIMET_MPI
void Scattering_matrix_ACA_FF_parallel(Geometry const &geometry, std::shared_ptr<Excitation const> incWave, std::vector<Matrix_ACA> &S_comp){

 
  auto const nMax = geometry.objects[0].nMax;
  auto const n = nMax * (nMax + 2);
  Matrix<t_complex> Tmatrix (2*n , 2*n), U, V;
  Matrix<t_complex> CoupSubm (2*n, 2*n);
  double distance, sizeMAT(0.0);
  int nobj = geometry.objects.size();
  
  mpi::Communicator communicator;
  int rank = communicator.rank();
  int size = communicator.size();
  int gran1, gran2, Ncp_proc;
  Vector<double> sizeMAT_vec(size);
  
  if (rank < (nobj % size)) {
    gran1 = rank * (nobj/size + 1);
    gran2 = gran1 + nobj/size + 1;
    } else {
    gran1 = rank * (nobj/size) + (nobj % size);
    gran2 = gran1 + (nobj/size);
    }

   Ncp_proc = nobj*(gran2 - gran1);
   S_comp.resize(Ncp_proc);
   
   int brojac(0);
   
 
  for(int ii = gran1; ii < gran2; ++ii) {
  
     geometry.objects[ii].getTLocal(Tmatrix, incWave->omega(), geometry.bground);

      for(int jj = 0; jj != nobj; ++jj) {
       
    
      if(ii == jj) {
        
        S_comp[brojac].S_sub= Matrix<t_complex>::Identity(2 * n, 2 * n);
        S_comp[brojac].dim = 2 * n;
        sizeMAT = sizeMAT + (S_comp[brojac].S_sub.size())*(16.0/1e6);
        
      } else {
      
        distance = Tools::findDistance(geometry.objects[ii].vR, geometry.objects[jj].vR);
      
        Coupling const AB(geometry.objects[ii].vR - geometry.objects[jj].vR, incWave->waveK, nMax);
        
        CoupSubm.block(0, 0, n, n) = AB.diagonal.transpose();
        CoupSubm.block(n, n, n, n) = AB.diagonal.transpose();
        CoupSubm.block(0, n, n, n) = AB.offdiagonal.transpose();
        CoupSubm.block(n, 0, n, n) = AB.offdiagonal.transpose();
        CoupSubm = - Tmatrix * CoupSubm;
        

        if (distance >= 4.0*(geometry.objects[ii].radius + geometry.objects[jj].radius)){ //admissibility criterion for ACA
        
        ACA_compression(U , V, CoupSubm);
        
        S_comp[brojac].U = U;
        S_comp[brojac].V = V;
        S_comp[brojac].dim = 2 * n;
        sizeMAT = sizeMAT + (S_comp[brojac].U.size())*(16.0/1e6) + (S_comp[brojac].V.size())*(16.0/1e6); 
        }
        
        else{
        S_comp[brojac].S_sub =  CoupSubm;
        S_comp[brojac].dim = 2 * n;
        sizeMAT = sizeMAT + (S_comp[brojac].S_sub.size())*(16.0/1e6);
        } 
        

   }

  brojac++;
  } 

 }
// sum all the partial sizes of matrices
MPI_Gather(&sizeMAT, 1, MPI_DOUBLE, &sizeMAT_vec(0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

if(rank==0)
std::cout<<"The size of the FF scattering matrix in MB is"<< sizeMAT_vec.sum()<<std::endl;
 
}
#endif
void Scattering_matrix_ACA_FF(Geometry const &geometry, std::shared_ptr<Excitation const> incWave, std::vector<Matrix_ACA> &S_comp){

 
  auto const nMax = geometry.objects[0].nMax;
  auto const n = nMax * (nMax + 2);
  Matrix<t_complex> Tmatrix (2*n , 2*n), U, V;
  Matrix<t_complex> CoupSubm (2*n, 2*n);
  double distance, sizeMAT(0.0);
  int nobj = geometry.objects.size();
  S_comp.resize(nobj*nobj);
  
  for(int ii = 0; ii != nobj; ++ii) {
      
     geometry.objects[ii].getTLocal(Tmatrix, incWave->omega(), geometry.bground);
     
    
    for(int jj = 0; jj != nobj; ++jj) {

    
      if(ii == jj) {
        S_comp[nobj*ii + jj].S_sub= Matrix<t_complex>::Identity(2 * n, 2 * n);
        S_comp[nobj*ii + jj].dim = 2 * n;
        sizeMAT = sizeMAT + (S_comp[ii*nobj + jj].S_sub.size())*(16.0/1e6);
        
      } else {
      
        distance = Tools::findDistance(geometry.objects[ii].vR, geometry.objects[jj].vR);
      
        Coupling const AB(geometry.objects[ii].vR - geometry.objects[jj].vR, incWave->waveK, nMax);
        
        CoupSubm.block(0, 0, n, n) = AB.diagonal.transpose();
        CoupSubm.block(n, n, n, n) = AB.diagonal.transpose();
        CoupSubm.block(0, n, n, n) = AB.offdiagonal.transpose();
        CoupSubm.block(n, 0, n, n) = AB.offdiagonal.transpose();
        CoupSubm = - (Tmatrix) * CoupSubm;
        

        if (distance >= 4.0*(geometry.objects[ii].radius + geometry.objects[jj].radius)){ //admissibility criterion for ACA
        
        ACA_compression(U , V, CoupSubm);
        
        S_comp[nobj*ii + jj].U = U;
        S_comp[nobj*ii + jj].V = V;
        S_comp[nobj*ii + jj].dim = 2 * n; 
        sizeMAT = sizeMAT + (S_comp[ii*nobj + jj].U.size())*(16.0/1e6) + (S_comp[ii*nobj + jj].V.size())*(16.0/1e6);
        }
        
        else{
        S_comp[nobj*ii + jj].S_sub =  CoupSubm; 
        S_comp[nobj*ii + jj].dim = 2 * n;
        sizeMAT = sizeMAT + (S_comp[ii*nobj + jj].S_sub.size())*(16.0/1e6);
        }
        

    }
  
      }
  
        }
        
  std::cout<<"The size of the FF matrix in MB is"  <<sizeMAT<< std::endl;      
  
}



Matrix<t_complex>
preconditioned_scattering_matrixSH(std::vector<Scatterer>::const_iterator const &first,
                                 std::vector<Scatterer>::const_iterator const &end_first,
                                 std::vector<Scatterer>::const_iterator const &second,
                                 std::vector<Scatterer>::const_iterator const &end_second,
                                 ElectroMagnetic const &bground,
                                 std::shared_ptr<Excitation const> incWave) {
                                 
                              
  auto const nMaxS = first->nMaxS;
  auto const n = nMaxS * (nMaxS + 2);
  Matrix<t_complex> TmatrixSH (2*n , 2*n);
  Matrix<t_complex> resultSH;

  if (first->scatterer_type == "sphere"){

  resultSH.resize(2 * n * (end_first - first), 2 * n * (end_second - second));

  if(first == end_first or second == end_second)
  return Matrix<t_complex>::Zero(2 * n * (end_first - first), 2 * n * (end_second - second));
        
    size_t x(0);
    for(auto iteri(first); iteri != end_first; ++iteri, x += 2 * n) {
    iteri->getTLocalSH(TmatrixSH, incWave->omega(), bground);

    size_t y(0);
    for(auto iterj(second); iterj != end_second; ++iterj, y += 2 * n) {

      if(iteri == iterj) {

        resultSH.block(x, y, 2 * n, 2 * n) = Matrix<t_complex>::Identity(2 * n, 2 * n);
  
      } 
      
      else {


        Coupling const AB(iteri->vR - iterj->vR, 2.0 * incWave->waveK, nMaxS);
 
        resultSH.block(x, y, n, n) = AB.diagonal.transpose();
        resultSH.block(x + n, y + n, n, n) = AB.diagonal.transpose();
        resultSH.block(x, y + n, n, n) = AB.offdiagonal.transpose();
        resultSH.block(x + n, y, n, n) = AB.offdiagonal.transpose();
        resultSH.block(x, y, 2 * n, 2 * n) = - (TmatrixSH)*resultSH.block(x, y, 2 * n, 2 * n);
      
     }
     
   }
   
 }
  
} //if
  
  return resultSH;
  
}

#ifdef OPTIMET_MPI
void Scattering_matrix_ACA_SH_parallel(Geometry const &geometry, std::shared_ptr<Excitation const> incWave, std::vector<Matrix_ACA> &S_comp){

  auto const nMaxS = geometry.objects[0].nMaxS;
  auto const n = nMaxS * (nMaxS + 2);
  Matrix<t_complex> TmatrixSH (2*n , 2*n), U, V;
  Matrix<t_complex> CoupSubmSH (2*n, 2*n);
  double distance, sizeMAT(0.0);
  int nobj = geometry.objects.size();

  mpi::Communicator communicator;
  int rank = communicator.rank();
  int size = communicator.size();
  int gran1, gran2, Ncp_proc;
  Vector<double> sizeMAT_vec(size);

  if (rank < (nobj % size)) {
    gran1 = rank * (nobj/size + 1);
    gran2 = gran1 + nobj/size + 1;
    } else {
    gran1 = rank * (nobj/size) + (nobj % size);
    gran2 = gran1 + (nobj/size);
    }

   Ncp_proc = nobj*(gran2 - gran1);
   S_comp.resize(Ncp_proc);
   int brojac(0);
  
 
  for(int ii = gran1; ii < gran2; ++ii) {
      
     geometry.objects[ii].getTLocalSH(TmatrixSH, incWave->omega(), geometry.bground);
     
    
    for(int jj = 0; jj != nobj; ++jj) {

    
      if(ii == jj) {
        S_comp[brojac].S_sub= Matrix<t_complex>::Identity(2 * n, 2 * n);
        S_comp[brojac].dim = 2 * n;
        sizeMAT = sizeMAT + (S_comp[brojac].S_sub.size())*(16.0/1e6);
        
      } else {
      
        distance = Tools::findDistance(geometry.objects[ii].vR, geometry.objects[jj].vR);
      
        Coupling const AB(geometry.objects[ii].vR - geometry.objects[jj].vR, 2.0 * incWave->waveK, nMaxS);
        
        CoupSubmSH.block(0, 0, n, n) = AB.diagonal.transpose();
        CoupSubmSH.block(n, n, n, n) = AB.diagonal.transpose();
        CoupSubmSH.block(0, n, n, n) = AB.offdiagonal.transpose();
        CoupSubmSH.block(n, 0, n, n) = AB.offdiagonal.transpose();
        CoupSubmSH = - (TmatrixSH) * CoupSubmSH;
        

        if (distance >= 4.0*(geometry.objects[ii].radius + geometry.objects[jj].radius)){ //admissibility criterion for ACA

        ACA_compression(U , V, CoupSubmSH);
        
        S_comp[brojac].U = U;
        S_comp[brojac].V = V;
        S_comp[brojac].dim = 2 * n;
        sizeMAT = sizeMAT + (S_comp[brojac].U.size())*(16.0/1e6) + (S_comp[brojac].V.size())*(16.0/1e6);
         
        }
        
        else{
        S_comp[brojac].S_sub =  CoupSubmSH;
        S_comp[brojac].dim = 2 * n;
        sizeMAT = sizeMAT + (S_comp[brojac].S_sub.size())*(16.0/1e6);   
        } 
        

   }
  brojac++;
  }
  
 }

MPI_Gather(&sizeMAT, 1, MPI_DOUBLE, &sizeMAT_vec(0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

if(rank==0)
std::cout<<"The size of the SH scattering matrix in MB is"<<sizeMAT_vec.sum()<<std::endl;
  
}
#endif

void Scattering_matrix_ACA_SH(Geometry const &geometry, std::shared_ptr<Excitation const> incWave, std::vector<Matrix_ACA> &S_comp){

  auto const nMaxS = geometry.objects[0].nMaxS;
  auto const n = nMaxS * (nMaxS + 2);
  Matrix<t_complex> TmatrixSH (2*n , 2*n), U, V;
  Matrix<t_complex> CoupSubmSH (2*n, 2*n);
  double distance, sizeMAT(0.0);
  int nobj = geometry.objects.size();
  S_comp.resize(nobj*nobj);
  
 
  for(int ii = 0; ii != nobj; ++ii) {
      
     geometry.objects[ii].getTLocalSH(TmatrixSH, incWave->omega(), geometry.bground);
     
    
    for(int jj = 0; jj != nobj; ++jj) {

    
      if(ii == jj) {
        S_comp[nobj*ii + jj].S_sub= Matrix<t_complex>::Identity(2 * n, 2 * n);
        S_comp[nobj*ii + jj].dim = 2 * n;
        sizeMAT = sizeMAT + (S_comp[ii*nobj + jj].S_sub.size())*(16.0/1e6);
        
      } else {
      
        distance = Tools::findDistance(geometry.objects[ii].vR, geometry.objects[jj].vR);
      
        Coupling const AB(geometry.objects[ii].vR - geometry.objects[jj].vR, 2.0 * incWave->waveK, nMaxS);
        
        CoupSubmSH.block(0, 0, n, n) = AB.diagonal.transpose();
        CoupSubmSH.block(n, n, n, n) = AB.diagonal.transpose();
        CoupSubmSH.block(0, n, n, n) = AB.offdiagonal.transpose();
        CoupSubmSH.block(n, 0, n, n) = AB.offdiagonal.transpose();
        CoupSubmSH = - (TmatrixSH) * CoupSubmSH;
        

        if (distance >= 4.0*(geometry.objects[ii].radius + geometry.objects[jj].radius)){ //admissibility criterion for ACA

        ACA_compression(U , V, CoupSubmSH);
        
        S_comp[nobj*ii + jj].U = U;
        S_comp[nobj*ii + jj].V = V;
        S_comp[nobj*ii + jj].dim = 2 * n;
        sizeMAT = sizeMAT + (S_comp[ii*nobj + jj].U.size())*(16.0/1e6) + (S_comp[ii*nobj + jj].V.size())*(16.0/1e6);
        }
        
        else{
        S_comp[nobj*ii + jj].S_sub =  CoupSubmSH; 
        S_comp[nobj*ii + jj].dim = 2 * n;
        sizeMAT = sizeMAT + (S_comp[ii*nobj + jj].S_sub.size())*(16.0/1e6);
        }

    }
  
      }
  
        }
        
 std::cout<<"The size of the SH matrix in MB is"  << sizeMAT<< std::endl;         
}



void ACA_compression(Matrix<t_complex> &U, Matrix<t_complex> &V, Matrix<t_complex> &CoupMat){

 int kmax = CoupMat.cols(); // square coupling matrices

 U.resize(kmax , 1);
 V.resize(1 , kmax);
 Matrix<t_complex>  sum11(1 , kmax);
 Vector<int> I(1), J(1);
 Vector<t_complex> RowCol, Row, Col, Row1, Col1, Row2, Col2, sum22;
 Vector<double> NORMA(1);
 Matrix<t_complex> R = Matrix<t_complex>::Zero(kmax, kmax);
 double eps_ACA = 1e-3; // compression tolerance
 
 for (int k = 0; k != kmax; ++k) {
 
 if(k == 0){
 I(k) = 0;
 R.row(I(k)) = CoupMat.row(I(k));
 RowCol = R.row(I(k));
 J(k) = getMaxInd(RowCol, J, kmax);
 
 V.row(k) = R.row(I(k)) / R(I(k) , J(k));
 Row = V.row(k);
 
 R.col(J(k)) = CoupMat.col(J(k));
 U.col(k) = R.col(J(k));
 Col = U.col(k);
 
 NORMA(k) = std::pow(Col.norm(), 2) * std::pow(Row.norm(), 2);
 
 RowCol = R.col(J(0));
 I.conservativeResize(2);
 I(1) = getMaxInd(RowCol, I, kmax);
 
 }
 else
 {
 
 Matrix<t_complex> sum1 = Matrix<t_complex>::Zero(1, kmax);
 
 for (int p = 0; p <= (k-1); ++p) {
 sum11 = U(I(k) , p) * V.row(p);
 
 sum1 = sum1 + sum11;
 }
 
 R.row(I(k)) = CoupMat.row(I(k)) - sum1;
 RowCol = R.row(I(k));
 J.conservativeResize(k+1);
 J(k) = getMaxInd(RowCol, J, kmax);

 V.conservativeResize(k+1 , kmax);
 V.row(k) = R.row(I(k)) / R(I(k) , J(k));
 
 Vector<t_complex> sum2 = Vector<t_complex>::Zero(kmax);
 for (int p = 0; p <= (k-1); ++p) {
 sum22 = V(p , J(k)) * U.col(p);
 sum2 = sum2 + sum22;
 }
 
 R.col(J(k)) = CoupMat.col(J(k)) - sum2;
 
 U.conservativeResize(kmax , k+1);
 U.col(k) = R.col(J(k));
 
 double sum = 0.0;
 Col2 = U.col(k);
 Row2= V.row(k);
 for (int p = 0; p != (k-1); ++p) {
 Col1 = U.col(p);
 Row1= V.row(p);

 t_complex pom1(0.0, 0.0), pom2(0.0, 0.0);
 
  for (int tt = 0; tt != k; ++tt){
  pom1 = pom1 + Col1(tt)*Col2(tt);
  pom2 = pom2 + Row1(tt)*Row2(tt);
  }

 sum = sum + abs(pom1) * abs(pom2);
 }
 
 Row = V.row(k);
 Col = U.col(k);
 
 NORMA.conservativeResize(k+1);
 NORMA(k) = NORMA(k-1) + std::pow(Col.norm(), 2) * std::pow(Row.norm(), 2) + 2.0 * sum;
 
 if (eps_ACA * sqrt(NORMA(k)) >= Col.norm() * Row.norm())
 break;
 
 RowCol = R.col(J(k));
 I.conservativeResize(k+2);
 I(k+1) = getMaxInd(RowCol, I, kmax);
 
  }
 
 }
 
}

int getMaxInd(Vector<t_complex> &RowCol, Vector<int> &K, int kmax){

double max = 0.0;
int imax, flag_same = 0;

for (int i = 0; i != kmax; ++i) {

for (int j = 0; j < (K.size()-1); ++j) {

 if(i==K(j))
 flag_same = 1;
 
 }
 
 if (flag_same==0){
 
 if (abs (RowCol(i)) > max){
   max = abs (RowCol(i));
   imax = i;
   }
 
 } 
 flag_same = 0;
}
return imax;
}

// gmres solver for ACA compressed matrices, compressed blocks are always square
Vector<t_complex> Gmres_Zcomp(std::vector<Matrix_ACA>const &S_comp, Vector<t_complex>const &Y, double tol, int maxit, int no_rest, Geometry const &geometry){

int N = Y.size(); // right hand side
int n(0), brojac(0);
mpi::Communicator communicator;
int rank = communicator.rank();

Vector<t_complex> vn(N) , w(N) , vt(N), res(N), gipom, ym, x, gi;
x = Vector<t_complex>::Zero(N);
Vector<double> err(1);
err(0) = 1;

Eigen::SparseMatrix<t_complex, Eigen::ColMajor> v(N , maxit+1);
Eigen::SparseVector<t_complex>  w_sps(N), res_sps(N);
Matrix<t_complex> H = Matrix<t_complex>::Zero(maxit + 1 , maxit);
Matrix<t_complex> Rigi, Ri, Ripom;

double beta;
double abs_y = Y.norm();

for (int rest = 1;  rest <= no_rest; ++rest) {

if(err(n)<=tol)
break;

#ifdef OPTIMET_MPI
w = matvec_parallel(S_comp , x, geometry);
#else
w = matvec(S_comp , x, geometry);
#endif

res = Y - w;
beta = res.norm();
res_sps = res.sparseView();
v.col(0) = res_sps / beta;

n = 0; 
while ((n < maxit) && (err(n) > tol)){

vn= v.col(n);

#ifdef OPTIMET_MPI
w = matvec_parallel(S_comp , vn, geometry);// matrix-vector product for compressed matrices
#else
w = matvec(S_comp , vn, geometry);
#endif

for (int t = 0; t <= n; ++t) {
vt = v.col(t);
H(t , n) = vt.adjoint() * w;
w = w - H(t , n) * vt;
}

H(n+1 , n) = w.norm();
w_sps = w.sparseView();
v.col(n+1) = w_sps / H(n+1,n);
 
Rigi.conservativeResize(n+2 , n+2);

Rigi = det_approx (beta , n , H);

Ri.conservativeResize(n+2 , n+1);

Ri = Rigi.block(0 , 0 , n+2 , n+1);

gi.conservativeResize(n+2);

gi = Rigi.col(n+1);

err.conservativeResize(n+2);
err(n+1) = abs(gi(n+1))/abs_y;
n = n + 1;
brojac++;
}

Ripom = Ri.block(0,0,n,n);
gipom = gi.segment(0 , n);
ym = Ripom.colPivHouseholderQr().solve(gipom);

for (int j = 0; j != ym.size(); ++j) {
vn =v.col(j);
x = x + ym(j) * vn;
}

}// for restart

if(rank==0){
std::cout<<"GMRES converged at iteration"<<'\t'<<brojac<<std::endl;
std::cout<<"The relative residual is"<<'\t'<<err(n)<<std::endl;
}


return x;
}

#ifdef OPTIMET_MPI
// matrix - vector product in parallel
Vector<t_complex> matvec_parallel(std::vector<Matrix_ACA>const &S_comp, Vector<t_complex> &J, Geometry const &geometry){
// here we go through the matrix blocks and check if it is compressed or not
auto const nobj = geometry.objects.size();
int gran1, gran2, Ncp_proc, N;
mpi::Communicator communicator;
  int rank = communicator.rank();
  int size = communicator.size();
 
   if(rank==0)
   N = S_comp[0].dim;
   
  MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD); 
 
    if (rank < (nobj % size)) {
    gran1 = rank * (nobj/size + 1);
    gran2 = gran1 + nobj/size + 1;
    } else {
    gran1 = rank * (nobj/size) + (nobj % size);
    gran2 = gran1 + (nobj/size);
    }

Ncp_proc = N*(gran2 - gran1);


Vector<t_complex> Y_proc = Vector<t_complex>::Zero(Ncp_proc);; // process solution
Vector<t_complex> Y_fin = Vector<t_complex>::Zero(nobj*N);; // final solution
double distance;
Matrix<t_complex> U , V;
int brojac1(0), brojac2(0);

for(int ii = gran1; ii < gran2; ii++)  {

  for(int jj = 0; jj != nobj; jj++)  {

  distance = Tools::findDistance(geometry.objects[ii].vR, geometry.objects[jj].vR);

   if (distance >= 4.0*(geometry.objects[ii].radius + geometry.objects[jj].radius)){ //admissibility criterion for ACA
      
        Y_proc.segment(brojac1*N , N) = Y_proc.segment(brojac1*N , N) + (S_comp[brojac2].U)*(S_comp[brojac2].V * J.segment(jj*N , N));
        
        }
        
    else{
    Y_proc.segment(brojac1*N , N) = Y_proc.segment(brojac1*N , N) + S_comp[brojac2].S_sub * J.segment(jj*N , N);    
     
    }
    brojac2++;
    }
    brojac1++;
  } 

// now gather all the partial solutions

   Vector<int> sizesProc(size), disps(size);

    MPI_Allgather (&Ncp_proc, 1, MPI_INT, &sizesProc(0), 1, MPI_INT, MPI_COMM_WORLD);

   for (int kk = 0; kk < size; kk++)
   disps(kk) = (kk > 0) ? (disps(kk-1) + sizesProc(kk-1)) : 0;

    MPI_Gatherv (&Y_proc(0), Ncp_proc, MPI_DOUBLE_COMPLEX, &Y_fin(0), &sizesProc(0), &disps(0), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    MPI_Bcast(&Y_fin(0), nobj*N, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

return Y_fin;
}
#endif

// matrix-vector product in serial
Vector<t_complex> matvec (std::vector<Matrix_ACA>const &S_comp, Vector<t_complex> &J, Geometry const &geometry){
auto const nobj = geometry.objects.size();
int N = S_comp[0].dim;
Vector<t_complex> Y = Vector<t_complex>::Zero(nobj*N); // solution vector 
double distance;
Matrix<t_complex> U , V;

for(int ii = 0; ii != nobj; ii++)  {

  for(int jj = 0; jj != nobj; jj++)  {

  distance = Tools::findDistance(geometry.objects[ii].vR, geometry.objects[jj].vR);

   if (distance >= 2.0*(geometry.objects[ii].radius + geometry.objects[jj].radius)){ //ACA admissibility
      
        Y.segment(ii*N , N) = Y.segment(ii*N , N) + (S_comp[ii*nobj + jj].U)*(S_comp[ii*nobj + jj].V * J.segment(jj*N , N));

        }
        
    else{
    Y.segment(ii*N , N) = Y.segment(ii*N , N) + S_comp[ii*nobj + jj].S_sub * J.segment(jj*N , N);           
    }
    }
   
}

return Y;
}

Matrix<t_complex> det_approx (double beta, int n, Matrix<t_complex> &H){

Matrix<t_complex> Rigi, Ri, W, POM(2,2);

Rigi.conservativeResize(n+2 , n+2);
Ri.conservativeResize(n+2 , n+1);


Vector<t_complex> gi = Vector<t_complex>::Zero(n+2);
t_complex hi1, hi2, temp, c, s;

Ri = H.block(0,0, n+2,n+1);
gi(0) = beta;

for (int i = 0; i <= n; ++i) {

hi1 = Ri(i,i);
hi2 = Ri(i+1,i);

if (abs(hi2)>abs(hi1)){
      temp = hi1/hi2; 
      s = 1.0 / sqrt(1.0 + std::pow(abs(temp),2)); 
      c = -temp*s;
      }
      
 else{
     temp = hi2/hi1; 
      c = 1.0 / sqrt(1.0 + std::pow(abs(temp),2)); 
      s = -temp*c;
 } 

 POM(0,0) = std::conj(c);
 POM(0,1) = std::conj(-s);
 POM(1,0) = s;
 POM(1,1) = c;
 W = Matrix<t_complex>::Identity(n+2 , n+2);
 W.block(i,i, 2, 2)  = POM;
 
 Ri = W*Ri;
 gi = W*gi; 
}

Rigi.block(0,0, n+2,n+1) = Ri;
Rigi.col(n+1) = gi; 
 
return Rigi;
}



Matrix<t_complex> preconditioned_scattering_matrix(std::vector<Scatterer> const &objects,
                                                   ElectroMagnetic const &bground,
                                                   std::shared_ptr<Excitation const> incWave) {
                                                                                                  
  return preconditioned_scattering_matrix(objects.begin(), objects.end(), objects.begin(),
                                          objects.end(), bground, incWave);
}

 Matrix<t_complex> preconditioned_scattering_matrixSH(std::vector<Scatterer> const &objects,
                                                   ElectroMagnetic const &bground,
                                                   std::shared_ptr<Excitation const> incWave) {
                                                                                                    
  return preconditioned_scattering_matrixSH(objects.begin(), objects.end(), objects.begin(),
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


Matrix<t_complex> preconditioned_scattering_matrixSH(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave) {
                                                   
                                                  
  if(geometry.objects.size() == 0)
    return Matrix<t_complex>(0, 0);
  // Check nMax is same accross all objects
    auto const nMaxS = geometry.objects.front().nMaxS;
    for(auto const &scatterer : geometry.objects)
     if(scatterer.nMaxS != nMaxS)
     throw std::runtime_error("All objects must have same number of SH harmonics"); 
     return preconditioned_scattering_matrixSH(geometry.objects, geometry.bground, incWave);
                      
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

Matrix<t_complex> preconditioned_scattering_matrix_SH(Geometry const &geometry,
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

  auto const nMaxS = geometry.objects.front().nMaxS;
  auto const remainder = linear_context.is_valid() ? nobj % linear_context.size() : 0;
  auto const nloc = linear_context.is_valid() ? nobj / linear_context.size() : 0;
  t_uint n;

  n = nMaxS * (nMaxS + 2);

  scalapack::Sizes const non_cyclic{linear_context.is_valid() ? nobj * n * 2 : 1,
                                    linear_context.is_valid() ? nloc * n * 2 : 1};
  scalapack::Matrix<t_complex> linear_matrix(linear_context, {nobj * n * 2, nobj * n * 2},
                                             non_cyclic);

  if(linear_context.is_valid()) {
    assert(geometry.objects.end() > geometry.objects.begin() + nloc * linear_context.col());
    assert(geometry.objects.end() >= geometry.objects.begin() + nloc * (1 + linear_context.col()));
    assert(nloc * 2 * n > 0);
    linear_matrix.local().leftCols(nloc * 2 * n) = preconditioned_scattering_matrixSH(
        geometry.objects.begin(), geometry.objects.end(),
        geometry.objects.begin() + nloc * linear_context.col(),
        geometry.objects.begin() + nloc * (1 + linear_context.col()), geometry.bground, incWave);
  }

  if(remainder > 0 and linear_context.is_valid()) {
    auto const remainder_context = linear_context.subcontext(rank_map.leftCols(remainder));
    scalapack::Matrix<t_complex> remainder_matrix(
        remainder_context, {nobj * n * 2, remainder * n * 2}, {nobj * n * 2, 2 * n});
    if(remainder_context.is_valid())
      remainder_matrix.local() = preconditioned_scattering_matrixSH(
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

Matrix<t_complex> preconditioned_scattering_matrix_SH(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave,
                                                   scalapack::Context const &,
                                                   scalapack::Sizes const &) {
  return preconditioned_scattering_matrixSH(geometry, incWave);


}
#endif

Vector<t_complex> source_vector(std::vector<Scatterer>::const_iterator first,
                                std::vector<Scatterer>::const_iterator const &last,
                                std::shared_ptr<Excitation const> incWave, Geometry const &geometry) {
  if(first == last)
    return Vector<t_complex>::Zero(0);
  auto const nMax = first->nMax;
  auto const flatMax = nMax * (nMax + 2);
  Matrix<t_complex> TmatrixFF (2 * flatMax , 2 * flatMax);
  Vector<t_complex> result(2 * flatMax * (last - first));

  for(size_t i(0); first != last; ++first, i += 2 * flatMax){

    incWave->getIncLocal(first->vR, result.data() + i, nMax);
    first->getTLocal(TmatrixFF, incWave->omega(), geometry.bground);   
    result.segment(i , 2 * flatMax) = TmatrixFF * result.segment(i , 2 * flatMax);

    }
  return result;
}

Vector<t_complex> source_vectorSH(Geometry &geometry, std::vector<Scatterer>::const_iterator first,
                                std::vector<Scatterer>::const_iterator const &last,
                                std::shared_ptr<Excitation const> incWave, Vector<t_complex> &internalCoef_FF_, 
                                  Vector<t_complex> &scatteredCoef_FF_, std::vector<double *> CGcoeff) {
                                
  if(first == last)
  return Vector<t_complex>::Zero(0);

  auto const nobj = geometry.objects.size();
  auto const nMaxS = first->nMaxS;
  auto const flatMax = nMaxS * (nMaxS + 2);
  int objectIndex_=0;
  auto const k_b_SH = 2 * incWave->omega() * std::sqrt(geometry.bground.epsilon * geometry.bground.mu);
  Vector<t_complex> result, result1, result3, resultAna;
  Matrix<t_complex> TmatrixSH (2*flatMax, 2*flatMax);
  
  if (first->scatterer_type == "sphere"){

  result.resize((2 * flatMax * (last - first)));
  resultAna.resize((4 * flatMax * (last - first)));

  for(size_t i(0); first != last; ++first, i += 4 * flatMax){
  
  geometry.getIncLocalSH(CGcoeff, objectIndex_, incWave, internalCoef_FF_, nMaxS, resultAna.data()+i);
  
  objectIndex_++;
    
    }

   size_t i(0);
   size_t ii(0);
 
   for(int kk = 0; kk != nobj; kk++)  {
 
    result.segment(ii, 2*flatMax) =
        (geometry.objects[kk].getTLocalSH1_outer(incWave->omega(), geometry.bground).array() * resultAna.segment(i, 2*flatMax).array()) +
        (geometry.objects[kk].getTLocalSH2_outer(incWave->omega(), geometry.bground).array() * resultAna.segment(i + 2*flatMax, 2*flatMax).array());
       
  i += 4 * flatMax;
  ii += 2 * flatMax;

  }
    }//if sphere

  return result;
  
}

Vector<t_complex> source_vectorSH_K1ana(Geometry &geometry,std::shared_ptr<Excitation const> incWave, 
Vector<t_complex> &internalCoef_FF_, Vector<t_complex> &scatteredCoef_FF_, std::vector<double *> CGcoeff) {
                                
  
  auto const nobj = geometry.objects.size();
  auto const nMaxS = geometry.objects[0].nMaxS;
  auto const flatMax = nMaxS * (nMaxS + 2);
  Vector<t_complex> result, resultAna;

  if (geometry.objects[0].scatterer_type == "sphere"){
  
  result.resize((2 * flatMax * nobj));
  resultAna.resize((4 * flatMax * nobj));
  
  int i(0);
  int ii(0);
  
  for(int objIndex_ = 0; objIndex_ != nobj; objIndex_++)  {
  
  geometry.getIncLocalSH(CGcoeff, objIndex_, incWave, internalCoef_FF_, nMaxS, resultAna.data()+i);
  
  i += 4 * flatMax;
    
    }
    
  i = 0;
 // here I need only vpp and upp
 for(int kk = 0; kk != nobj; kk++)  {
 
    result.segment(ii, 2 *flatMax) =
    
    geometry.objects[kk].getIauxSH2(incWave->omega(), geometry.bground).array() * resultAna.segment(i + 2*flatMax, 2 * flatMax).array();
  
  i += 4 * flatMax;
  ii += 2 * flatMax;
  }
    
    }
    
    return result;
  
}

#ifdef OPTIMET_MPI
Vector<t_complex> source_vectorSH_parallel(Geometry &geometry, int gran1, int gran2,
                                std::shared_ptr<Excitation const> incWave, Vector<t_complex> &internalCoef_FF_, std::vector<double *> CGcoeff) {
  
                              
  if(gran1 == gran2)
  return Vector<t_complex>::Zero(0);

  auto const nMaxS = geometry.objects.front().nMaxS;
  
  Vector<t_complex> resultProc(4*(gran2 - gran1));
  
  geometry.getIncLocalSH_parallel(CGcoeff, gran1, gran2, incWave, internalCoef_FF_, nMaxS, resultProc.data());
    

  return resultProc;
  
}
#endif

Vector<t_complex> source_vector(std::vector<Scatterer> const &objects, std::shared_ptr<Excitation const> incWave, Geometry const &geometry) {
  return source_vector(objects.begin(), objects.end(), incWave, geometry);
}



Vector<t_complex> source_vector(Geometry const &geometry, std::shared_ptr<Excitation const> incWave) {
  if(geometry.objects.size() == 0)
    return Vector<t_complex>(0, 0);
  // Check nMax is same accross all objects
  auto const nMax = geometry.objects.front().nMax;
  for(auto const &scatterer : geometry.objects)
    if(scatterer.nMax != nMax)
      throw std::runtime_error("All objects must have same number of harmonics");
  return source_vector(geometry.objects, incWave, geometry);
}




Vector<t_complex> source_vectorSH(Geometry &geometry, std::vector<Scatterer> const &objects, std::shared_ptr<Excitation const> incWave, 
Vector<t_complex> &internalCoef_FF_, Vector<t_complex> &scatteredCoef_FF_, std::vector<double *> CGcoeff) {
  return source_vectorSH(geometry, objects.begin(), objects.end(), incWave, internalCoef_FF_, scatteredCoef_FF_, CGcoeff);
}




Vector<t_complex> source_vectorSH(Geometry &geometry, std::shared_ptr<Excitation const> incWave, Vector<t_complex> &internalCoef_FF_, Vector<t_complex> &scatteredCoef_FF_, std::vector<double *> CGcoeff) {
  if(geometry.objects.size() == 0)
    return Vector<t_complex>(0, 0);
  // Check nMax is same accross all objects
  auto const nMaxS = geometry.objects.front().nMaxS;
  for(auto const &scatterer : geometry.objects)
    if(scatterer.nMaxS != nMaxS)
      throw std::runtime_error("All objects must have same number of harmonics");
  return source_vectorSH(geometry, geometry.objects, incWave, internalCoef_FF_, scatteredCoef_FF_, CGcoeff);
}

}
