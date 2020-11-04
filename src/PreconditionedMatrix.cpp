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

Vector<t_complex> distributed_source_vector_SH_Mnode(Geometry &geometry,
                                           std::shared_ptr<Excitation const> incWave,
                                           Vector<t_complex> &X_int_, std::vector<double *> CGcoeff) {
auto const nobj = geometry.objects.size();
  if(nobj == 0)
     return Vector<t_complex>::Zero(0);

  int gran, gran1, gran2, objInd, kk;
  mpi::Communicator communicator;
  int rank = communicator.rank();
  int size = communicator.size();
  auto const nMaxS = geometry.objects.front().nMaxS;
  t_uint const pMax = nMaxS * (nMaxS + 2);
  int TMax = nobj * pMax;
  Vector<t_complex> result (4 * nobj * pMax);
  Vector<t_complex> resultK (4 * nobj * pMax);
  
  Vector<t_complex> X_int_proc;
  int sizeFFint;
    
    if (rank==0){
    sizeFFint =X_int_.size();
    X_int_proc = X_int_;
    }
    MPI_Bcast(&sizeFFint, 1, MPI_INT, 0, MPI_COMM_WORLD);
    X_int_proc.resize(sizeFFint);
    MPI_Bcast(&X_int_proc(0), sizeFFint, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  

      gran1 = (TMax / (size))*(rank);

      if(rank<(size-1)) {

      gran2 = (rank+1)*(TMax / (size));}

       	else { gran2 = TMax;}

    int sizeVec = 4 * (gran2 - gran1);
          
    Vector<t_complex> resultProc(sizeVec);

    Vector<int> sizesProc(size), disps(size);
  
    MPI_Allgather (&sizeVec, 1, MPI_INT, &sizesProc(0), 1, MPI_INT, MPI_COMM_WORLD);
  
   for (int kk = 0; kk < size; kk++)
   disps(kk) = (kk > 0) ? (disps(kk-1) + sizesProc(kk-1)) : 0;
      
    resultProc = source_vectorSH_parallel(geometry, gran1, gran2, incWave, X_int_proc, CGcoeff);

    MPI_Gatherv (&resultProc(0), sizeVec, MPI_DOUBLE_COMPLEX, &result(0), &sizesProc(0), &disps(0), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD); 


  if (rank == 0){ 
  
  for (int ranki = 0; ranki < size; ranki++){

     gran1 = (TMax / (size))*(ranki);

      if(ranki<(size-1)) {

      gran2 = (ranki+1)*(TMax / (size));}

        else { gran2 = TMax;}

      sizeVec = 4 * (gran2 - gran1);

  int brojac (0);
  for (int ii = gran1; ii < gran2; ii++){

  objInd = ii / pMax;

  kk = ii - objInd * pMax; 

  gran = objInd * 4 * pMax + kk;

  resultK(gran) = result(disps(ranki) + brojac);
  resultK(gran + pMax) = result(disps(ranki) + (sizeVec/4) + brojac);
  resultK(gran + 2*pMax) = result(disps(ranki) + 2*(sizeVec/4) + brojac);
  resultK(gran + 3*pMax) = result(disps(ranki) + 3*(sizeVec/4) + brojac);

  brojac++; 
  
  } // for ii

  }  // for ranki
 
  } // if

return resultK;
}


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


Matrix<t_complex>
preconditioned_scattering_matrixSH(std::vector<Scatterer>::const_iterator const &first,
                                 std::vector<Scatterer>::const_iterator const &end_first,
                                 std::vector<Scatterer>::const_iterator const &second,
                                 std::vector<Scatterer>::const_iterator const &end_second,
                                 ElectroMagnetic const &bground,
                                 std::shared_ptr<Excitation const> incWave) {
                                 
                              
  auto const nMaxS = first->nMaxS;
  auto const n = nMaxS * (nMaxS + 2);
  

  if(first == end_first or second == end_second)
  return Matrix<t_complex>::Zero(4 * n * (end_first - first), 4 * n * (end_second - second));

 Matrix<t_complex> resultSH = Matrix<t_complex>::Zero(4 * n * (end_first - first), 4 * n * (end_second - second));
        
  size_t y(0);

  for(auto iterj(second); iterj != end_second; ++iterj, y += 4 * n) {

    Vector<t_complex> const factorSH1 = -iterj->getTLocalSH1_outer(incWave->omega(), bground);  // bp and ap coefficients
    Vector<t_complex> const factorSH2 = -iterj->getTLocalSH2_outer(incWave->omega(), bground);  // bpp and app coefficients
 
    size_t x(0);
    
    for(auto iteri(first); iteri != end_first; ++iteri, x += 4 * n) {

      if(iteri == iterj) {

        resultSH.block(x, y, 4 * n, 4 * n) = Matrix<t_complex>::Identity(4 * n, 4 * n);
  
      } 
      
      else {


        Coupling const AB(iteri->vR - iterj->vR, 2.0 * incWave->waveK, nMaxS);
 
        resultSH.block(x, y, n, n) = AB.diagonal.transpose();
        resultSH.block(x + n, y + n, n, n) = AB.diagonal.transpose();
        resultSH.block(x, y + n, n, n) = AB.offdiagonal.transpose();
        resultSH.block(x + n, y, n, n) = AB.offdiagonal.transpose();
        resultSH.block(x, y, 2 * n, 2 * n).array().transpose().colwise() *= factorSH1.array();
        
        
        resultSH.block(x + 2*n, y + 2*n, n, n) = AB.diagonal.transpose();
        resultSH.block(x + 2*n + n, y + 2*n + n, n, n) = AB.diagonal.transpose();
        resultSH.block(x + 2*n, y + 2*n + n, n, n) = AB.offdiagonal.transpose();
        resultSH.block(x + 2*n + n, y + 2*n, n, n) = AB.offdiagonal.transpose();
        resultSH.block(x + 2*n, y + 2*n, 2 * n, 2 * n).array().transpose().colwise() *= factorSH2.array();
        
        
     }
     
   }
   
 }
  
    
  return resultSH;
  
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
  t_uint const n = 2 * nMaxS * (nMaxS + 2);
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

Matrix<t_complex> preconditioned_scattering_matrix_SH(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave,
                                                   scalapack::Context const &,
                                                   scalapack::Sizes const &) {
  return preconditioned_scattering_matrixSH(geometry, incWave);


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
  for(size_t i(0); first != last; ++first, i += 2 * flatMax){
    incWave->getIncLocal(first->vR, result.data() + i, nMax);
    }
  return result;
}

Vector<t_complex> source_vectorSH(Geometry &geometry, std::vector<Scatterer>::const_iterator first,
                                std::vector<Scatterer>::const_iterator const &last,
                                std::shared_ptr<Excitation const> incWave, Vector<t_complex> &internalCoef_FF_, std::vector<double *> CGcoeff) {
                                
  if(first == last)
  return Vector<t_complex>::Zero(0);
  auto const nMaxS = first->nMaxS;
  auto const flatMax = nMaxS * (nMaxS + 2);
  int objectIndex_=0;
  
  Vector<t_complex> result(4 * flatMax * (last - first));
  
  for(size_t i(0); first != last; ++first, i += 4 * flatMax){
  
  
  geometry.getIncLocalSH(CGcoeff, objectIndex_, incWave, internalCoef_FF_, nMaxS, result.data()+i);
  
  objectIndex_++;
    
    }

  return result;
  
}


Vector<t_complex> source_vectorSH_parallel(Geometry &geometry, int gran1, int gran2,
                                std::shared_ptr<Excitation const> incWave, Vector<t_complex> &internalCoef_FF_, std::vector<double *> CGcoeff) {
  
                              
  if(gran1 == gran2)
  return Vector<t_complex>::Zero(0);

  auto const nMaxS = geometry.objects.front().nMaxS;
  
  Vector<t_complex> resultProc(4*(gran2 - gran1));
  
  geometry.getIncLocalSH_parallel(CGcoeff, gran1, gran2, incWave, internalCoef_FF_, nMaxS, resultProc.data());
    

  return resultProc;
  
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




Vector<t_complex>
source_vectorSH(Geometry &geometry, std::vector<Scatterer> const &objects,
                std::shared_ptr<Excitation const> incWave, Vector<t_complex> &internalCoef_FF_, std::vector<double *> CGcoeff) {
  return source_vectorSH(geometry, objects.begin(), objects.end(), incWave, internalCoef_FF_, CGcoeff);
}




Vector<t_complex>
source_vectorSH(Geometry &geometry, std::shared_ptr<Excitation const> incWave, Vector<t_complex> &internalCoef_FF_, std::vector<double *> CGcoeff) {
  if(geometry.objects.size() == 0)
    return Vector<t_complex>(0, 0);
  // Check nMax is same accross all objects
  auto const nMaxS = geometry.objects.front().nMaxS;
  for(auto const &scatterer : geometry.objects)
    if(scatterer.nMaxS != nMaxS)
      throw std::runtime_error("All objects must have same number of harmonics");
  return source_vectorSH(geometry, geometry.objects, incWave, internalCoef_FF_, CGcoeff);
}



}



