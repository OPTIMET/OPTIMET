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
using namespace std::chrono;

namespace optimet {
#ifdef OPTIMET_MPI
Vector<t_complex> distributed_vector_SH_AR1(Geometry &geometry,
                                           std::shared_ptr<Excitation const> incWave,
                                           Vector<t_complex> &X_int_, Vector<t_complex> &X_sca_) {
auto const nobj = geometry.objects.size();
  if(nobj == 0)
     return Vector<t_complex>::Zero(0); 
  int gran1, gran2, gran11, gran22, sizeVec2;
  mpi::Communicator communicator;
  int rank = communicator.rank();
  int size = communicator.size();
  auto const nMaxS = geometry.objects.front().nMaxS;

  t_uint const pMax = nMaxS * (nMaxS + 2);

  Vector<t_complex> result1(2*pMax), resultK1(2*pMax), resultK1MNY(2*nobj*pMax);
 
  Vector<t_complex> X_int_proc, X_sca_proc;

  int sizeFFint;
    
    if (rank==0){
    sizeFFint =X_int_.size();
    X_int_proc = X_int_;
    X_sca_proc = X_sca_;
    }

    MPI_Bcast(&sizeFFint, 1, MPI_INT, 0, MPI_COMM_WORLD);
    X_int_proc.resize(sizeFFint); // broadcasting internal and external FF field coeff
    X_sca_proc.resize(sizeFFint);
    MPI_Bcast(&X_int_proc(0), sizeFFint, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&X_sca_proc(0), sizeFFint, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    if (rank < (pMax % size)) {
    gran1 = rank * (pMax/size + 1);
    gran2 = gran1 + pMax/size + 1;
    } else {
    gran1 = rank * (pMax/size) + (pMax % size);
    gran2 = gran1 + (pMax/size);
    }


  if (geometry.objects[0].scatterer_type == "arbitrary.shape"){

   for (int objIndex = 0; objIndex < nobj; objIndex++){

    int sizeVec = 2 * (gran2 - gran1);

     resultK1.setZero();

    Vector<t_complex> resultProc1(sizeVec);

    Vector<int> sizesProc(size), disps(size);

   MPI_Allgather (&sizeVec, 1, MPI_INT, &sizesProc(0), 1, MPI_INT, MPI_COMM_WORLD);

   for (int kk = 0; kk < size; kk++)
   disps(kk) = (kk > 0) ? (disps(kk-1) + sizesProc(kk-1)) : 0;

    
    resultProc1 = source_vectorSH_parallelAR1(geometry, gran1, gran2, incWave, X_int_proc, X_sca_proc, objIndex);

    MPI_Gatherv (&resultProc1(0), sizeVec, MPI_DOUBLE_COMPLEX, &result1(0), &sizesProc(0), &disps(0), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    if (rank == 0){

    for (int ranki = 0; ranki < size; ranki++){

     if (ranki < (pMax % size)) {
    gran11 = ranki * (pMax/size + 1);
    gran22 = gran11 + pMax/size + 1;
    } else {
    gran11 = ranki * (pMax/size) + (pMax % size);
    gran22 = gran11 + (pMax/size);
    }


      sizeVec2 = 2 * (gran22 - gran11);

    int brojac (0);

    for (int ii = gran11; ii < gran22; ii++){
 
   resultK1(ii) = result1(disps(ranki) + brojac);
   resultK1(ii + pMax) = result1(disps(ranki) + (sizeVec2/2) + brojac);
 
   brojac++;
 
   } // for ii
 
   }  // for ranki
  
   resultK1MNY.segment(objIndex*2*pMax, 2*pMax) = resultK1;
   }// rank0

     MPI_Barrier(MPI_COMM_WORLD);
 
    }// loop through the objects

      MPI_Bcast(&resultK1MNY(0), 2*nobj*pMax, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

     }// end if arb.shapes
 
       return resultK1MNY;
 
       }
 
Vector<t_complex> distributed_source_vector_SH_Mnode(Geometry &geometry,
                                           std::shared_ptr<Excitation const> incWave,
                                           Vector<t_complex> &X_int_, Vector<t_complex> &X_sca_, 
                                           Matrix<t_complex> &TmatrixSH) {
auto const nobj = geometry.objects.size();
  if(nobj == 0)
     return Vector<t_complex>::Zero(0);

  int gran1, gran2, gran11, gran22, kk, sizeVec2;
  mpi::Communicator communicator;
  int rank = communicator.rank();
  int size = communicator.size();
  auto const nMaxS = geometry.objects.front().nMaxS;
  auto const k_b_SH = 2.0 * incWave->omega() * std::sqrt(geometry.bground.epsilon * geometry.bground.mu);

  t_uint const pMax = nMaxS * (nMaxS + 2);

  Vector<t_complex> resultK(2*pMax), result3(2*pMax), result1(2*pMax), resultK3(2*pMax), resultK1(2*pMax), resultKMNY(2*nobj*pMax);
  Vector<t_complex> X_int_proc, X_sca_proc;
  Matrix<t_complex> TmatrixSINGLE;

  int sizeFFint;
    
    if (rank==0){
    sizeFFint =X_int_.size();
    X_int_proc = X_int_;
    X_sca_proc = X_sca_;
    }
    
    MPI_Bcast(&sizeFFint, 1, MPI_INT, 0, MPI_COMM_WORLD);
    X_int_proc.resize(sizeFFint); // broadcasting internal and external FF field coeff
    X_sca_proc.resize(sizeFFint);
    MPI_Bcast(&X_int_proc(0), sizeFFint, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&X_sca_proc(0), sizeFFint, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    
    if (rank < (pMax % size)) {
    gran1 = rank * (pMax/size + 1);
    gran2 = gran1 + pMax/size + 1;
    } else {
    gran1 = rank * (pMax/size) + (pMax % size);
    gran2 = gran1 + (pMax/size);
    }
    
  
// arbitrary shapes
 if (geometry.objects[0].scatterer_type == "arbitrary.shape"){

  for (int objIndex = 0; objIndex < nobj; objIndex++){
    int sizeVec = 2 * (gran2 - gran1);
     TmatrixSINGLE = TmatrixSH.block(0, objIndex*2*pMax, 2*pMax, 2*pMax);
     resultK.setZero();
     resultK1.setZero();
     resultK3.setZero();

    Vector<t_complex> resultProc3(sizeVec), resultProc1(sizeVec);

    Vector<int> sizesProc(size), disps(size);
   
   MPI_Allgather (&sizeVec, 1, MPI_INT, &sizesProc(0), 1, MPI_INT, MPI_COMM_WORLD);
   
   for (int kk = 0; kk < size; kk++)
   disps(kk) = (kk > 0) ? (disps(kk-1) + sizesProc(kk-1)) : 0;
 
    resultProc3 = source_vectorSH_parallelAR3(geometry, gran1, gran2, incWave, X_int_proc, X_sca_proc, objIndex);
    resultProc1 = source_vectorSH_parallelAR1(geometry, gran1, gran2, incWave, X_int_proc, X_sca_proc, objIndex); 
     
    MPI_Gatherv (&resultProc3(0), sizeVec, MPI_DOUBLE_COMPLEX, &result3(0), &sizesProc(0), &disps(0), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Gatherv (&resultProc1(0), sizeVec, MPI_DOUBLE_COMPLEX, &result1(0), &sizesProc(0), &disps(0), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
      
    if (rank == 0){

    for (int ranki = 0; ranki < size; ranki++){

     if (ranki < (pMax % size)) {
    gran11 = ranki * (pMax/size + 1);
    gran22 = gran11 + pMax/size + 1;
    } else {
    gran11 = ranki * (pMax/size) + (pMax % size);
    gran22 = gran11 + (pMax/size);
    }

      sizeVec2 = 2 * (gran22 - gran11);

  int brojac (0);
  for (int ii = gran11; ii < gran22; ii++){

  resultK3(ii) = result3(disps(ranki) + brojac);
  resultK3(ii + pMax) = result3(disps(ranki) + (sizeVec2/2) + brojac);
  resultK1(ii) = result1(disps(ranki) + brojac);
  resultK1(ii + pMax) = result1(disps(ranki) + (sizeVec2/2) + brojac);

  brojac++;

  } // for ii

  }  // for ranki
        
  resultK = (-consCi * k_b_SH ) * resultK1 + (consCi * k_b_SH) * TmatrixSINGLE * resultK3;
  resultKMNY.segment(objIndex*2*pMax, 2*pMax) = resultK;
  }// rank0

  MPI_Barrier(MPI_COMM_WORLD);


  } // end loop on objects

 MPI_Bcast(&resultKMNY(0), 2*nobj*pMax, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
}// end if arb.shapes

return resultKMNY;

}


 Vector<t_complex> getQmatrix_FF(Geometry const &geometry, ElectroMagnetic const &bground,
                                 std::shared_ptr<Excitation const> incWave, int gran1, int gran2, int objIndex) {
  
  auto const nMax = geometry.objects.front().nMax;                    
  auto const n = nMax * (nMax + 2);
  int numRow = gran2 - gran1; 
  
  Vector<t_complex> Qmatrix (4 * n * numRow);   

  geometry.objects[objIndex].getQLocal(Qmatrix, incWave->omega(), bground, gran1, gran2);       
  
  return Qmatrix;
}

 Vector<t_complex> getRgQmatrix_FF(Geometry const &geometry, ElectroMagnetic const &bground,
                                 std::shared_ptr<Excitation const> incWave, int gran1, int gran2, int objIndex) {

  auto const nMax = geometry.objects.front().nMax;
  auto const n = nMax * (nMax + 2);
  int numRow = gran2 - gran1;

  Vector<t_complex> RgQmatrix (4 * n * numRow);

  geometry.objects[objIndex].getRgQLocal(RgQmatrix, incWave->omega(), bground, gran1, gran2);


  return RgQmatrix;
}


 Vector<t_complex> getQmatrix_SH(Geometry const &geometry,
                                 ElectroMagnetic const &bground,
                                 std::shared_ptr<Excitation const> incWave, int gran1, int gran2, int objIndex) {
                                 
  auto const nMaxS = geometry.objects.front().nMaxS;        
  auto const n = nMaxS * (nMaxS + 2);
  int numRow = gran2 - gran1;

  Vector<t_complex> QmatrixSH (4 * n * numRow);

  geometry.objects[objIndex].getQLocalSH(QmatrixSH, incWave->omega(), bground, gran1, gran2);
     
  return QmatrixSH;
  
}

 Vector<t_complex> getRgQmatrix_SH(Geometry const &geometry,
                                 ElectroMagnetic const &bground,
                                 std::shared_ptr<Excitation const> incWave, int gran1, int gran2, int objIndex) {

  auto const nMaxS = geometry.objects.front().nMaxS;
  auto const n = nMaxS * (nMaxS + 2);
  int numRow = gran2 - gran1;

  Vector<t_complex> RgQmatrixSH (4 * n * numRow);

  geometry.objects[objIndex].getRgQLocalSH(RgQmatrixSH, incWave->omega(), bground, gran1, gran2);
 
    
  return RgQmatrixSH;
  
}


//Preconditioned scattering matrix for FF
 Matrix<t_complex> ScatteringMatrixFF(Matrix<t_complex> &TMatrixFF, Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave) {

  auto const nMax = geometry.objects.front().nMax;
  auto const n = nMax * (nMax + 2);
  auto const nobj = geometry.objects.size();
  Matrix<t_complex> TmatrixSingle (2*n , 2*n);  

  Matrix<t_complex> result(2 * n * nobj, 2 * n * nobj);

  size_t y(0);
  for(int jj = 0; jj != nobj; ++jj, y += 2 * n) {

  TmatrixSingle = TMatrixFF.block(0, jj*2*n, 2 * n, 2 * n);

  size_t x(0);
  for(int ii = 0; ii != nobj; ++ii, x += 2 * n) {


      if(ii == jj) {
        result.block(x, y, 2 * n, 2 * n) = Matrix<t_complex>::Identity(2 * n, 2 * n);


      } else {

        Coupling const AB(geometry.objects[ii].vR - geometry.objects[jj].vR, 1.0 * incWave->waveK, nMax);

        result.block(x, y, n, n) = AB.diagonal.transpose();
        result.block(x + n, y + n, n, n) = AB.diagonal.transpose();
        result.block(x, y + n, n, n) = AB.offdiagonal.transpose();
        result.block(x + n, y, n, n) = AB.offdiagonal.transpose();
        result.block(x, y, 2 * n, 2 * n) = - result.block(x, y, 2 * n, 2 * n) * (TmatrixSingle);


  }

    }

     }
return result;   
       }

// Scattering matrix for SH
 Matrix<t_complex> ScatteringMatrixSH(Matrix<t_complex> &TMatrixSH, Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave) {

  auto const nMaxS = geometry.objects.front().nMaxS;
  auto const n = nMaxS * (nMaxS + 2);
  auto const nobj = geometry.objects.size();
  auto const k_b_SH = 2 * incWave->omega() * std::sqrt(geometry.bground.epsilon * geometry.bground.mu);
  Matrix<t_complex> TmatrixSingle (2*n , 2*n);

  Matrix<t_complex> result(2 * n * nobj, 2 * n * nobj);

  size_t x(0);
  for(int ii = 0; ii != nobj; ++ii, x += 2 * n) {

  TmatrixSingle = TMatrixSH.block(0, ii*2*n, 2 * n, 2 * n);

  size_t y(0);
  for(int jj = 0; jj != nobj; ++jj, y += 2 * n) {


      if(ii == jj) {
        result.block(x, y, 2 * n, 2 * n) = Matrix<t_complex>::Identity(2 * n, 2 * n);


      } else {

        Coupling const AB(geometry.objects[ii].vR - geometry.objects[jj].vR, 2.0 * incWave->waveK, nMaxS);

        result.block(x, y, n, n) = AB.diagonal.transpose();
        result.block(x + n, y + n, n, n) = AB.diagonal.transpose();
        result.block(x, y + n, n, n) = AB.offdiagonal.transpose();
        result.block(x + n, y, n, n) = AB.offdiagonal.transpose();
        result.block(x, y, 2 * n, 2 * n) = (TmatrixSingle) * result.block(x, y, 2 * n, 2 * n);


    }
}

}
return result;
}



Matrix<t_complex> getTRgQmatrix_FF_parr(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave
                                                    ) {
 auto const nobj = geometry.objects.size();
 int gran1, gran2, gran11, gran22, rank, size, pMax, sizeVec2;
 mpi::Communicator communicator;
 rank = communicator.rank();
 size = communicator.size();
 auto const nMax = geometry.objects.front().nMax;

   pMax = nMax * (nMax + 2);

   if (rank < (pMax % size)) {
    gran1 = rank * (pMax/size + 1);
    gran2 = gran1 + pMax/size + 1;
    } else {
    gran1 = rank * (pMax/size) + (pMax % size);
    gran2 = gran1 + (pMax/size);
    }

 int numRow = gran2 - gran1;

 int sizeVec = 4 * numRow * pMax;

 Vector<t_complex> QmatrixFF_proc (sizeVec);
 Vector<t_complex> RgQmatrixFF_proc (sizeVec);
 Vector<t_complex> resultQ (4 * pMax * pMax);
 Vector<t_complex> resultRgQ (4 * pMax * pMax);
 Matrix<t_complex> QmatrixFF(2*pMax, 2*pMax), RgQmatrixFF(2*pMax, 2*pMax), TmatrixFF(2*pMax, 2*pMax), TRgQmatrixFF(2*pMax, 4*nobj*pMax);
 TRgQmatrixFF.setZero();

for (int objIndex = 0; objIndex < nobj; objIndex++){

  resultQ.setZero();
  resultRgQ.setZero();

  QmatrixFF_proc = getQmatrix_FF(geometry, geometry.bground, incWave, gran1, gran2, objIndex);
  RgQmatrixFF_proc = getRgQmatrix_FF(geometry, geometry.bground, incWave, gran1, gran2, objIndex);

    Vector<int> sizesProc(size), disps(size);

   MPI_Allgather (&sizeVec, 1, MPI_INT, &sizesProc(0), 1, MPI_INT, MPI_COMM_WORLD);

   for (int kk = 0; kk < size; kk++)
   disps(kk) = (kk > 0) ? (disps(kk-1) + sizesProc(kk-1)) : 0; // displacements

  MPI_Gatherv (&QmatrixFF_proc(0), sizeVec, MPI_DOUBLE_COMPLEX, &resultQ(0), &sizesProc(0), &disps(0), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Gatherv (&RgQmatrixFF_proc(0), sizeVec, MPI_DOUBLE_COMPLEX, &resultRgQ(0), &sizesProc(0), &disps(0), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
 
  MPI_Bcast(&resultQ(0), 4*pMax*pMax, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(&resultRgQ(0), 4*pMax*pMax, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  
 // rearranging vectors to matrices 

    for (int ranki = 0; ranki < size; ranki++){

    if (ranki < (pMax % size)) {
    gran11 = ranki * (pMax/size + 1);
    gran22 = gran11 + pMax/size + 1;
    } else {
    gran11 = ranki * (pMax/size) + (pMax % size);
    gran22 = gran11 + (pMax/size);
    }

    sizeVec2 = pMax * (gran22 - gran11);

  int brojac (0);
 
   for (int ii = gran11; ii < gran22; ii++){
    for (int jj = 0; jj < pMax; jj++){
 
  QmatrixFF(ii, jj) = resultQ(disps(ranki) + brojac);
  QmatrixFF(ii, jj + pMax) = resultQ(disps(ranki) + sizeVec2 + brojac);
  QmatrixFF(ii + pMax, jj) = resultQ(disps(ranki) + 2*sizeVec2 + brojac);
  QmatrixFF(ii + pMax, jj + pMax) = resultQ(disps(ranki) + 3*sizeVec2 + brojac);

  RgQmatrixFF(ii, jj) = resultRgQ(disps(ranki) + brojac);
  RgQmatrixFF(ii, jj + pMax) = resultRgQ(disps(ranki) + sizeVec2 + brojac);
  RgQmatrixFF(ii + pMax, jj) = resultRgQ(disps(ranki) + 2*sizeVec2 + brojac);
  RgQmatrixFF(ii + pMax, jj + pMax) = resultRgQ(disps(ranki) + 3*sizeVec2 + brojac);
 
   brojac++;
 
   } // for jj
   } // for ii
 
   }  // for ranki
 
  
  TmatrixFF = - RgQmatrixFF * QmatrixFF.inverse();  
  TRgQmatrixFF.block(0, objIndex*2*pMax, 2*pMax, 2*pMax) = TmatrixFF;
  TRgQmatrixFF.block(0, nobj*2*pMax + objIndex*2*pMax, 2*pMax, 2*pMax) = RgQmatrixFF;

 
  MPI_Barrier(MPI_COMM_WORLD);
}// loop on many targets

               
 return TRgQmatrixFF;
}

Matrix<t_complex> getTRgQmatrix_SH_parr(Geometry const &geometry,
                                                  std::shared_ptr<Excitation const> incWave) {
auto const nobj = geometry.objects.size();
 int gran1, gran2, gran11, gran22, rank, size, pMax, sizeVec2;
 mpi::Communicator communicator;
 rank = communicator.rank();
 size = communicator.size();
 auto const nMaxS = geometry.objects.front().nMaxS;

 pMax = nMaxS * (nMaxS + 2);

    if (rank < (pMax % size)) {
    gran1 = rank * (pMax/size + 1);
    gran2 = gran1 + pMax/size + 1;
    } else {
    gran1 = rank * (pMax/size) + (pMax % size);
    gran2 = gran1 + (pMax/size);
    }


 int numRow = gran2 - gran1;

 int sizeVec = 4 * numRow * pMax;

 Vector<t_complex> QmatrixSH_proc (sizeVec);
 Vector<t_complex> RgQmatrixSH_proc (sizeVec);
 Vector<t_complex> resultQ (4 * pMax * pMax);
 Vector<t_complex> resultRgQ (4 * pMax * pMax);
 Matrix<t_complex> QmatrixSH(2*pMax, 2*pMax), RgQmatrixSH(2*pMax, 2*pMax), TmatrixSH(2*pMax, 2*pMax), TRgQmatrixSH(2*pMax, 4*nobj*pMax);
 TRgQmatrixSH.setZero();

for (int objIndex = 0; objIndex < nobj; objIndex++){

 QmatrixSH_proc = getQmatrix_SH(geometry, geometry.bground, incWave, gran1, gran2, objIndex);
 RgQmatrixSH_proc = getRgQmatrix_SH(geometry, geometry.bground, incWave, gran1, gran2, objIndex);

    Vector<int> sizesProc(size), disps(size);

   MPI_Allgather (&sizeVec, 1, MPI_INT, &sizesProc(0), 1, MPI_INT, MPI_COMM_WORLD);

   for (int kk = 0; kk < size; kk++)
   disps(kk) = (kk > 0) ? (disps(kk-1) + sizesProc(kk-1)) : 0; // displacements

  MPI_Gatherv (&QmatrixSH_proc(0), sizeVec, MPI_DOUBLE_COMPLEX, &resultQ(0), &sizesProc(0), &disps(0), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Gatherv (&RgQmatrixSH_proc(0), sizeVec, MPI_DOUBLE_COMPLEX, &resultRgQ(0), &sizesProc(0), &disps(0), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  
  MPI_Bcast(&resultQ(0), 4*pMax*pMax, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(&resultRgQ(0), 4*pMax*pMax, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

 // rearranging vectors to matrices 
 
    for (int ranki = 0; ranki < size; ranki++){

     if (ranki < (pMax % size)) {
    gran11 = ranki * (pMax/size + 1);
    gran22 = gran11 + pMax/size + 1;
    } else {
    gran11 = ranki * (pMax/size) + (pMax % size);
    gran22 = gran11 + (pMax/size);
    }

    sizeVec2 = pMax * (gran22 - gran11);

  int brojac (0);
 //single target only
 
for (int ii = gran11; ii < gran22; ii++){
    for (int jj = 0; jj < pMax; jj++){

  QmatrixSH(ii, jj) = resultQ(disps(ranki) + brojac);
  QmatrixSH(ii, jj + pMax) = resultQ(disps(ranki) + sizeVec2 + brojac);
  QmatrixSH(ii + pMax, jj) = resultQ(disps(ranki) + 2*sizeVec2 + brojac);
  QmatrixSH(ii + pMax, jj + pMax) = resultQ(disps(ranki) + 3*sizeVec2 + brojac);

  RgQmatrixSH(ii, jj) = resultRgQ(disps(ranki) + brojac);
  RgQmatrixSH(ii, jj + pMax) = resultRgQ(disps(ranki) + sizeVec2 + brojac);
  RgQmatrixSH(ii + pMax, jj) = resultRgQ(disps(ranki) + 2*sizeVec2 + brojac);
  RgQmatrixSH(ii + pMax, jj + pMax) = resultRgQ(disps(ranki) + 3*sizeVec2 + brojac);

   brojac++;

   } // for jj
   } // for ii

}  // for ranki

  TmatrixSH = RgQmatrixSH * QmatrixSH.inverse();
  TRgQmatrixSH.block(0, objIndex*2*pMax, 2*pMax, 2*pMax) = TmatrixSH;
  TRgQmatrixSH.block (0, nobj*2*pMax + objIndex*2*pMax, 2*pMax, 2*pMax) = RgQmatrixSH;

  MPI_Barrier(MPI_COMM_WORLD);

}// loop over many targets

 return TRgQmatrixSH;
 
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

#ifdef OPTIMET_MPI
Vector<t_complex> source_vectorSH_parallelAR3(Geometry &geometry, int gran1, int gran2,
                                std::shared_ptr<Excitation const> incWave, Vector<t_complex> &internalCoef_FF_, 
                                 Vector<t_complex> &scatteredCoef_FF_, int objIndex) {

if(gran1 == gran2)
  return Vector<t_complex>::Zero(0);

  Vector<t_complex> resultProc(2*(gran2 - gran1));
  
  geometry.getEXCvecSH_ARB3_parall(resultProc, incWave, scatteredCoef_FF_, internalCoef_FF_, gran1, gran2, objIndex);

  return resultProc;
}


Vector<t_complex> source_vectorSH_parallelAR1(Geometry &geometry, int gran1, int gran2,
                                std::shared_ptr<Excitation const> incWave, Vector<t_complex> &internalCoef_FF_,
                                 Vector<t_complex> &scatteredCoef_FF_, int objIndex) {
if(gran1 == gran2)
  return Vector<t_complex>::Zero(0);

  Vector<t_complex> resultProc(2*(gran2 - gran1));
  
if (geometry.objects[0].scatterer_type == "arbitrary.shape"){

geometry.getEXCvecSH_ARB1_parall(resultProc, incWave, scatteredCoef_FF_, internalCoef_FF_, gran1, gran2, objIndex);

}

return resultProc;

}
#endif

Vector<t_complex> source_vector(std::vector<Scatterer> const &objects, std::shared_ptr<Excitation const> incWave) {
  return source_vector(objects.begin(), objects.end(), incWave);
}



Vector<t_complex> source_vector(Geometry const &geometry, std::shared_ptr<Excitation const> incWave) {
  if(geometry.objects.size() == 0)
    return Vector<t_complex>(0, 0);
  // Check nMax is same accross all objects
  auto const nMax = geometry.objects.front().nMax;
  for(auto const &scatterer : geometry.objects)
    if(scatterer.nMax != nMax)
      throw std::runtime_error("All objects must have same number of harmonics");
  return source_vector(geometry.objects, incWave);
}

}
