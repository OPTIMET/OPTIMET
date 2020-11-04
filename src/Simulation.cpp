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

#include "Simulation.h"

#include "Aliases.h"
#include "CompoundIterator.h"
#include "Output.h"
#include "Reader.h"
#include "Result.h"
#include "Run.h"
#include "Solver.h"
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace std::chrono;

namespace optimet {
int Simulation::run() {
  
  // Read the case file
  auto run = simulation_input(caseFile + ".xml");
  
#ifdef OPTIMET_MPI
  run.parallel_params.grid = scalapack::squarest_largest_grid(communicator().size());
  run.communicator = communicator();
#endif

  // Initialize the solver
  auto const solver = optimet::solver::factory(run);

  switch(run.outputType) {
  case 0:
    field_simulation(run, solver);
    break;
  case 11:
    scan_wavelengths(run, solver);
    break;
  case 12:
    radius_scan(run, solver);
    break;
  case 112:
    radius_and_wavelength_scan(run, solver);
    break;
  case 2:
    coefficients(run, solver);
    break;
  default:
    std::cerr << "Nothing to do?\n";
    return 1;
  }
  return 0;
}

void Simulation::field_simulation(Run &run, std::shared_ptr<solver::AbstractSolver> solver) {
  // Determine the simulation type and proceed accordingly

  Result result(run.geometry, run.excitation);

   int nMax = run.geometry->nMax();
  int nMaxS = run.geometry->nMaxS();
  int flatMax = nMax * (nMax + 2);
  int flatMaxS = nMaxS * (nMaxS + 2);
  int sizeCF = flatMaxS * flatMax * flatMax;
  
  int gran1CG, gran2CG;
  int rank_pc = communicator().rank();
  int size_r = communicator().size();

      gran1CG = (sizeCF / (size_r))*(rank_pc);

      if(rank_pc < (size_r-1)) {

      gran2CG = (rank_pc + 1)*(sizeCF / (size_r)); }

        else { gran2CG = sizeCF;}

   int sizeCF_par = gran2CG - gran1CG;



  std::vector<double> C_10m1_par(sizeCF_par), C_11m1_par(sizeCF_par), C_00m1_par(sizeCF_par), C_01m1_par(sizeCF_par);
  std::vector<double> W_m1m1_par(sizeCF_par), W_11_par(sizeCF_par), W_00_par(sizeCF_par), W_10_par(sizeCF_par), W_01_par(sizeCF_par);

  std::vector<double> C_10m1(sizeCF), C_11m1(sizeCF), C_00m1(sizeCF), C_01m1(sizeCF);
  std::vector<double> W_m1m1(sizeCF), W_11(sizeCF), W_00(sizeCF), W_10(sizeCF), W_01(sizeCF);

  std::vector<double *> CLGcoeff_par = {&C_10m1_par[0], &C_11m1_par[0], &C_00m1_par[0], &C_01m1_par[0],
                                       &W_m1m1_par[0], &W_11_par[0], &W_00_par[0], &W_10_par[0], &W_01_par[0]};

  std::vector<double *> CLGcoeff = {&C_10m1[0], &C_11m1[0], &C_00m1[0], &C_01m1[0],
                                       &W_m1m1[0], &W_11[0], &W_00[0], &W_10[0], &W_01[0]};

  run.geometry->Coefficients(nMax, nMaxS, CLGcoeff_par, gran1CG, gran2CG);

  All2all(CLGcoeff, CLGcoeff_par, sizeCF_par);

  solver->solve(result.scatter_coef, result.internal_coef, result.scatter_coef_SH, result.internal_coef_SH, CLGcoeff);

  std::vector<double> Rr, Rthe, Rphi;

  int size, sizeField;

  if(communicator().rank() == communicator().root_id()) {

     OutputGrid oEGrid_FF(O3DCartesianRegular, run.params);
    OutputGrid oHGrid_FF(O3DCartesianRegular, run.params);

    OutputGrid oEGrid_SH(O3DCartesianRegular, run.params);
    OutputGrid oHGrid_SH(O3DCartesianRegular, run.params);
   
  Spherical<double> Rloc;
  
  int gran = oEGrid_FF.gridPoints / (communicator().size()-1);

  int brojac(1), numproc(1);  
   
  while(brojac<=oEGrid_FF.gridPoints) {
    Rloc = oEGrid_FF.getPoint();
    oHGrid_FF.getPoint();
    oEGrid_SH.getPoint();
    oHGrid_SH.getPoint();
   
    Rr.push_back(Rloc.rrr);
    Rthe.push_back(Rloc.the);
    Rphi.push_back(Rloc.phi);
    
    if ((brojac == numproc * gran)&&(numproc < (communicator().size()-1))){
       size = Rr.size();    
       MPI_Send(&Rr[0], Rr.size(), MPI_DOUBLE, numproc, 1, MPI_COMM_WORLD);
       MPI_Send(&Rthe[0], Rthe.size(), MPI_DOUBLE, numproc, 2, MPI_COMM_WORLD);
       MPI_Send(&Rphi[0], Rphi.size(), MPI_DOUBLE, numproc, 3, MPI_COMM_WORLD);           
       MPI_Send(&size, 1, MPI_INT, numproc, 4, MPI_COMM_WORLD);
     numproc++;

     Rr.clear();
     Rthe.clear();
     Rphi.clear();

     }

    else if(brojac == oEGrid_FF.gridPoints){
       size = Rr.size();
       MPI_Send(&Rr[0], Rr.size(), MPI_DOUBLE, numproc, 1, MPI_COMM_WORLD);
       MPI_Send(&Rthe[0], Rthe.size(), MPI_DOUBLE, numproc, 2, MPI_COMM_WORLD);
       MPI_Send(&Rphi[0], Rphi.size(), MPI_DOUBLE, numproc, 3, MPI_COMM_WORLD);
       MPI_Send(&size, 1, MPI_INT, numproc, 4, MPI_COMM_WORLD);
    }

 brojac++;

    oEGrid_FF.gotoNext();
    oHGrid_FF.gotoNext();
    oEGrid_SH.gotoNext();
    oHGrid_SH.gotoNext();
} //while     
 
   } // root_id 0 proccess

if(communicator().rank() != communicator().root_id()) {
     MPI_Recv(&size, 1, MPI_INT, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     Rr.resize(size);
     Rthe.resize(size);
     Rphi.resize(size);
   MPI_Recv(&Rr[0], Rr.size(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   MPI_Recv(&Rthe[0], Rthe.size(), MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   MPI_Recv(&Rphi[0], Rphi.size(), MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  
result.setFields(Rr, Rthe, Rphi, run.projection, CLGcoeff);
}

if(communicator().rank() == communicator().root_id()) {

    Output oFile_FF(caseFile + "_FF.h5");
    Output oFile_SH(caseFile + "_SH.h5");

    OutputGrid oEGrid_FF2(O3DCartesianRegular, run.params, oFile_FF.getHandle("Field_E"));
    OutputGrid oHGrid_FF2(O3DCartesianRegular, run.params, oFile_FF.getHandle("Field_H"));

    OutputGrid oEGrid_SH2(O3DCartesianRegular, run.params, oFile_SH.getHandle("Field_E"));
    OutputGrid oHGrid_SH2(O3DCartesianRegular, run.params, oFile_SH.getHandle("Field_H"));

    int NOpoints =  oEGrid_FF2.gridPoints;

    std::vector<SphericalP<std::complex<double>>> EField_FFvec (NOpoints);
    std::vector<SphericalP<std::complex<double>>> HField_FFvec (NOpoints);
    std::vector<SphericalP<std::complex<double>>> EField_SHvec (NOpoints);
    std::vector<SphericalP<std::complex<double>>> HField_SHvec (NOpoints);


    SphericalP<std::complex<double>> EField_FF;
    SphericalP<std::complex<double>> HField_FF;
    SphericalP<std::complex<double>> EField_SH;
    SphericalP<std::complex<double>> HField_SH;


    std::vector<std::complex<double>> EField_FF_x;
    std::vector<std::complex<double>> EField_FF_y;
    std::vector<std::complex<double>> EField_FF_z;
    std::vector<std::complex<double>> HField_FF_x;
    std::vector<std::complex<double>> HField_FF_y;
    std::vector<std::complex<double>> HField_FF_z;

    std::vector<std::complex<double>> EField_SH_x;
    std::vector<std::complex<double>> EField_SH_y;
    std::vector<std::complex<double>> EField_SH_z;
    std::vector<std::complex<double>> HField_SH_x;
    std::vector<std::complex<double>> HField_SH_y;
    std::vector<std::complex<double>> HField_SH_z;

    int kk = 0;

   for(int rank = 1; rank < communicator().size(); rank++){
    
   MPI_Recv(&sizeField, 1, MPI_INT, rank, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

   EField_FF_x.resize(sizeField);
   EField_FF_y.resize(sizeField);
   EField_FF_z.resize(sizeField);
   HField_FF_x.resize(sizeField);
   HField_FF_y.resize(sizeField);
   HField_FF_z.resize(sizeField);

   EField_SH_x.resize(sizeField);
   EField_SH_y.resize(sizeField);
   EField_SH_z.resize(sizeField);
   HField_SH_x.resize(sizeField);
   HField_SH_y.resize(sizeField);
   HField_SH_z.resize(sizeField);

   MPI_Recv(&EField_FF_x[0], EField_FF_x.size(), MPI_DOUBLE_COMPLEX, rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   MPI_Recv(&EField_FF_y[0], EField_FF_y.size(), MPI_DOUBLE_COMPLEX, rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   MPI_Recv(&EField_FF_z[0], EField_FF_z.size(), MPI_DOUBLE_COMPLEX, rank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   MPI_Recv(&HField_FF_x[0], HField_FF_x.size(), MPI_DOUBLE_COMPLEX, rank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   MPI_Recv(&HField_FF_y[0], HField_FF_y.size(), MPI_DOUBLE_COMPLEX, rank, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   MPI_Recv(&HField_FF_z[0], HField_FF_z.size(), MPI_DOUBLE_COMPLEX, rank, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

   MPI_Recv(&EField_SH_x[0], EField_SH_x.size(), MPI_DOUBLE_COMPLEX, rank, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   MPI_Recv(&EField_SH_y[0], EField_SH_y.size(), MPI_DOUBLE_COMPLEX, rank, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   MPI_Recv(&EField_SH_z[0], EField_SH_z.size(), MPI_DOUBLE_COMPLEX, rank, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   MPI_Recv(&HField_SH_x[0], HField_SH_x.size(), MPI_DOUBLE_COMPLEX, rank, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   MPI_Recv(&HField_SH_y[0], HField_SH_y.size(), MPI_DOUBLE_COMPLEX, rank, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   MPI_Recv(&HField_SH_z[0], HField_SH_z.size(), MPI_DOUBLE_COMPLEX, rank, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
   for (int ii = 0; ii < EField_FF_x.size(); ii++){

   EField_FF.rrr = EField_FF_x[ii];   
   EField_FF.the = EField_FF_y[ii];
   EField_FF.phi = EField_FF_z[ii];
   HField_FF.rrr = HField_FF_x[ii];
   HField_FF.the = HField_FF_y[ii];
   HField_FF.phi = HField_FF_z[ii];

   EField_SH.rrr = EField_SH_x[ii];
   EField_SH.the = EField_SH_y[ii];
   EField_SH.phi = EField_SH_z[ii];
   HField_SH.rrr = HField_SH_x[ii];
   HField_SH.the = HField_SH_y[ii];
   HField_SH.phi = HField_SH_z[ii];

    EField_FFvec[kk] = EField_FF;
    HField_FFvec[kk] = HField_FF;
    EField_SHvec[kk] = EField_SH;
    HField_SHvec[kk] = HField_SH;
 
   kk++; 
       
   }

   }

    kk = 0;

    while(!((oEGrid_FF2.gridDone)&&(oEGrid_SH2.gridDone))) {
    oEGrid_FF2.getPoint();
    oHGrid_FF2.getPoint();
    oEGrid_SH2.getPoint();
    oHGrid_SH2.getPoint(); 

    oHGrid_FF2.pushDataNext(HField_FFvec[kk]);
    oEGrid_FF2.pushDataNext(EField_FFvec[kk]);
    oHGrid_SH2.pushDataNext(HField_SHvec[kk]);
    oEGrid_SH2.pushDataNext(EField_SHvec[kk]);
   
    kk++;
  }

    oEGrid_FF2.close();
    oHGrid_FF2.close();
    oFile_FF.close();
    oEGrid_SH2.close();
    oHGrid_SH2.close();
    oFile_SH.close();
  }
}

void Simulation::All2all (std::vector<double *> CLGcoeff, std::vector<double *> CLGcoeff_par, int sizeVec){

    int size = communicator().size();
    Vector<int> sizesProc(size), dispr(size);

    MPI_Allgather(&sizeVec, 1, MPI_INT, &sizesProc(0), 1, MPI_INT, MPI_COMM_WORLD);
 
  
   for (int kk = 0; kk < size; kk++)
   dispr(kk) = (kk > 0) ? (dispr(kk-1) + sizesProc(kk-1)) : 0;
 
  
    for (int i = 0; i < 9; i++)
    MPI_Allgatherv(CLGcoeff_par[i], sizeVec, MPI_DOUBLE, CLGcoeff[i],
                                       &sizesProc(0), &dispr(0), MPI_DOUBLE, MPI_COMM_WORLD);


}


void Simulation::scan_wavelengths(Run &run, std::shared_ptr<solver::AbstractSolver> solver) {
  std::ofstream outASec_FF, outESec_FF, outSSec_SH, outASec_SH;
 
  int nMax = run.geometry->nMax();
  int nMaxS = run.geometry->nMaxS();
  int flatMax = nMax * (nMax + 2);
  int flatMaxS = nMaxS * (nMaxS + 2);
  int sizeCF = flatMaxS * flatMax * flatMax;
  int gran1, gran2, gran1AC, gran2AC, gran1CG, gran2CG;
  int rank = communicator().rank();
  int size = communicator().size(); 
 
      gran1CG = (sizeCF / (size))*(rank);

      if(rank<(size-1)) {

      gran2CG = (rank+1)*(sizeCF / (size)); }

       	else { gran2CG = sizeCF;}

   int sizeCF_par = gran2CG - gran1CG;
 
  
  
  std::vector<double> C_10m1_par(sizeCF_par), C_11m1_par(sizeCF_par), C_00m1_par(sizeCF_par), C_01m1_par(sizeCF_par);
  std::vector<double> W_m1m1_par(sizeCF_par), W_11_par(sizeCF_par), W_00_par(sizeCF_par), W_10_par(sizeCF_par), W_01_par(sizeCF_par);

  std::vector<double> C_10m1(sizeCF), C_11m1(sizeCF), C_00m1(sizeCF), C_01m1(sizeCF);
  std::vector<double> W_m1m1(sizeCF), W_11(sizeCF), W_00(sizeCF), W_10(sizeCF), W_01(sizeCF);
  
  std::vector<double *> CLGcoeff_par = {&C_10m1_par[0], &C_11m1_par[0], &C_00m1_par[0], &C_01m1_par[0],
                                       &W_m1m1_par[0], &W_11_par[0], &W_00_par[0], &W_10_par[0], &W_01_par[0]};

  std::vector<double *> CLGcoeff = {&C_10m1[0], &C_11m1[0], &C_00m1[0], &C_01m1[0],
                                       &W_m1m1[0], &W_11[0], &W_00[0], &W_10[0], &W_01[0]};

  auto start1 = high_resolution_clock::now(); 

  run.geometry->Coefficients(nMax, nMaxS, CLGcoeff_par, gran1CG, gran2CG);

      
  All2all(CLGcoeff, CLGcoeff_par, sizeCF_par);
 
    auto stop1 = high_resolution_clock::now();
    auto duration1 = duration_cast<microseconds>(stop1 - start1);

   if(communicator().rank() == 0) {
    //  std::cout << "CG coeff asse-" << std::endl;
    //  std::cout << duration1.count()/1e6 <<"e-0"<< std::endl;
    }

  
  if(communicator().rank() == communicator().root_id()) {
  
 
    outASec_FF.open(caseFile + "_AbsorptionCS_FF.dat");
    outESec_FF.open(caseFile + "_ExtinctionCS_FF.dat");
    
   
    outSSec_SH.open(caseFile + "_ScatteringCS_SH.dat");
    outASec_SH.open(caseFile + "_AbsorptionCS_SH.dat");
    
  }

  // Now scan over the wavelengths given in params
  double lami = run.params[0];
  double lamf = run.params[1];
  int steps = run.params[2];

  double lam;
  double lams;

  double absCS_SH(0.0), scaCS_SH(0.0), absCS_FF(0.0), extCS_FF(0.0);
  int NO = run.geometry->objects.size();
  int sizet = NO;

  Vector<double> absCS_SH_vec(size), scaCS_SH_vec(size), absCS_FF_vec(size), extCS_FF_vec(size);
  
  lams = (lamf - lami) / (steps - 1); 
  
   for(int i = 0; i < steps; i++) {
    
    lam = lami + i * lams;

    run.excitation->updateWavelength(lam);
    run.geometry->update(run.excitation);

    auto start2 = high_resolution_clock::now();
    solver->update(run);
   
    int nMaxS = run.nMaxS;
    t_uint const pMax = nMaxS * (nMaxS + 2);
    int TMax = NO * pMax;
    

    Result result(run.geometry, run.excitation);

    if(communicator().rank() == communicator().root_id()) {

      std::cout << "Solving for Lambda = " << lam << std::endl;

    }
   
  solver->solve(result.scatter_coef, result.internal_coef, result.scatter_coef_SH, result.internal_coef_SH, CLGcoeff);


    auto stop2 = high_resolution_clock::now();
    auto duration2 = duration_cast<microseconds>(stop2 - start2);

    if(communicator().rank() == 0) {
   // std::cout << "FF and SH assembly and solve-" << std::endl;
   // std::cout << duration2.count()/1e6<<"e-0"<< std::endl; 
    }


      if (size <= NO) {  // if the number of processes is less or eq numb of part

      gran1 = (NO / (size))*(rank);

      if(rank<(size-1)) {

      gran2 = (rank+1)*(NO / (size));}

	else { gran2 = NO;}

    gran1AC = (TMax / (size))*(rank);

      if(rank<(size-1)) {

      gran2AC = (rank+1)*(TMax / (size)); }

       	else { gran2AC = TMax;}

   
     absCS_SH = result.getAbsorptionCrossSection_SH(CLGcoeff, gran1AC, gran2AC);
     scaCS_SH = result.getScatteringCrossSection_SH(gran1, gran2);
     absCS_FF = result.getAbsorptionCrossSection(gran1, gran2);
     extCS_FF = result.getExtinctionCrossSection(gran1, gran2);

    MPI_Gather(&absCS_SH, 1, MPI_DOUBLE, &absCS_SH_vec(0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);    
    MPI_Gather(&scaCS_SH, 1, MPI_DOUBLE, &scaCS_SH_vec(0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);     
    MPI_Gather(&absCS_FF, 1, MPI_DOUBLE, &absCS_FF_vec(0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&extCS_FF, 1, MPI_DOUBLE, &extCS_FF_vec(0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

 
      if(communicator().rank() == communicator().root_id()) {    
      
     //  std::cout<<absCS_FF_vec.sum()<<std::endl;    
     //    std::cout<<absCS_SH_vec.sum()<<std::endl;
              
      outASec_FF << lam << "\t" << absCS_FF_vec.sum() << std::endl;
      outESec_FF << lam << "\t" << extCS_FF_vec.sum() << std::endl;

      outSSec_SH << lam << "\t" << scaCS_SH_vec.sum() << std::endl;
      outASec_SH << lam << "\t" << absCS_SH_vec.sum() << std::endl;
  }
    
  }// if


 else if ((size > NO))  {  // if the number of processes is more than numb of part

    gran1AC = (TMax / (size))*(rank);

      if(rank<(size-1)) {

      gran2AC = (rank+1)*(TMax / (size));}

       	else { gran2AC = TMax;}
    
   //  absCS_SH = result.getAbsorptionCrossSection_SH(CLGcoeff, gran1AC, gran2AC);
    
   //  MPI_Gather(&absCS_SH, 1, MPI_DOUBLE, &absCS_SH_vec(0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

     if(communicator().rank() == communicator().root_id()) {

       // std::cout<<absCS_SH_vec.sum()<<std::endl;
  
      outASec_SH << lam << "\t" << absCS_SH_vec.sum() << std::endl;

      }

     if (rank<NO){     
     
     gran1 = rank;

     gran2 = rank + 1;
 
     scaCS_SH = result.getScatteringCrossSection_SH(gran1, gran2);
     absCS_FF = result.getAbsorptionCrossSection(gran1, gran2);
     extCS_FF = result.getExtinctionCrossSection(gran1, gran2);

   }//if
 
    MPI_Gather(&scaCS_SH, 1, MPI_DOUBLE, &scaCS_SH_vec(0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&absCS_FF, 1, MPI_DOUBLE, &absCS_FF_vec(0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&extCS_FF, 1, MPI_DOUBLE, &extCS_FF_vec(0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


      if(communicator().rank() == communicator().root_id()) {

      std::cout<<scaCS_SH_vec.sum()<<std::endl;
      outASec_FF << lam << "\t" << absCS_FF_vec.sum() << std::endl;
      outESec_FF << lam << "\t" << extCS_FF_vec.sum() << std::endl;
      outSSec_SH << lam << "\t" << scaCS_SH_vec.sum() << std::endl;

      }

      
}// else if

}// for


  if(communicator().rank() == communicator().root_id()) {

    outASec_FF.close();
    outESec_FF.close();
    
    outSSec_SH.close();
    outASec_SH.close();
  }
}

void Simulation::radius_scan(Run &run, std::shared_ptr<solver::AbstractSolver> solver) {
  std::ofstream outASec, outESec;

  if(communicator().rank() == communicator().root_id()) {
    outASec.open(caseFile + "_AbsorptionCS.dat");
    outESec.open(caseFile + "_ExtinctionCS.dat");
  }

  // Now scan over the wavelengths given in params
  double radi = run.params[3];
  double radf = run.params[4];
  int radsteps = run.params[5];

  double rad;
  double rads;

  int gran1, gran2;

  rads = (radf - radi) / (radsteps - 1);

  for(int i = 0; i < radsteps; i++) {
    rad = radi + i * rads;

    std::cout << "Solving for R = " << rad << std::endl;

    for(size_t k = 0; k < run.geometry->objects.size(); k++) {
      run.geometry->updateRadius(rad, k);
    }

    if(run.geometry->structureType == 1) {
      run.geometry->rebuildStructure();
    }

    if(!run.geometry->is_valid()) {
      std::cerr << "Geometry no longer valid!";
      exit(1);
    }

    solver->update(run);

    Result result(run.geometry, run.excitation);
    solver->solve(result.scatter_coef, result.internal_coef, result.scatter_coef_SH, result.internal_coef_SH, result.CLGcoeff);

    if(communicator().rank() == communicator().root_id()) {
      outASec << rad << "\t" << result.getAbsorptionCrossSection(gran1, gran2) << std::endl;
      outESec << rad << "\t" << result.getExtinctionCrossSection(gran1, gran2) << std::endl;
    }
  }

  if(communicator().rank() == communicator().root_id()) {
    outASec.close();
    outESec.close();
  }
}

void Simulation::radius_and_wavelength_scan(Run &run,
                                            std::shared_ptr<solver::AbstractSolver> solver) {
  std::ofstream outASec, outESec, outParams;

  if(communicator().rank() == communicator().root_id()) {
    outASec.open(caseFile + "_AbsorptionCS.dat");
    outESec.open(caseFile + "_ExtinctionCS.dat");
    outParams.open(caseFile + "_RadiusLambda.dat");
  }

  // Now scan over the wavelengths given in params
  double lami = run.params[0];
  double lamf = run.params[1];
  int lamsteps = run.params[2];

  double radi = run.params[3];
  double radf = run.params[4];
  int radsteps = run.params[5];

  double lam, lams, rad, rads;

  lams = (lamf - lami) / (lamsteps - 1);
  rads = (radf - radi) / (radsteps - 1);

  for(int i = 0; i < lamsteps; i++) {
    lam = lami + i * lams;

    for(int j = 0; j < radsteps; j++) {
      rad = radi + j * rads;

      std::cout << "Solving for Lambda = " << lam << " and R =" << rad << std::endl;

      run.excitation->updateWavelength(lam);
      run.geometry->update(run.excitation);
      for(size_t k = 0; k < run.geometry->objects.size(); k++) {
        run.geometry->updateRadius(rad, k);
      }

      if(run.geometry->structureType == 1) {
        run.geometry->rebuildStructure();
      }

      if(!run.geometry->is_valid()) {
        std::cerr << "Geometry no longer valid!";
        exit(1);
      }

      solver->update(run);

      int gran1, gran2;

      Result result(run.geometry, run.excitation);
      solver->solve(result.scatter_coef, result.internal_coef, result.scatter_coef_SH, result.internal_coef_SH, result.CLGcoeff);

      if(communicator().rank() == communicator().root_id()) {
        outASec << result.getAbsorptionCrossSection(gran1, gran2) << "\t";
        outESec << result.getExtinctionCrossSection(gran1, gran2) << "\t";
        outParams << "(" << rad * 1e9 << " , " << lam * 1e9 << ")"
                  << "\t";
      }
    }

    if(communicator().rank() == communicator().root_id()) {
      outASec << std::endl;
      outESec << std::endl;
      outParams << std::endl;
    }
  }

  if(communicator().rank() == communicator().root_id()) {
    outASec.close();
    outESec.close();
    outParams.close();
  }
}

void Simulation::coefficients(Run &run, std::shared_ptr<solver::AbstractSolver> solver) {
  // Scattering coefficients requests

  Result result(run.geometry, run.excitation);
  solver->solve(result.scatter_coef, result.internal_coef, result.scatter_coef_SH, result.internal_coef_SH, result.CLGcoeff);

  if(communicator().rank() == communicator().root_id()) {
    std::ofstream outPCoef(caseFile + "_pCoefficients.dat");
    std::ofstream outQCoef(caseFile + "_qCoefficients.dat");

    for(CompoundIterator p = 0; p < p.max(run.nMax); p++) {
      outPCoef << p.first << "\t" << p.second << "\t"
               << abs(result.scatter_coef(static_cast<int>(p))) << std::endl;
      outQCoef << p.first << "\t" << p.second << "\t"
               << abs(result.scatter_coef(static_cast<int>(p.compound) + p.max(run.nMax)))
               << std::endl;
    }

    outPCoef.close();
    outQCoef.close();
  }
}

int Simulation::done() {
  // Placeholder method. Not needed at the moment.
  return 0;
}
}
