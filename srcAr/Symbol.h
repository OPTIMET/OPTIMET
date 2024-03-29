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

/**
 * The Symbol class contains special symbol routines.
 * Symbol is released under the GSL. To view a copy
 * of the licence, look in the documentation.
 */

#ifndef OPTIMET_SYMBOL_H
#define OPTIMET_SYMBOL_H

#include "Scatterer.h"
#include "ElectroMagnetic.h"
#include "CompoundIterator.h"

namespace optimet {
namespace symbol {
 
double Wigner3j(int j1, int j2, int j3, int m1, int m2, int m3);
                           
std::complex<double> CXm1(CompoundIterator &kk, double *W_m1m1, double *W_00, double *W_11, int nMax,
                            optimet::Vector<optimet::t_complex> &internalCoef_FF_,
                            double r,
                            int objectIndex_, double omega,
                            const Scatterer &object);
                            
std::complex<double> CXp1(CompoundIterator &kk, double *W_m1m1, double *W_00, double *W_11, int nMax,
                            optimet::Vector<optimet::t_complex> &internalCoef_FF_,
                            double r,
                            int objectIndex_, double omega,
                            const Scatterer &object);                            
#ifdef OPTIMET_MPI                                                       
void C_10m1coeff (double *C_10m1, int nMax, int nMaxS, int gran1, int gran2); 
void C_11m1coeff (double *C_11m1, int nMax, int nMaxS, int gran1, int gran2);  
void C_00m1coeff (double *C_00m1, int nMax, int nMaxS, int gran1, int gran2);  
void C_01m1coeff (double *C_01m1, int nMax, int nMaxS, int gran1, int gran2);  
void W_m1m1coeff  (double *W_m1m1, int nMax, int nMaxS, int gran1, int gran2);    
void W_11coeff  (double *W_11, int nMax, int nMaxS, int gran1, int gran2); 
void W_00coeff  (double *W_00, int nMax, int nMaxS, int gran1, int gran2);
void W_10coeff  (double *W_10, int nMax, int nMaxS, int gran1, int gran2);  
void W_01coeff  (double *W_01, int nMax, int nMaxS, int gran1, int gran2);
#endif  
} // namespace symbol
} // namespace optimet

#endif /* OPTIMET_SYMBOL_H */
