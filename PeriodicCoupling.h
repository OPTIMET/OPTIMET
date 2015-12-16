/*
 * Periodic_Couplin.h
 *
 *  Created on: June 27, 2013
 *      Author: uceealj
 */

#ifndef Periodic_Coupling_H_
#define Periodic_Coupling_H_

#include "Spherical.h"
#include "Cartesian.h"
#include <complex>

class PeriodicCoupling
{
  private:

  public:
    static void compute_AlBe_nmlk(Spherical<double> R, std::complex<double> waveK, int BHreg, int n_max, std::complex<double> ****AlBe_nmlk);

    static std::complex<double> compute_Bnm(int n1, int n2, int n3, int m1, int m2, int m3);

    static int compute_Znmlk(double a1, double a2, std::complex<double> waveK, int Rmax, int n_max, std::complex<double> ****Z_nmlk);

    static int compute_OMEGAnmlk(std::complex<double> waveK, int n_max, std::complex<double> **dataOMEGA_1pq, std::complex<double> **dataOMEGA_2pq);

    static int compute_Ap_vkg(int n_max, Cartesian<double> vecKg, std::complex<double> waveK, double Ao, Cartesian<std::complex<double> > *dataAp_vkg);

    static double a_nm_p(double n, double m);
    static double a_nm_m(double n, double m);
    static double b_nm_p(double n, double m);
    static double b_nm_m(double n, double m);

    static double compute_alpha_nm(int n, int m);
    static double compute_beta_nm(int n, int m);
};

#endif /*Periodic_Coupling_H_ */
