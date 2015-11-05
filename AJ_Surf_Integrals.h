/*
 * AJ_Surf_Integrals.h
 *
 *  Created on: Mar 4, 2013
 *      Author: uceealj
 */

#ifndef AJ_SURF_INTEGRALS_H_
#define AJ_SURF_INTEGRALS_H_

#include "AuxCoefficients.h"

// Decalre a four members structure : complex<double> --------------------------------
class Qclass{
    public:
    complex<double> Q11, Q12, Q21, Q22;
};
//------------------------------------------------------------------------------------

Qclass calcQ_rthephi(int n, int m, int np, int mp, int N, int BHreg, complex<double> wavek1, complex<double> wavek2, double a, double b, double c);

Qclass calcQ_rthe(int n, int m, int np, int mp, int N, int BHreg, complex<double> wavek1, complex<double> wavek2, double a, double b);

Qclass calcTsphere(int n, int m, int np, int mp, complex<double> wavek1, complex<double> wavek2, double a);

#endif /* AJ_SURF_INTEGRALS_H_ */
