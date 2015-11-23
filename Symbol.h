/**
 * The Symbol class contains special symbol routines.
 * Symbol is released under the GSL. To view a copy
 * of the licence, look in the documentation.
 */

#ifndef SPECFUN_H_
#define SPECFUN_H_

#include "Scatterer.h"
#include "ElectroMagnetic.h"

class Symbol
{
private:
	//
public:
	Symbol();
	virtual ~Symbol();

	/**
	 * Calculates the Wigner 3j symbol: \n
	 * (j1 j2 j3 \n
	 *  m1 m2 m3).
	 * @param j1 coefficient.
	 * @param j2 coefficient.
	 * @param j3 coefficient.
	 * @param m1 coefficient.
	 * @param m2 coefficient.
	 * @param m3 coefficient.
	 * @return the Wigner 3j symbol.
	 */
	static double Wigner3j(int j1, int j2, int j3, int m1, int m2, int m3);

	/**
	 * Calculates the Wigner 6j symbol: \n
	 * {j1 j2 j3 \n
	 *  j4 j5 j6}
	 * @param j1 coefficient.
	 * @param j2 coefficient.
	 * @param j3 coefficient.
	 * @param j4 coefficient.
	 * @param j5 coefficient.
	 * @param j6 coefficient.
	 * @return the Wigner 6j symbol.
	 */
	static double Wigner6j(int j1, int j2, int j3, int j4, int j5, int j6);

	/**
	 * Calculates the Wigner 9j symbol: \n
	 * {j11 j12 j13 \n
	 *  j21 j22 j23 \n
	 *  j31 j32 j33}
	 * @param j11 coefficient.
	 * @param j12 coefficient.
	 * @param j13 coefficient.
	 * @param j21 coefficient.
	 * @param j22 coefficient.
	 * @param j23 coefficient.
	 * @param j31 coefficient.
	 * @param j32 coefficient.
	 * @param j33 coefficient.
	 * @return the Wigner 9j symbol.
	 */
	static double Wigner9j(int j11, int j12, int j13, int j21, int j22, int j23, int j31, int j32, int j33);

	/**
	 * Calculates the Clebsch-Gordan coefficient:
	 * \f$C_{j_1 m_1 j_2 m_2}^{j m}\f$
	 * (i.e. the decomposition of |j,m> in terms of |j1, m1> |j2, m2>)
	 * @param j coefficient.
	 * @param m coefficient.
	 * @param j1 coefficient.
	 * @param m1 coefficient.
	 * @param j2 coefficient.
	 * @param m2 coefficient.
	 * @return the \f$C_{j_1 m_1 j_2 m_2}^{j m}\f$ coefficient.
	 */
	static double CleGor(int j, int m, int j1, int m1, int j2, int m2);


	// AJ ------------------------------------------------------------
	/**
	 * Calculate C^(0, 1, -1) number
	 * This is obtained from CleGor & Wigner9j
	 * @param J1 coefficient.
	 * @param M1 coefficient.
	 * @param J2 coefficient.
	 * @param M2 coefficient.
	 * @param J coefficient.
	 * @param M coefficient.
	 * @return the C^(0,1,-1)_(J1,M1,J2,M2,J,M)
	 */
	static double C_01m1(int J1, int M1, int J2, int M2, int J, int M);

	/**
	 * Calculate C^(1, 0, -1) number
	 * This is obtained from CleGor & Wigner9j
	 * @param J1 coefficient.
	 * @param M1 coefficient.
	 * @param J2 coefficient.
	 * @param M2 coefficient.
	 * @param J coefficient.
	 * @param M coefficient.
	 * @return the C^(1,0,-1)_(J1,M1,J2,M2,J,M)
	 */
	static double C_10m1(int J1, int M1, int J2, int M2, int J, int M);

	/**
	 * Calculate C^(0, 0, -1) number
	 * This is obtained from CleGor & Wigner9j
	 * @param J1 coefficient.
	 * @param M1 coefficient.
	 * @param J2 coefficient.
	 * @param M2 coefficient.
	 * @param J coefficient.
	 * @param M coefficient.
	 * @return the C^(0,0,-1)_(J1,M1,J2,M2,J,M)
	 */
	static double C_00m1(int J1, int M1, int J2, int M2, int J, int M);

	/**
	 * Calculate C^(1, 1, -1) number
	 * This is obtained from CleGor & Wigner9j
	 * @param J1 coefficient.
	 * @param M1 coefficient.
	 * @param J2 coefficient.
	 * @param M2 coefficient.
	 * @param J coefficient.
	 * @param M coefficient.
	 * @return the C^(1,1,-1)_(J1,M1,J2,M2,J,M)
	 */
	static double C_11m1(int J1, int M1, int J2, int M2, int J, int M);

	/**
	 * Calculate A^(0)_mn number
	 * @param m coefficient.
	 * @param n coefficient.
	 * @param r sphere radius.
	 * @param waveK_i wave number inside sphere.
	 * @param cmn_1 internal expansion coefficient at the fundamental frequency.
	 * @return the A^(0)_(m,n)
	 */
	static complex<double> A_0(int m, int n, double R, complex<double> waveK_i, complex<double> cmn_1);

	/**
	 * Calculate A^(1)_mn number
	 * @param m coefficient.
	 * @param n coefficient.
	 * @param r sphere radius.
	 * @param waveK_i wave number inside sphere.
	 * @param cmn_1 internal expansion coefficient at the fundamental frequency.
	 * @return the A^(0)_(m,n)
	 */
	static complex<double> A_1(int m, int n, double R, complex<double> waveK_i, complex<double> dmn_1);

	/**
	 * Calculate A^(m1)_mn number
	 * @param m coefficient.
	 * @param n coefficient.
	 * @param r sphere radius.
	 * @param waveK_i wave number inside sphere.
	 * @param cmn_1 internal expansion coefficient at the fundamental frequency.
	 * @return the A^(0)_(m,n)
	 */
	static complex<double> A_m1(int m, int n, double R, complex<double> waveK_i, complex<double> dmn_1);

	/**
	 * Calculate W^(L1, J1, M1, L2, J2, M2)_LM number
	 * @param L1 coefficient.
	 * @param J1 coefficient.
	 * @param M1 coefficient.
	 * @param L2 coefficient.
	 * @param J2 coefficient.
	 * @param M2 coefficient.
	 * @param L coefficient.
	 * @param M coefficient.
	 * @return the W^(L1, J1, M1, L2, J2, M2)_LM number
	 */
	static double W(int L1, int J1, int M1, int L2, int J2, int M2, int L, int M);

	static complex<double> up_mn(int m, int n, int nMax_, complex<double> cmn_1,complex<double> dmn_1,  double omega_, Scatterer *object_, ElectroMagnetic bground_);

	static complex<double> vp_mn(int m, int n, int nMax_, complex<double> cmn_1, complex<double> dmn_1, double omega_, Scatterer *object_, ElectroMagnetic bground_);

	static complex<double> upp_mn(int m, int n, int nMax_, complex<double> cmn_1, complex<double> dmn_1, double omega_, Scatterer *object_, ElectroMagnetic bground_);

	// ---------------------------------------------------------------
};

#endif /* SPECFUN_H_ */
