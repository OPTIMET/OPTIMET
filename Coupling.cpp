#include "Coupling.h"

#include "Tools.h"
#include "CompoundIterator.h"
#include "PeriodicCoupling.h"
#include "constants.h"

#include <cmath>
#include <iostream>

double Coupling::a_nm_p(double n, double m)
{
  return -sqrt( ((n+m+1.)*(n-m+1.))/((2.*n+1.)*(2.*n+3.)) );
}

double Coupling::a_nm_m(double n, double m)
{
  return +sqrt( ((n+m)*(n-m))/((2.*n+1.)*(2.*n-1.)) );
}

double Coupling::b_nm_p(double n, double m)
{
  return +sqrt( ((n+m+2.)*(n+m+1.))/((2.*n+1.)*(2.*n+3.)) );
}

double Coupling::b_nm_m(double n, double m)
{
  return +sqrt( ((n-m)*(n-m-1.))/((2.*n+1.)*(2.*n-1.)) );
}

void Coupling::TransferCoefficients(Spherical<double> R, std::complex<double> waveK, int BHreg,
      int n_max, std::complex<double> **dataApq, std::complex<double> **dataBpq)
{
  // Translation coefficients calculation  -----------------------------------------
  // -------------------------------------------------------------------------------
  // temp variables
  std::complex<double> c_temp(0., 0.);
  // working variables
  int i(0), j(0);
  int ii(0), jj(0);
  int n(0), m(0), l(0), k(0);         // duale subscript variables
  double d_n(0.), d_m(0.), d_l(0.), d_k(0.);  // duale subscript variables
  // auxiliary variables
  std::complex<double> AlBe2(0., 0.);        // auxiliary coefficients
  std::complex<double> AlBe3(0., 0.);        // auxiliary coefficients
  std::complex<double> AlBe4(0., 0.);        // auxiliary coefficients
  double A1(0.), A2(0.), A3(0.);        // auxiliary coefficients
  double B1(0.), B2(0.), B3(0.);        // auxiliary coefficients
  int p=(0), q(0);              // compound like iterator
  CompoundIterator pl, ql;          // Create a compound iterator

  // Assign corresponding input values ---------------------------------------------
  std::complex<double> wkRconj(0., 0.);
  std::complex<double> exp_iwk_theji(0., 0.);    // function of m & theji

  // II - Global matrix - Alpha-Beta (n, m, l, k) - AlBe_nmlk ---------------------
  // prepare for computing transfer matrices ---------------------------------------
  // Matrices sizes ----------------------------------------------------------------
  // based on n_max
  int n_Matsize, m_Matsize;         // dependent on (n_max)
  n_Matsize = (n_max)+1;            // up to and including n_max    : indexed from 1
  m_Matsize = 2*(n_max)+1;          // up to and including m_max    : indexed from 0 + 1 for m==0
  // based on n_max+e             // (e) is required from recurrence relations
  int e=7;                  // required for extra values in 'AlBe_00lk'
  int n_Matsize1(0), m_Matsize1(0);     // dependent on (n_max+e)
  n_Matsize1 = (n_max+e)+1;         // up to and including (n_max+e)  : indexed from 1
  m_Matsize1 = 2*(n_max+e)+1;         // up to and including (n_max+e)  : indexed from 0 + 1 for m==0
  // scalar translation-addition theorem ------------------------------------------
  std::complex<double> ****AlBe_nmlk;
  AlBe_nmlk = Tools::Get_4D_c_double(n_Matsize1, m_Matsize1, n_Matsize1, m_Matsize1);
  PeriodicCoupling::compute_AlBe_nmlk(R, waveK, BHreg, n_max, AlBe_nmlk);

  // III - Global matrix - Anmlk & Bnmlk (n, m, l, k) -----------------------------
  // vector wave translation-addition theorem -------------------------------------
  std::complex<double> ****Anmlk, ****Bnmlk;
  Anmlk = Tools::Get_4D_c_double(n_Matsize1, m_Matsize1, n_Matsize1, m_Matsize1);
  Bnmlk = Tools::Get_4D_c_double(n_Matsize1, m_Matsize1, n_Matsize1, m_Matsize1);
  for(i=0; i<n_Matsize1; i++){
    for(j=0; j<m_Matsize1; j++){
      for(ii=0; ii<n_Matsize1; ii++){
        for(jj=0; jj<m_Matsize1; jj++){
          Anmlk[i][j][ii][jj]=std::complex<double>(0., 0.);
          Bnmlk[i][j][ii][jj]=std::complex<double>(0., 0.);
        }
      }
    }
  }

  // A - evaluate Anmlk & Bnmlk ------------------------------------------------------
  for(i=1, n=1; i<n_Matsize1-e; i++, n++){              // start at n==1 up to and including n_max
    for(j=e, m=-n_max; j<m_Matsize1-e; j++, m++){         // increment by the padded e value

      d_n=double(n);  d_m=double(m);

      if(std::abs(m)<=n){

      // for all m values -------------------------------------------------------
      if(1){

        for(ii=1, l=1; ii<n_Matsize1-e; ii++, l++){
          for(jj=e, k=-n_max; jj<m_Matsize1-e; jj++, k++){

            d_l=double(l);  d_k=double(k);

            if(std::abs(k) > l){
              Anmlk[i][j][ii][jj]=std::complex<double>(0., 0.);
              Bnmlk[i][j][ii][jj]=std::complex<double>(0., 0.);
            }

            else{
              // Anmlk ---------------------------------------------------
              // obtain three coefficients -------------------------------
              A1 = 0.5* sqrt( 1./(d_l*(d_l+1.)*d_n*(d_n+1.)) );
              A2 =    sqrt( (d_n-d_m)*(d_n+d_m+1.)*(d_l-d_k)*(d_l+d_k+1.) );
              A3 =    sqrt( (d_n+d_m)*(d_n-d_m+1.)*(d_l+d_k)*(d_l-d_k+1.) );

              Anmlk[i][j][ii][jj] = 2.*d_k*d_m*AlBe_nmlk[i][j][ii][jj]
                        + A2*AlBe_nmlk[i][j+1][ii][jj+1]
                        + A3*AlBe_nmlk[i][j-1][ii][jj-1];
              Anmlk[i][j][ii][jj]*=A1;

              // Bnmlk ---------------------------------------------------
              // obtain three coefficients -------------------------------
              B1 = 0.5* sqrt( ((2.*d_l+1.)/(2.*d_l-1.))*(1./(d_l*(d_l+1.)*d_n*(d_n+1.))) );
              B2 =    sqrt( (d_n-d_m)*(d_n+d_m+1.)*(d_l-d_k)*(d_l-d_k-1.) );
              B3 =    sqrt( (d_n+d_m)*(d_n-d_m+1.)*(d_l+d_k)*(d_l+d_k-1.) );
              Bnmlk[i][j][ii][jj] = 2.*d_m*sqrt((d_l-d_k)*(d_l+d_k))*AlBe_nmlk[i][j][ii-1][jj]
                          + B2*AlBe_nmlk[i][j+1][ii-1][jj+1]
                          - B3*AlBe_nmlk[i][j-1][ii-1][jj-1];
              Bnmlk[i][j][ii][jj]*=std::complex<double>(0., -B1);


            }
          } // jj
        } // ii

      } // if(m>=0 && m!=n)

      } // if(1)

    } // j
  } // i

  // B - ---------------------------------------------------------------------------------
  // store all Anmlk & Bnmlk into a compound iterator like 2D matrix ---------------------
  int n_check(0), m_check(0), l_check(0), k_check(0);
  //p=0;
  for(i=1, n=1; i<n_Matsize1-e; i++, n++){              // start at n==1 up to and including n_max
    for(j=e, m=-n_max; j<m_Matsize1-e; j++, m++){         // increment by the padded e value

      if(std::abs(m)<=n){

        // double check current n & m -----------------------------------------------
        pl.init(n,m);
        p=pl;
        n_check = int(sqrt(double(p)+1.));
        m_check = -(p+1)+n*(n+1);
        //q=0;

        // for all m values --------------------------------------------------------
        for(ii=1, l=1; ii<n_Matsize1-e; ii++, l++){
          for(jj=e, k=-n_max; jj<m_Matsize1-e; jj++, k++){

            if(std::abs(k)<=l){
              // double check current l & k -----------------------------------
              ql.init(l,k);
              q=ql;
              l_check = int(sqrt(double(q)+1.));
              k_check = -(q+1)+l*(l+1);
              // map 4D matrix --> 2D matrix ----------------------------------
              dataApq[p][q]=Anmlk[i][j][ii][jj];
              dataBpq[p][q]=Bnmlk[i][j][ii][jj];
              //q++;                    // increment q
            }

          } // jj
        } // ii

        //p++;                          // increment p

      } // if(abs(m)<=n)

    } // j
  } // i

  // C - delete all other matrices --------------------------------------------------------
  for(i=0; i<n_Matsize1; i++){
    for(j=0; j<m_Matsize1; j++){
      for(ii=0; ii<n_Matsize1; ii++){
        delete[] AlBe_nmlk[i][j][ii];
        delete[] Anmlk[i][j][ii];   delete[] Bnmlk[i][j][ii];
  } } }

  for(i=0; i<n_Matsize1; i++){
    for(j=0; j<m_Matsize1; j++){
        delete[] AlBe_nmlk[i][j];
        delete[] Anmlk[i][j];     delete[] Bnmlk[i][j];
  } }

  for(i=0; i<n_Matsize1; i++){
        delete[] AlBe_nmlk[i];
        delete[] Anmlk[i];        delete[] Bnmlk[i];
  }

  delete[] AlBe_nmlk;
  delete[] Anmlk;   delete[] Bnmlk;
}

Coupling::Coupling()
{
  initDone = false;
}

Coupling::Coupling(Spherical<double> relR_, std::complex<double> waveK_,
    int regular_, long nMax_)
{
  init(relR_, waveK_, regular_, nMax_);
}

Coupling::~Coupling()
{
  if(initDone)
  {
    for(int i=0; i<Tools::iteratorMax(nMax); i++)
    {
      delete [] dataApq[i];
      delete [] dataBpq[i];
    }
    delete [] dataApq;
    delete [] dataBpq;
  }
}

void Coupling::init(Spherical<double> relR_, std::complex<double> waveK_,
    int regular_, long nMax_)
{
  relR = relR_;
  waveK = waveK_;
  nMax = nMax_;

  regular = regular_;

  if(regular_ == 0)
    regular = 1;
  else
    regular = 0;

  if(initDone)
  {
    for(int i=0; i<Tools::iteratorMax(nMax); i++)
    {
      delete [] dataApq[i];
      delete [] dataBpq[i];
    }
    delete [] dataApq;
    delete [] dataBpq;

    initDone = false;
  }

  if(!initDone)
  {
    dataApq = new std::complex<double> *[Tools::iteratorMax(nMax)];
    dataBpq = new std::complex<double> *[Tools::iteratorMax(nMax)];

    for(int i=0; i<Tools::iteratorMax(nMax); i++)
    {
      dataApq[i] = new std::complex<double>[Tools::iteratorMax(nMax)];
      dataBpq[i] = new std::complex<double>[Tools::iteratorMax(nMax)];
    }
  }

  for(int i=0; i<Tools::iteratorMax(nMax); i++)
    for(int j=0; i<Tools::iteratorMax(nMax); i++)
    {
      dataApq[i][j] = std::complex<double>(0.0, 0.0);
      dataBpq[i][j] = std::complex<double>(0.0, 0.0);
    }

  initDone = true;
}

int Coupling::populate()
{
  if(!initDone)
  {
    std::cerr << "Coupling object not initialized.";
    return 1;
  }

  if(std::abs(relR.rrr) < errEpsilon) //Check for NO translation case
    for(int i=0; i<Tools::iteratorMax(nMax); i++)
      dataApq[i][i] = std::complex<double>(1.0, 0.0);
  else
    TransferCoefficients(relR, waveK, regular, nMax, dataApq, dataBpq);

  return 0;
}
