#include "Solver.h"

#include "CompoundIterator.h"
#include "Tools.h"
#include "AlgebraS.h"
#include "Algebra.h"
#include "constants.h"
#include "aliases.h"
#include <iostream>
#include <cstdlib>

using std::cerr;
using std::cout;
using std::endl;

Solver::Solver()
{
  initDone = false;
  flagSH = false;
}

Solver::Solver(Geometry* geometry_, Excitation* incWave_, int method_, long nMax_)
{
  init(geometry_, incWave_, method_, nMax_);
}

Solver::~Solver()
{
  if(initDone)
  {
    CompoundIterator p;

    for(p=0; p<(int)(2*p.max(nMax)* geometry->noObjects); p++)
    {
      if(S[p])
        delete [] S[p];
    }

    delete [] S;

    delete [] Q;
  }
}

void Solver::init(Geometry* geometry_, Excitation* incWave_, int method_, long nMax_)
{
  if(initDone)
  {
    cerr << "Solver object initialized previously! Use update()!";
    exit(1);
  }

  geometry = geometry_;
  incWave = incWave_;
  nMax = nMax_;

  flagSH = false;

  solverMethod = method_;

  CompoundIterator p;

  S = new complex<double>* [2 * p.max(nMax) * geometry->noObjects];
  for(p=0; p<(int)(2*p.max(nMax)* geometry->noObjects); p++)
  {
    S[p] = new complex<double>[2 * p.max(nMax) * geometry->noObjects];
  }

  Q = new complex<double>[2 * p.max(nMax) * geometry->noObjects];

  result_FF = NULL;

  initDone = true;
}

int Solver::populate()
{
  if(!initDone)
  {
    cerr << "Solver not initialized!";
    return 1;
  }

  if(solverMethod == O3DSolverDirect)
    populateDirect();
  else
    if(solverMethod == O3DSolverIndirect)
      populateIndirect();
    else //Default
      populateDirect();

  return 0;
}

int Solver::populateDirect()
{
  if(!initDone)
  {
    cerr << "Solver not initialized!";
    return 1;
  }

  CompoundIterator p;
  CompoundIterator q;

  //Local matrices for storing T, T_fin = -T*A and T_AB = AB to be multiplied by BLAS.
  //Also local matrix Q_local to be added to S.

  complex<double> **T = new complex<double>*[2*p.max(nMax)];
  complex<double> **Ti = new complex<double>*[2*p.max(nMax)];
  complex<double> **T_fin = new complex<double>*[2*p.max(nMax)];
  complex<double> **T_AB = new complex<double>*[2*(p.max(nMax))];

  for(p=0; p<(int)(2*p.max(nMax)); p++)
  {
    T[p] = new complex<double>[2*p.max(nMax)];
    Ti[p] = new complex<double>[2*p.max(nMax)];
    T_fin[p] = new complex<double>[2*p.max(nMax)];
    T_AB[p] = new complex<double>[2*p.max(nMax)];
  }

  if(flagSH) //if SH simulation, set the local source first
  {
    geometry->setSourcesSingle(incWave, result_FF->internal_coef, nMax);
  }

  int i, j;
  for(i=0; i<geometry->noObjects; i++)
  {
    int pMax = p.max(nMax);
    int qMax = q.max(nMax);

    complex<double> *Q_local = new complex<double>[2*pMax];

    //Get the T and IncLocal matrices first
    geometry->getTLocal(incWave->omega, i, nMax, T);
    if(flagSH) //we are in the SH case -> get the local sources from the geometry
    {
      geometry->getSourceLocal(i, incWave, result_FF->internal_coef, nMax, Q_local);
    }
    else //we are in the FF case -> get the incoming excitation from the geometry
    {
      incWave->getIncLocal(geometry->objects[i].vR, Q_local, nMax);
    }

    //Allocate space for a "result" matrix C, multiply T [A] with Q_local and push to Q.
    complex<double> *C = new complex<double>[2*pMax];
    Algebra::multiplyVectorMatrix(T, 2*pMax, 2*pMax, Q_local, C, consC1, consC0);

    for(j=0; j<2*pMax; j++)
    {
      Q[j+i*2*pMax] = C[j];
    }

    delete [] Q_local;
    delete [] C;


    for(j=0; j<geometry->noObjects; j++)
    {
      if(i == j)
      {
        //Build and push an I matrix
        Tools::makeUnitMatrix(2*pMax, Ti);
        Tools::pushToMatrix(Ti, 2*pMax, 2*pMax, S, i*2*pMax, i*2*pMax);
      }
      else
      {
        //Build the T_AB matrix
        AB.init(geometry->objects[i].vR - geometry->objects[j].vR, incWave->waveK, 0, nMax);
        AB.populate();

        for(p=0; p<pMax; p++)
          for(q=0; q<qMax; q++)
          {
            T_AB[p][q] = AB.dataApq[q][p];
            T_AB[p+pMax][q+qMax] = AB.dataApq[q][p];
            T_AB[p+pMax][q] = AB.dataBpq[q][p];
            T_AB[p][q+qMax] = AB.dataBpq[q][p];
          }

        //Multiply T and T_AB to get T_fin = -T*T_AB
        Algebra::multiplyMatrixMatrix(T, 2*pMax, 2*pMax, T_AB, 2*pMax, 2*pMax, T_fin, consCm1, consC0);

          //Push T into the S matrix
        Tools::pushToMatrix(T_fin, 2*pMax, 2*pMax, S, i*2*pMax, j*2*pMax);
      }
    }
  }

  //Clean up memory

  for(p=0; p<(int)(2*p.max(nMax)); p++)
  {
    delete [] T[p];
    delete [] Ti[p];
    delete [] T_fin[p];
    delete [] T_AB[p];
  }

  delete [] T;
  delete [] Ti;
  delete [] T_fin;
  delete [] T_AB;

  return 0;
}

int Solver::solve(complex<double> *X_sca_, complex<double> *X_int_)
{
  if(solverMethod == O3DSolverDirect)
    solveScatteredDirect(X_sca_);
  else
    if(solverMethod == O3DSolverIndirect)
      solveScatteredIndirect(X_sca_);
    else //Default
      solveScatteredDirect(X_sca_);

  solveInternal(X_sca_, X_int_);

  return 0;
}

int Solver::solveScatteredDirect(complex<double>* X_sca_)
{
  if(!initDone)
  {
    cerr << "Solver object not initialized!";
    return 1;
  }

  AlgebraS::solveMatrixVector(S, 2*Tools::iteratorMax(nMax)*geometry->noObjects, 2*Tools::iteratorMax(nMax)*geometry->noObjects,
        Q, X_sca_);

  return 0;
}

int Solver::solveScatteredIndirect(complex<double>* X_sca_)
{
  if(!initDone)
  {
    cerr << "Solver object not initialized!";
    return 1;
  }

  CompoundIterator p;

  complex<double> *X_sca_local = new complex<double>[2*p.max(nMax)*geometry->noObjects];
  complex<double> *X_sca_part = new complex<double>[2*p.max(nMax)*geometry->noObjects];
  complex<double> *X_sca_part_fin = new complex<double>[2*p.max(nMax)*geometry->noObjects];
  complex<double> **T = new complex<double>*[2*p.max(nMax)];

  for(p=0; p<(int)(2*p.max(nMax)); p++)
  {
    T[p] = new complex<double>[2*p.max(nMax)];
  }

  //Solve the equation, here Q and S correspond to Eq. 10 in (Stout2002). Store result in X_sca_local
  AlgebraS::solveMatrixVector(S, 2*Tools::iteratorMax(nMax)*geometry->noObjects, 2*Tools::iteratorMax(nMax)*geometry->noObjects,
        Q, X_sca_local);

  for(int i=0; i<geometry->noObjects; i++)
  {
    //Get the local scattering matrix for object i
    geometry->getTLocal(incWave->omega, i, nMax, T);

    //Build a partial X_sca_part vector
    for(p=0; p<(int)(2*p.max(nMax)); p++)
    {
      X_sca_part[p] = X_sca_local[i*2*p.max(nMax) + p.compound];
    }

    //Multiply T with X_sca_part to get X_sca_fin
    Algebra::multiplyVectorMatrix(T, 2*p.max(nMax), 2*p.max(nMax), X_sca_part, X_sca_part_fin, consC1, consC0);

    //Push new X_sca_fin into the X_sca_ final solution
    for(p=0; p<(int)(2*p.max(nMax)); p++)
    {
      X_sca_[i*2*p.max(nMax) + p.compound] = X_sca_part_fin[p];
    }

  }

  //Clean up the memory
  delete [] X_sca_local;
  delete [] X_sca_part;
  delete [] X_sca_part_fin;
  for(p=0; p<(int)(2*p.max(nMax)); p++)
  {
    delete [] T[p];
  }
  delete [] T;

  return 0;
}

int Solver::switchSH(Excitation* incWave_, Result* result_FF_, long nMax_)
{
  if(!initDone)
  {
    cerr << "Solver object was not created and initialized in the FF case!";
    return 1;
  }

  incWave = incWave_;
  nMax = nMax_;

  result_FF = result_FF_;

  flagSH = true;

  return 0;
}

int Solver::solveInternal(complex<double>* X_sca_, complex<double> *X_int_)
{
  if(!initDone)
  {
    cerr << "Solver object not initialized!";
    return 1;
  }

  CompoundIterator p,q;

  complex<double> *Iaux = new complex<double>[2*Tools::iteratorMax(nMax)];


  for(int j=0; j<geometry->noObjects; j++)
  {
    int pMax = p.max(nMax);
    geometry->getIaux(incWave->omega, j, nMax, Iaux);

    for(p=0; p<pMax; p++)
    {
      // scattered
      X_int_[j*2*pMax+p.compound] = X_sca_[j*2*pMax+p.compound] * Iaux[p];
      X_int_[pMax+j*2*pMax+p.compound] = X_sca_[pMax+j*2*pMax+p.compound] * Iaux[p+pMax];
    }
  }

  delete[] Iaux;

  return 0;
}

int Solver::populateIndirect()
{
  if(!initDone)
  {
    cerr << "Geometry not initialized!";
    return 1;
  }

  CompoundIterator p;
  CompoundIterator q;

  //Local matrices for storing T, T_fin = -T*A and T_AB = AB to be multiplied by BLAS.
  //Also local matrix Q_local to be added to S.

  complex<double> **T = new complex<double>*[2*p.max(nMax)];
  complex<double> **Ti = new complex<double>*[2*p.max(nMax)];
  complex<double> **T_fin = new complex<double>*[2*p.max(nMax)];
  complex<double> **T_AB = new complex<double>*[2*(p.max(nMax))];

  for(p=0; p<(int)(2*p.max(nMax)); p++)
  {
    T[p] = new complex<double>[2*p.max(nMax)];
    Ti[p] = new complex<double>[2*p.max(nMax)];
    T_fin[p] = new complex<double>[2*p.max(nMax)];
    T_AB[p] = new complex<double>[2*p.max(nMax)];
  }

  if(flagSH) //if SH simulation, set the local source first
  {
    geometry->setSourcesSingle(incWave, result_FF->internal_coef, nMax);
  }

  int i, j;
  for(i=0; i<geometry->noObjects; i++)
  {
    int pMax = p.max(nMax);
    int qMax = q.max(nMax);

    complex<double> *Q_local = new complex<double>[2*pMax];

    //Get the IncLocal matrices
    //These correspond directly to the Beta*a in Stout2002 Eq. 10 as
    // they are already translated.
    if(flagSH) //we are in the SH case -> get the local sources from the geometry
    {
      geometry->getSourceLocal(i, incWave, result_FF->internal_coef, nMax, Q_local);
    }
    else //we are in the FF case -> get the incoming excitation from the geometry
    {
      incWave->getIncLocal(geometry->objects[i].vR, Q_local, nMax);
    }

    //Push Q_local into Q (direct equivalence) (NOTE: j is unused at this point)
    for(j=0; j<2*pMax; j++)
    {
      Q[j+i*2*pMax] = Q_local[j];
    }

    delete [] Q_local;

    for(j=0; j<geometry->noObjects; j++)
    {
      if(i == j)
      {
        //Build and push an I matrix
        Tools::makeUnitMatrix(2*pMax, Ti);
        Tools::pushToMatrix(Ti, 2*pMax, 2*pMax, S, i*2*pMax, i*2*pMax);
      }
      else
      {
        //Build the T_AB matrix (non-regular corresponding to alpha(i,j))
        AB.init(geometry->objects[i].vR - geometry->objects[j].vR, incWave->waveK, 0, nMax);
        AB.populate();

        for(p=0; p<pMax; p++)
          for(q=0; q<qMax; q++)
          {
            T_AB[p][q] = AB.dataApq[q][p];
            T_AB[p+pMax][q+qMax] = AB.dataApq[q][p];
            T_AB[p+pMax][q] = AB.dataBpq[q][p];
            T_AB[p][q+qMax] = AB.dataBpq[q][p];
          }

        //In this case, T needs to correspond to the SECOND object; so get it here
        geometry->getTLocal(incWave->omega, j, nMax, T);

        //Multiply T_AB and T to get T_fin = -T_AB*T
        Algebra::multiplyMatrixMatrix(T_AB, 2*pMax, 2*pMax, T, 2*pMax, 2*pMax, T_fin, consCm1, consC0);

        //Push T_fin into the S matrix
        Tools::pushToMatrix(T_fin, 2*pMax, 2*pMax, S, i*2*pMax, j*2*pMax);
      }
    }
  }

  //Clean up memory

  for(p=0; p<(int)(2*p.max(nMax)); p++)
  {
    delete [] T[p];
    delete [] Ti[p];
    delete [] T_fin[p];
    delete [] T_AB[p];
  }

  delete [] T;
  delete [] Ti;
  delete [] T_fin;
  delete [] T_AB;

  return 0;
}

void Solver::update(Geometry* geometry_, Excitation* incWave_, long nMax_)
{
  geometry = geometry_;
  incWave = incWave_;
  nMax = nMax_;

  result_FF = NULL;

  populate();
}
