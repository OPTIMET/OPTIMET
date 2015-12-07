#include "Result.h"
#include "Coupling.h"
#include "Algebra.h"
#include "constants.h"
#include "Tools.h"
#include "AuxCoefficients.h"

#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>

using std::ofstream;

using std::cerr;
using std::cout;
using std::endl;
using std::real;
using std::conj;

Result::Result()
{
  initDone = false;
  flagSH = false;
  result_FF = NULL;
}

Result::Result(Geometry *geometry_, Excitation *excitation_, int nMax_)
{
  init(geometry_, excitation_, nMax_);
}

Result::Result(Geometry *geometry_, Excitation *excitation_, Result *result_FF_, int nMax_)
{
  init(geometry_, excitation_, result_FF_, nMax_);
}

Result::~Result()
{
  if(initDone)
  {
    delete[] scatter_coef;
    delete[] internal_coef;
    delete[] c_scatter_coef;
  }
}

void Result::init(Geometry *geometry_, Excitation *excitation_, int nMax_)
{
  if(initDone)
  {
    cerr << "Result object initialized previously! Use update()!";
    exit(1);
  }

  geometry = geometry_;
  nMax = nMax_;
  excitation = excitation_;
  waveK = excitation->waveK;
  flagSH = false;
  result_FF = NULL;

  scatter_coef = new complex<double>[2*Tools::iteratorMax(nMax)*geometry->noObjects];
  internal_coef = new complex<double>[2*Tools::iteratorMax(nMax)*geometry->noObjects];
  c_scatter_coef = new complex<double>[2*Tools::iteratorMax(nMax)];

  initDone = true;
}

void Result::update(Geometry *geometry_, Excitation *excitation_, int nMax_)
{
  geometry = geometry_;
  nMax = nMax_;
  excitation = excitation_;
  waveK = excitation->waveK;
}

void Result::init(Geometry *geometry_, Excitation *excitation_, Result *result_FF_, int nMax_)
{
  geometry = geometry_;
  nMax = nMax_;
  excitation = excitation_;
  waveK = excitation->waveK;
  flagSH = true;
  result_FF = result_FF_;

  if(!initDone)
  {
    scatter_coef = new complex<double>[2*Tools::iteratorMax(nMax)*geometry->noObjects];
    internal_coef = new complex<double>[2*Tools::iteratorMax(nMax)*geometry->noObjects];
    c_scatter_coef = new complex<double>[2*Tools::iteratorMax(nMax)];
  }

  initDone = true;
}

void Result::getEHFieldsModal(Spherical<double> R_, SphericalP<complex<double> > &EField_, SphericalP<complex<double> > &HField_, int projection_, CompoundIterator p_, int singleComponent_)
{
  SphericalP<complex<double> > Efield = SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
  SphericalP<complex<double> > Einc = SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
  SphericalP<complex<double> > Hfield = SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
  SphericalP<complex<double> > Hinc = SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));

  Spherical<double> Rrel;

  complex<double> iZ = (consCmi / sqrt(geometry->bground.mu/geometry->bground.epsilon));

  int intInd = geometry->checkInner(R_);

  if(intInd < 0) //Outside a sphere
  {
    if(!flagSH) //this a fundamental frequency result - calculate the incoming field
    {
      //Incoming field
      AuxCoefficients aCoefInc;
      aCoefInc.init(R_, waveK, 1, nMax);
      aCoefInc.populate();

      if(singleComponent_ == 1) //Only TE part
      {
        Einc = Einc + aCoefInc.dataMp[p_] * excitation->dataIncAp[p_];
        Hinc = Hinc + aCoefInc.dataMp[p_] * excitation->dataIncBp[p_] * iZ;
      }

      if(singleComponent_ == 2) //Only TM part
      {
        Einc = Einc + aCoefInc.dataNp[p_] * excitation->dataIncBp[p_];
        Hinc = Hinc + aCoefInc.dataNp[p_] * excitation->dataIncAp[p_] * iZ;
      }

      if(!singleComponent_) //Both TE and TM part
      {
        Einc = Einc + (aCoefInc.dataMp[p_] * excitation->dataIncAp[p_] + aCoefInc.dataNp[p_] * excitation->dataIncBp[p_]);
        Hinc = Hinc + (aCoefInc.dataNp[p_] * excitation->dataIncAp[p_] + aCoefInc.dataMp[p_] * excitation->dataIncBp[p_]) * iZ;
      }
    }
    else //this a second harmonic frequency result - calculate the source fields (save it in Einc for convenience)
    {
      //Source fields
      for(int j=0; j<geometry->noObjects; j++)
      {
        Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
        AuxCoefficients aCoef;
        aCoef.init(Rrel, waveK, 0, nMax);
        aCoef.populate();

        Einc = Einc + (aCoef.dataMp[p_] * geometry->objects[j].sourceCoef[p_] + aCoef.dataNp[p_] * geometry->objects[j].sourceCoef[p_ + p_.max(nMax)]);
      }
    }

    //Scattered field
    for(int j=0; j<geometry->noObjects; j++)
    {
      SphericalP<complex<double> > Efield_local = SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));

      Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
      AuxCoefficients aCoef;
      aCoef.init(Rrel, waveK, 0, nMax);
      aCoef.populate();

      if(singleComponent_ == 1) //TE Part only
      {
        Efield = Efield + aCoef.dataMp[p_]* scatter_coef[j*2*p_.max(nMax)+p_.compound];
        Hfield = Hfield + aCoef.dataMp[p_] * scatter_coef[p_.max(nMax)+j*2*p_.max(nMax)+p_.compound] * iZ;
      }

      if(singleComponent_ == 2) //TM Part only
      {
        Efield = Efield + aCoef.dataNp[p_] * scatter_coef[p_.max(nMax)+j*2*p_.max(nMax)+p_.compound];
        Hfield = Hfield + aCoef.dataNp[p_]* scatter_coef[j*2*p_.max(nMax)+p_.compound] * iZ;
      }

      if(!singleComponent_)
      {
        Efield = Efield + aCoef.dataMp[p_]* scatter_coef[j*2*p_.max(nMax)+p_.compound] + aCoef.dataNp[p_] * scatter_coef[p_.max(nMax)+j*2*p_.max(nMax)+p_.compound];
        Hfield = Hfield + (aCoef.dataNp[p_]* scatter_coef[j*2*p_.max(nMax)+p_.compound]  + aCoef.dataMp[p_] * scatter_coef[p_.max(nMax)+j*2*p_.max(nMax)+p_.compound]) * iZ;
      }
    }
  }
  else //Inside a sphere
  {
    Rrel = Tools::toPoint(R_, geometry->objects[intInd].vR);
    AuxCoefficients aCoef;
    aCoef.init(Rrel, waveK * sqrt(geometry->objects[intInd].elmag.epsilon_r * geometry->objects[intInd].elmag.mu_r), 1, nMax);
    aCoef.populate();

    complex<double> iZ_object = (consCmi / sqrt(geometry->objects[intInd].elmag.mu / geometry->objects[intInd].elmag.epsilon));

    if(singleComponent_ == 1) //TE Part Only
    {
      Efield = Efield + aCoef.dataMp[p_] * internal_coef[intInd*2*p_.max(nMax)+p_.compound];
      Hfield = Hfield + aCoef.dataMp[p_] * internal_coef[p_.max(nMax)+intInd*2*p_.max(nMax)+p_.compound] * iZ_object;
    }

    if(singleComponent_ == 2) //TM Part Only
    {
      Efield = Efield + aCoef.dataNp[p_] * internal_coef[p_.max(nMax)+intInd*2*p_.max(nMax)+p_.compound];
      Hfield = Hfield + aCoef.dataNp[p_] * internal_coef[intInd*2*p_.max(nMax)+p_.compound] * iZ_object;
    }

    if(!singleComponent_) //Both TE and TM
    {
      Efield = Efield + (aCoef.dataMp[p_] * internal_coef[intInd*2*p_.max(nMax)+p_.compound] + aCoef.dataNp[p_] * internal_coef[p_.max(nMax)+intInd*2*p_.max(nMax)+p_.compound]);
      Hfield = Hfield + (aCoef.dataNp[p_] * internal_coef[intInd*2*p_.max(nMax)+p_.compound] + aCoef.dataMp[p_] * internal_coef[p_.max(nMax)+intInd*2*p_.max(nMax)+p_.compound]) * iZ_object;
    }
  }

  if(projection_)
  {
    SphericalP<complex<double> > SphEField;
    SphericalP<complex<double> > SphHField;
    Rrel = Tools::toPoint(R_, geometry->objects[0].vR);

    SphEField = Tools::fromProjection(Rrel, Einc+Efield);
    SphHField = Tools::fromProjection(Rrel, Hinc+Hfield);

    EField_ = SphEField;
    HField_ = SphHField;
  }
  else
  {
    EField_ = Einc + Efield;
    HField_ = Hinc + Hfield;
  }
}

void Result::getEHFields(Spherical<double> R_, SphericalP<complex<double> > &EField_, SphericalP<complex<double> > &HField_, int projection_)
{
  SphericalP<complex<double> > Efield = SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
  SphericalP<complex<double> > Einc = SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
  SphericalP<complex<double> > Hfield = SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
  SphericalP<complex<double> > Hinc = SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));

  Spherical<double> Rrel;

  complex<double> iZ = (consCmi / sqrt(geometry->bground.mu/geometry->bground.epsilon));

  int pMax = Tools::iteratorMax(nMax);

  CompoundIterator p;

  //Check for inner point and set to 0
  int intInd = geometry->checkInner(R_);

  if(intInd < 0) //Outside a sphere
  {
    if(!flagSH) //this a fundamental frequency result - calculate the incoming field
    {
      //Incoming field
      AuxCoefficients aCoefInc;
      aCoefInc.init(R_, waveK, 1, nMax);
      aCoefInc.populate();
      for(p=0; p<p.max(nMax); p++)
      {
        Einc = Einc + (aCoefInc.dataMp[p] * excitation->dataIncAp[p] + aCoefInc.dataNp[p] * excitation->dataIncBp[p]);
        Hinc = Hinc + (aCoefInc.dataNp[p] * excitation->dataIncAp[p] + aCoefInc.dataMp[p] * excitation->dataIncBp[p]) * iZ;
      }
    }
    else //this a second harmonic frequency result - calculate the source fields (save it in Einc for convenience)
    {
      //Source fields
      for(int j=0; j<geometry->noObjects; j++)
      {
        Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
        AuxCoefficients aCoef;
        aCoef.init(Rrel, waveK, 0, nMax);
        aCoef.populate();

        for(p=0; p<pMax; p++)
        {
          Einc = Einc + (aCoef.dataMp[p] * geometry->objects[j].sourceCoef[p] + aCoef.dataNp[p] * geometry->objects[j].sourceCoef[p + pMax]);
        }
      }
    }

    //Scattered field
    for(int j=0; j<geometry->noObjects; j++)
    {
      SphericalP<complex<double> > Efield_local = SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));

      Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
      AuxCoefficients aCoef;
      aCoef.init(Rrel, waveK, 0, nMax);
      aCoef.populate();

      for(p=0; p<p.max(nMax); p++)
      {
        Efield = Efield + aCoef.dataMp[p]* scatter_coef[j*2*pMax+p.compound] + aCoef.dataNp[p] * scatter_coef[pMax+j*2*pMax+p.compound];
        Hfield = Hfield + (aCoef.dataNp[p]* scatter_coef[j*2*pMax+p.compound]  + aCoef.dataMp[p] * scatter_coef[pMax+j*2*pMax+p.compound]) * iZ;
      }
    }
  }
  else //Inside a sphere
  {
    Rrel = Tools::toPoint(R_, geometry->objects[intInd].vR);
    AuxCoefficients aCoef;
    aCoef.init(Rrel, waveK * sqrt(geometry->objects[intInd].elmag.epsilon_r * geometry->objects[intInd].elmag.mu_r), 1, nMax);
    aCoef.populate();

    complex<double> iZ_object = (consCmi / sqrt(geometry->objects[intInd].elmag.mu / geometry->objects[intInd].elmag.epsilon));

    for(p=0; p<p.max(nMax); p++)
    {
      Efield = Efield + (aCoef.dataMp[p] * internal_coef[intInd*2*pMax+p.compound] + aCoef.dataNp[p] * internal_coef[pMax+intInd*2*pMax+p.compound]);
      Hfield = Hfield + (aCoef.dataNp[p] * internal_coef[intInd*2*pMax+p.compound] + aCoef.dataMp[p] * internal_coef[pMax+intInd*2*pMax+p.compound]) * iZ_object;
    }
  }

  if(projection_)
  {
    SphericalP<complex<double> > SphEField;
    SphericalP<complex<double> > SphHField;
    Rrel = Tools::toPoint(R_, geometry->objects[0].vR);

    SphEField = Tools::fromProjection(Rrel, Einc+Efield);
    SphHField = Tools::fromProjection(Rrel, Hinc+Hfield);

    EField_ = SphEField;
    HField_ = SphHField;
  }
  else
  {
    EField_ = Einc + Efield;
    HField_ = Hinc + Hfield;
  }
}

double Result::getExtinctionCrossSection()
{
  CompoundIterator p;
  int pMax = Tools::iteratorMax(nMax);

  double Cext(0.);
  complex<double> *Q_local = new complex<double>[2*pMax];

  for(int j=0; j<geometry->noObjects; j++)
  {
    excitation->getIncLocal(geometry->objects[j].vR, Q_local, nMax);
    for(p=0; p<pMax; p++)
    {
      Cext += real( conj(Q_local[p]) * scatter_coef[j*2*pMax+p.compound]
             +  conj(Q_local[p.compound + pMax]) * scatter_coef[pMax+j*2*pMax+p.compound]);
    }

  }

  delete [] Q_local;
  return (-1. / (real(waveK) * real(waveK)))*Cext;
}

double Result::getAbsorptionCrossSection()
{
  CompoundIterator p;
  int pMax = Tools::iteratorMax(nMax);

  double Cabs(0.);
  double temp1(0.), temp2(0.);
  double *Cabs_aux = new double[2*pMax];

  for(int j=0; j<geometry->noObjects; j++)
  {

    geometry->getCabsAux(excitation->omega, j, nMax, Cabs_aux);

    for(p=0; p<pMax; p++)
    {
      temp1 = abs(scatter_coef[j*2*pMax+p.compound]);
      temp1*=temp1;
      temp2 = abs(scatter_coef[pMax+j*2*pMax+p.compound]);
      temp2*=temp2;
      Cabs +=   temp1 * Cabs_aux[p.compound]
          +   temp2 * Cabs_aux[pMax+p.compound];
    }
  }
  delete [] Cabs_aux;
  return ( 1 / (real(waveK) * real(waveK))) * Cabs;
}

int Result::setFields(OutputGrid oEGrid_, OutputGrid oHGrid_, int projection_)
{
  if(!initDone)
  {
    cerr << "Result object not initialized!";
    return -1;
  }

  Spherical<double> Rloc;

  //centerScattering();

  //Calculate the fields
  while(!oEGrid_.gridDone)
  {
    Rloc = oEGrid_.getPoint();
    oHGrid_.getPoint();
    cout << "Calculating fields for point " << oEGrid_.iterator+1 << " out of " << oEGrid_.gridPoints << endl;

    SphericalP<complex<double> > EField;
    SphericalP<complex<double> > HField;

    getEHFields(Rloc, EField, HField, projection_);

    oHGrid_.pushDataNext(HField);
    oEGrid_.pushDataNext(EField);

  }

  return 0;
}

int Result::setFieldsModal(OutputGrid oEGrid_, OutputGrid oHGrid_, int projection_, CompoundIterator p_, int singleComponent_)
{
  if(!initDone)
  {
    cerr << "Result object not initialized!";
    return -1;
  }

  Spherical<double> Rloc;

  //Calculate the fields
  while(!oEGrid_.gridDone)
  {
    Rloc = oEGrid_.getPoint();
    oHGrid_.getPoint();
    cout << "Calculating fields for point " << oEGrid_.iterator+1 << " out of " << oEGrid_.gridPoints << endl;

    SphericalP<complex<double> > EField;
    SphericalP<complex<double> > HField;

    getEHFieldsModal(Rloc, EField, HField, projection_, p_, singleComponent_);

    oHGrid_.pushDataNext(HField);
    oEGrid_.pushDataNext(EField);
  }
  return 0;
}

void Result::centerScattering()
{
  CompoundIterator p,q;

  int pMax = p.max(nMax);
  int qMax = q.max(nMax);

  for(p=0; p<2*p.max(nMax); p++)
  {
    c_scatter_coef[p] = complex<double>(0.0, 0.0);
  }

  complex<double> **T_AB = new complex<double>*[2*(p.max(nMax))];
  complex<double> *scatter_aux = new complex<double>[2*p.max(nMax)];
  complex<double> *scatter_fin = new complex<double>[2*p.max(nMax)];

  for(p=0; p<(int)(2*p.max(nMax)); p++)
  {
    T_AB[p] = new complex<double>[2*p.max(nMax)];
  }


  for(int i=0; i<geometry->noObjects; i++)
  {
    Spherical<double> Rrel = Tools::toPoint(Spherical<double>(0.0, 0.0, 0.0), geometry->objects[i].vR);

    Coupling coupling;
    coupling.init(Rrel, excitation->waveK, 0, nMax);
    coupling.populate();

    for(p=0; p<pMax; p++)
      for(q=0; q<qMax; q++)
      {
        T_AB[p][q] = coupling.dataApq[p][q];
        T_AB[p+pMax][q+qMax] = coupling.dataApq[p][q];
        T_AB[p+pMax][q] = coupling.dataBpq[p][q];
        T_AB[p][q+qMax] = coupling.dataBpq[p][q];
      }

    for(p=0; p<2*p.max(nMax); p++)
    {
      scatter_aux[p] += scatter_coef[p.compound + i*2*pMax*geometry->noObjects];
    }

    Algebra::multiplyVectorMatrix(T_AB, 2*pMax, 2*pMax, scatter_aux, scatter_fin, consC1, consC0);

    for(p=0; p<2*pMax; p++)
    {
      c_scatter_coef[p] += scatter_fin[p];
    }
  }


  for(p=0; p<(int)(2*p.max(nMax)); p++)
  {
    delete [] T_AB[p];
  }

  delete [] T_AB;
  delete [] scatter_aux;
  delete [] scatter_fin;
}

CompoundIterator Result::getDominant()
{
  CompoundIterator p, q;

  q = 0;

  complex<double> TEMax = scatter_coef[0];
  complex<double> TMMax = scatter_coef[p.max(nMax)];

  for(p=0; p<p.max(nMax); p++)
  {
    if((abs(scatter_coef[p]) > abs(TEMax)) || (abs(scatter_coef[p.max(nMax) + p.compound]) > abs(TMMax)))
    {
      q = p;

      TEMax = scatter_coef[p];
      TMMax = scatter_coef[p.max(nMax) + p.compound];
    }
  }

  return q;
}

void Result::getEHFieldsContCheck(Spherical<double> R_, SphericalP<complex<double> > &EField_, SphericalP<complex<double> > &HField_, int projection_, int inside_)
{
  SphericalP<complex<double> > Efield = SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
  SphericalP<complex<double> > Einc = SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
  SphericalP<complex<double> > Hfield = SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
  SphericalP<complex<double> > Hinc = SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));

  Spherical<double> Rrel;

  complex<double> iZ = (consCmi / sqrt(geometry->bground.mu/geometry->bground.epsilon));

  int pMax = Tools::iteratorMax(nMax);

  CompoundIterator p;

  //Check for inner point and set to 0
  int intInd = geometry->checkInner(R_);
  intInd = inside_;

  if(intInd < 0) //Outside a sphere
  {
    if(!flagSH) //this a fundamental frequency result - calculate the incoming field
    {
      //Incoming field
      AuxCoefficients aCoefInc;
      aCoefInc.init(R_, waveK, 1, nMax);
      aCoefInc.populate();
      for(p=0; p<p.max(nMax); p++)
      {
        Einc = Einc + (aCoefInc.dataMp[p] * excitation->dataIncAp[p] + aCoefInc.dataNp[p] * excitation->dataIncBp[p]);
        Hinc = Hinc + (aCoefInc.dataNp[p] * excitation->dataIncAp[p] + aCoefInc.dataMp[p] * excitation->dataIncBp[p]) * iZ;
      }
    }
    else //this a second harmonic frequency result - calculate the source fields (save it in Einc for convenience)
    {
      //Source fields
      for(int j=0; j<geometry->noObjects; j++)
      {
        Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
        AuxCoefficients aCoef;
        aCoef.init(Rrel, waveK, 0, nMax);
        aCoef.populate();

        for(p=0; p<pMax; p++)
        {
          Einc = Einc + (aCoef.dataMp[p] * geometry->objects[j].sourceCoef[p] + aCoef.dataNp[p] * geometry->objects[j].sourceCoef[p + pMax]);
        }
      }
    }

    //Scattered field
    for(int j=0; j<geometry->noObjects; j++)
    {
      SphericalP<complex<double> > Efield_local = SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));

      Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
      AuxCoefficients aCoef;
      aCoef.init(Rrel, waveK, 0, nMax);
      aCoef.populate();

      for(p=0; p<p.max(nMax); p++)
      {
        Efield = Efield + aCoef.dataMp[p]* scatter_coef[j*2*pMax+p.compound] + aCoef.dataNp[p] * scatter_coef[pMax+j*2*pMax+p.compound];
        Hfield = Hfield + (aCoef.dataNp[p]* scatter_coef[j*2*pMax+p.compound]  + aCoef.dataMp[p] * scatter_coef[pMax+j*2*pMax+p.compound]) * iZ;
      }
    }
  }
  else //Inside a sphere
  {
    Rrel = Tools::toPoint(R_, geometry->objects[intInd].vR);
    AuxCoefficients aCoef;
    aCoef.init(Rrel, waveK * sqrt(geometry->objects[intInd].elmag.epsilon_r * geometry->objects[intInd].elmag.mu_r), 1, nMax);
    aCoef.populate();

    complex<double> iZ_object = (consCmi / sqrt(geometry->objects[intInd].elmag.mu / geometry->objects[intInd].elmag.epsilon));

    for(p=0; p<p.max(nMax); p++)
    {
      Efield = Efield + (aCoef.dataMp[p] * internal_coef[intInd*2*pMax+p.compound] + aCoef.dataNp[p] * internal_coef[pMax+intInd*2*pMax+p.compound]);
      Hfield = Hfield + (aCoef.dataNp[p] * internal_coef[intInd*2*pMax+p.compound] + aCoef.dataMp[p] * internal_coef[pMax+intInd*2*pMax+p.compound]) * iZ_object;
    }
  }

  if(projection_)
  {
    SphericalP<complex<double> > SphEField;
    SphericalP<complex<double> > SphHField;
    Rrel = Tools::toPoint(R_, geometry->objects[0].vR);

    SphEField = Tools::fromProjection(Rrel, Einc+Efield);
    SphHField = Tools::fromProjection(Rrel, Hinc+Hfield);

    EField_ = SphEField;
    HField_ = SphHField;
  }

  else{ // AJ - no spherical projection
    EField_ = Einc + Efield;
    HField_ = Hinc + Hfield;
  }
}

void Result::writeContinuityCheck(int objectIndex_)
{
    SphericalP<complex<double> > AnEField_in, AnEField_out;
    SphericalP<complex<double> > AnHField_in, AnHField_out;
    Spherical<double> APoint(0.0, 0.0, 0.0);
    int projection = 1;     //Spherical projection - True - projection is internally set to be evaluated w.r.t. object[0]
    int outside = -1;     //Forces result to be outside an object
    int inside = 0;       //Forces result to be inside an object
    double radius       = geometry->objects[objectIndex_].radius;;
    complex<double> eps_r = geometry->objects[objectIndex_].elmag.epsilon_r;
    complex<double> mu_r  = geometry->objects[objectIndex_].elmag.mu_r;

    // the-phi - 2D plot ----------------------------------------------------------------------------------------------------------------
    ofstream E1_err_mag("E1_err_mag");
    ofstream E2_err_mag("E2_err_mag");
    ofstream E3_err_mag("E3_err_mag");
    ofstream H1_err_mag("H1_err_mag");
    ofstream H2_err_mag("H2_err_mag");
    ofstream H3_err_mag("H3_err_mag");
    int max_ii = 180;         // theta observation range (increment by 1 degree) - [1, max_ii-1]
    int max_jj = 180;         // phi   observation range (increment by 1 degree) - [1, max_jj-1]
    for(int ii=1; ii<=max_ii-1; ii++){
      for(int jj=1; jj<=max_jj-1; jj++){

        APoint=Spherical<double>(radius, consPi*(double(ii)/180.), consPi*(double(jj)/180.));
        getEHFieldsContCheck(APoint, AnEField_out, AnHField_out, projection, outside);
        getEHFieldsContCheck(APoint, AnEField_in, AnHField_in, projection, inside);
        // op -------------------------------------------------------------------
        cout<<"Continuity check : computed "<<jj+((ii-1)*(max_ii-1))<<" out of a total of "<<(max_ii-1)*(max_jj-1)<<endl;
        // EF
        E1_err_mag<<(abs(AnEField_out.rrr)-abs(AnEField_in.rrr*eps_r))/abs(AnEField_in.rrr*eps_r)<<" ";
        E2_err_mag<<(abs(AnEField_out.the)-abs(AnEField_in.the))/abs(AnEField_in.the)<<" ";
        E3_err_mag<<(abs(AnEField_out.phi)-abs(AnEField_in.phi))/abs(AnEField_in.phi)<<" ";
        // EF
        H1_err_mag<<(abs(AnHField_out.rrr)-abs(AnHField_in.rrr*mu_r))/abs(AnHField_in.rrr*mu_r)<<" ";
        H2_err_mag<<(abs(AnHField_out.the)-abs(AnHField_in.the))/abs(AnHField_in.the)<<" ";
        H3_err_mag<<(abs(AnHField_out.phi)-abs(AnHField_in.phi))/abs(AnHField_in.phi)<<" ";
      }
      E1_err_mag<<endl; E2_err_mag<<endl; E3_err_mag<<endl;
      H1_err_mag<<endl; H2_err_mag<<endl; H3_err_mag<<endl;
    }
    E1_err_mag.flush(); E2_err_mag.flush(); E3_err_mag.flush();
    H1_err_mag.flush(); H2_err_mag.flush(); H3_err_mag.flush();
    E1_err_mag.close(); E2_err_mag.close(); E3_err_mag.close();
    H1_err_mag.close(); H2_err_mag.close(); H3_err_mag.close();
}
