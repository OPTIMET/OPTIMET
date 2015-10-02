#include "Legendre.h"

using std::cerr;
using std::pow;
using std::abs;

Legendre::Legendre()
{
	initDone = false;
}

Legendre::Legendre(double argument_, int nMax_)
{
	init(argument_, nMax_);
}

Legendre::~Legendre()
{
	if(initDone)
	{
		delete [] data;
		delete [] ddata;
	}
}

void Legendre::init(double argument_, int nMax_)
{
	argument = argument_;
	nMax = nMax_;
	data = new double[Tools::iteratorMax(nMax)];
	ddata = new double[Tools::iteratorMax(nMax)];
	initDone = true;
}

int Legendre::populate()
{
	if(!initDone)
	{
		cerr << "Legendre engine was not initialized!";
		return 1;
	}

	CompoundIterator q;

	for(int m=0; m<=nMax; m++)
	{
		int n;
		//Compute positives
		double *data_l = new double [nMax-m+1];
		double *ddata_l = new double [nMax-m+1];

		gsl_sf_legendre_Plm_deriv_array(nMax, m, argument, data_l, ddata_l);

		if(m == 0)
			n = 1;
		else
			n = m;

		for(int i=0; i<=(nMax-m); i++)
		{
			if((m == 0) && (i==0)) //For (0,0), skip first
			{
				i++;
			}

			q.init(n,m);
			data[q] = data_l[i];
			ddata[q] = ddata_l[i];

			//Make negatives
			if(m > 0)
			{
				q.init(n,-m);

				double neg_factor = pow(-1.0, m) *
						(gsl_sf_fact(n-m) /
						 gsl_sf_fact(n+m));

				data[q] = neg_factor * data_l[i];
				ddata[q] = neg_factor * ddata_l[i];
			}

			n++;
		}

		delete [] data_l;
		delete [] ddata_l;
	}

	return 0;
}


