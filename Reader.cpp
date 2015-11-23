#include "Reader.h"

#include "Geometry.h"
#include "Scatterer.h"
#include "Spherical.h"
#include "Cartesian.h"
#include "constants.h"
#include "Tools.h"
#include <iostream>
#include <cstring>

using std::cerr;
using std::endl;
using std::strcmp;

using namespace pugi;

Reader::Reader()
{
	initDone = false;
}

Reader::~Reader()
{
	//
}

int Reader::readGeometry()
{
	xml_node geo_node;	 		//the main work node
	unsigned long no_objects=0;	//the number of declared objects

	//Find the simulation node
	geo_node = inputFile.child("simulation");
	if (!geo_node)
	{
		cerr << "Simulation parameters not defined!" << endl;
		return 1;
	}

	run->nMax = geo_node.child("harmonics").attribute("nmax").as_int();

	//Find the geometry node
	geo_node = inputFile.child("geometry");
	if (!geo_node)
	{
		cerr << "Geometry not defined!" << endl;
		return 1;
	}

	//Check if a structure is defined
	if(geo_node.child("structure"))
	{
		return readStructure(geo_node);
	}

	run->geometry.structureType = 0;

	//Find all objects and count them
	for(xml_node node = geo_node.child("object"); node; node = node.next_sibling("object"))
	{
		no_objects++;
	}

	run->geometry.init(no_objects);	//initialize the geometry

	//Find all objects again and add them NOTE: the parsing is really fast so this is not a problem
	for(xml_node node = geo_node.child("object"); node; node = node.next_sibling("object"))
	{
		Scatterer work_object; 		//a work object used to push into the Geometry (not initialized)
		work_object.nMax = run->nMax;

		//Choose object type (since comparing strings, elif is better)

		if (!strcmp(node.attribute("type").value(), "sphere")) //Sphere object
		{
			//Assign coordinates to the Scatterer work_object
			if(node.child("cartesian")) //Cartesian coordinates
			{
				double aux_x, aux_y, aux_z;

				aux_x = node.child("cartesian").attribute("x").as_double() * consFrnmTom;
				aux_y = node.child("cartesian").attribute("y").as_double() * consFrnmTom;
				aux_z = node.child("cartesian").attribute("z").as_double() * consFrnmTom;

				work_object.vR = Tools::toSpherical(Cartesian<double> (aux_x, aux_y, aux_z));
			}
			else if(node.child("spherical")) //Spherical coordinates
			{
				double aux_r, aux_t, aux_p;

				aux_r = node.child("spherical").attribute("rrr").as_double() * consFrnmTom;
				aux_t = node.child("spherical").attribute("the").as_double();
				aux_p = node.child("spherical").attribute("phi").as_double();

				work_object.vR = Spherical<double>(aux_r, aux_t, aux_p);
			}

			//Assign properties to the Scatterer work_object
			if(node.child("properties").attribute("radius"))
			{
				work_object.radius = node.child("properties").attribute("radius").as_double() * consFrnmTom;
			}

			//Assign electromagnetic properties to the Scatterer
			if(node.child("epsilon") || node.child("mu"))
			{
				complex<double> aux_epsilon = complex<double>(1.0, 0.);
				complex<double> aux_mu = complex<double>(1.0, 0.);
				// Drude model
				complex<double> plasma_freq = complex<double>(1.0, 0.);
				complex<double> damping_freq = complex<double>(1.0, 0.);
//				double input_frequency = 1.0; //consC/wavelength;


				//Read mu first
				if(!strcmp(node.child("mu").attribute("type").value(), "relative"))
				{
					aux_mu.real(node.child("mu").attribute("value.real").as_double());
					aux_mu.imag(node.child("mu").attribute("value.imag").as_double());
				}

				//Now the two epsilon models

				if(!strcmp(node.child("epsilon").attribute("type").value(), "relative"))
				{
					// Static values
					aux_epsilon.real(node.child("epsilon").attribute("value.real").as_double());
					aux_epsilon.imag(node.child("epsilon").attribute("value.imag").as_double());
					work_object.elmag.init_r(aux_epsilon, aux_mu);
				}

				if(!strcmp(node.child("epsilon").attribute("type").value(), "DrudeModel"))
				{
					//Drude model
//					input_frequency = node.child("epsilon").child("parameters").attribute("input_frequency").as_double();
					plasma_freq.real(node.child("epsilon").child("parameters").attribute("plasma_frequency").as_double());
					damping_freq.real(node.child("epsilon").child("parameters").attribute("damping_frequency").as_double());
					damping_freq*=complex<double>(0., 1.);
	//				aux_epsilon = 1. - ( (plasma_freq*plasma_freq) / (input_freq*(input_freq+damping_freq)) );
		//			work_object.elmag.init_r(aux_epsilon, aux_mu);
					work_object.elmag.initDrudeModel_r(plasma_freq, damping_freq, aux_mu);

				}

				if(!strcmp(node.child("epsilon").attribute("type").value(), "sellmeier"))
				{
					//Sellmeier model
					double B1(0.), C1(0.), B2(0.), C2(0.), B3(0.), C3(0.), B4(0.), C4(0.), B5(0.), C5(0.);
					B1=node.child("epsilon").child("parameters").attribute("B1").as_double();
					C1=node.child("epsilon").child("parameters").attribute("C1").as_double();
					B2=node.child("epsilon").child("parameters").attribute("B2").as_double();
					C2=node.child("epsilon").child("parameters").attribute("C2").as_double();
					B3=node.child("epsilon").child("parameters").attribute("B3").as_double();
					C3=node.child("epsilon").child("parameters").attribute("C3").as_double();
					B4=node.child("epsilon").child("parameters").attribute("B4").as_double();
					C4=node.child("epsilon").child("parameters").attribute("C4").as_double();
					B5=node.child("epsilon").child("parameters").attribute("B5").as_double();
					C5=node.child("epsilon").child("parameters").attribute("C5").as_double();
					work_object.elmag.initSellmeier_r(	B1,
														C1,
														B2,
														C2,
														B3,
														C3,
														B4,
														C4,
														B5,
														C5,
														aux_mu.real());
				}
		}

		//Push the object into the geometry
		run->geometry.pushObject(work_object);
		}
	}

	//Add the background properties
	if(geo_node.child("background"))
    {
//        if(strcmp(geo_node.child("background").attribute("type").value(), "absolute"))
  //      {
    //        run->geometry.bground.init(geo_node.child("background").child("epsilon").attribute("value").as_double(), geo_node.child("background").child("mu").attribute("value").as_double());
      // }
		if(strcmp(geo_node.child("background").attribute("type").value(), "relative"))
        {
			complex<double> aux_epsilon = complex<double>(1.0, 0.);
			complex<double> aux_mu = complex<double>(1.0, 0.);
            aux_epsilon.real(geo_node.child("background").child("epsilon").attribute("value.real").as_double());
			aux_epsilon.imag(geo_node.child("background").child("epsilon").attribute("value.imag").as_double());
			aux_mu.real(geo_node.child("background").child("mu").attribute("value.real").as_double());
			aux_mu.imag(geo_node.child("background").child("mu").attribute("value.imag").as_double());
			run->geometry.bground.init_r(aux_epsilon, aux_mu);
       }
//        else if(strcmp(geo_node.child("background").attribute("type").value(), "relative"))
  //      {
	//	run->geometry.bground.init_r(geo_node.child("background").child("epsilon").attribute("value").as_double(), geo_node.child("background").child("mu").attribute("value").as_double());
		//}
    }

	//Validate the geometry in the return

	return run->geometry.validate();
}

int Reader::readExcitation()
{
    xml_node ext_node;  //the main work node

    //Find the source node
    ext_node = inputFile.child("source");
    if(!ext_node)
    {
        cerr << "Source not defined!" << endl;
        return 1;
    }

    int source_type;
    double wavelength;
    SphericalP<complex<double> > Einc(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
    Spherical<double> vKinc(0.0, 0.0, 0.0);

    //Determine source type
    if(!strcmp(ext_node.attribute("type").value(), "planewave"))
        source_type = 0;
    else //Default is always planewave
        source_type = 0;

    //Determine wavelength
    wavelength = ext_node.child("wavelength").attribute("value").as_double();
    wavelength *= 1e-9;

    //Determine propagation values
    vKinc = Spherical<double>(2*consPi/wavelength, ext_node.child("propagation").attribute("theta").as_double() * consPi/180.0, ext_node.child("propagation").attribute("phi").as_double() * consPi/180.0);

    //Determine polarisation (initial field values)
    SphericalP<complex<double> > Eaux;
    Spherical<double> vAux = Spherical<double>(0.0, vKinc.the, vKinc.phi);
    Eaux = SphericalP<complex<double> >(complex<double>(0.0, 0.0),
    									complex<double>(ext_node.child("polarization").attribute("Etheta.real").as_double(),
    											ext_node.child("polarization").attribute("Etheta.imag").as_double()),
    									complex<double>(ext_node.child("polarization").attribute("Ephi.real").as_double(),
    											ext_node.child("polarization").attribute("Ephi.imag").as_double()));
    Einc = Tools::toProjection(vAux, Eaux);

    //Initialize and populate the excitation
    run->excitation.init(source_type, Einc, vKinc, run->nMax);
    run->excitation.populate();

    //Update the geometry in case we had dynamic models
    run->geometry.update(&(run->excitation));

    return 0;
}

int Reader::readStructure(xml_node geo_node_)
{
	xml_node struct_node = geo_node_.child("structure");

	run->geometry.structureType = 1; //set the spiral structure flag

	if(!strcmp(struct_node.attribute("type").value(), "spiral"))
	{
		//Build a spiral
		double R, d; //radius (length) and distance between spheres
		int Np, No; //number of points/objects
		double Theta;

		Np = struct_node.child("properties").attribute("points").as_int();
		Theta = consPi/(Np-1); //Calculate the separation angle

		int arms = struct_node.attribute("arms").as_int();

		if(!strcmp(struct_node.child("properties").attribute("length").value(), ""))
		{
			//Distance based simulation;
			d = 2.0 * struct_node.child("properties").attribute("distance").as_double() * consFrnmTom;
			R = d / (4 * sin(Theta/2.0)); //Remember that the radius needs to be half the length
		}
		else
		{
			//Length based simulation
			R = struct_node.child("properties").attribute("length").as_double() * consFrnmTom;
			d = 2 * R * sin(Theta/2.0); //Not really used but may be useful
			R = R/2.0; //Remember that the radius needs to be half the length
		}

		No = (Np-1)*arms + 1; //Number of objects

		run->geometry.init(No);

		//Instantiate a work object and push the central sphere
		Scatterer work_object; 		//a work object used to push into the Geometry (not initialized)
		work_object.nMax = run->nMax;

		//Assign properties to the Scatterer work_object
		if(struct_node.child("object").child("properties").attribute("radius"))
		{
			work_object.radius = struct_node.child("object").child("properties").attribute("radius").as_double() * consFrnmTom;
			run->geometry.spiralSeparation = (d/2) - 2 * work_object.radius;
		}

		//Assign electromagnetic properties to the Scatterer
		if(struct_node.child("object").child("epsilon") || struct_node.child("object").child("mu"))
		{
			double aux_epsilon = 1.0;
			double aux_mu = 1.0;

			//Read mu first
			if(!strcmp(struct_node.child("object").child("mu").attribute("type").value(), "relative"))
			{
				aux_mu = struct_node.child("object").child("mu").attribute("value").as_double();
			}

			//Now the two epsilon models

			if(!strcmp(struct_node.child("object").child("epsilon").attribute("type").value(), "relative"))
			{
				aux_epsilon = struct_node.child("object").child("epsilon").attribute("value").as_double();
				work_object.elmag.init_r(complex<double>(aux_epsilon, 0.0), complex<double>(aux_mu, 0.0));
			}

			if(!strcmp(struct_node.child("object").child("epsilon").attribute("type").value(), "sellmeier"))
			{
				//Sellmeier model
					double B1(0.), C1(0.), B2(0.), C2(0.), B3(0.), C3(0.), B4(0.), C4(0.), B5(0.), C5(0.);
					B1=struct_node.child("epsilon").child("parameters").attribute("B1").as_double();
					C1=struct_node.child("epsilon").child("parameters").attribute("C1").as_double();
					B2=struct_node.child("epsilon").child("parameters").attribute("B2").as_double();
					C2=struct_node.child("epsilon").child("parameters").attribute("C2").as_double();
					B3=struct_node.child("epsilon").child("parameters").attribute("B3").as_double();
					C3=struct_node.child("epsilon").child("parameters").attribute("C3").as_double();
					B4=struct_node.child("epsilon").child("parameters").attribute("B4").as_double();
					C4=struct_node.child("epsilon").child("parameters").attribute("C4").as_double();
					B5=struct_node.child("epsilon").child("parameters").attribute("B5").as_double();
					C5=struct_node.child("epsilon").child("parameters").attribute("C5").as_double();
					work_object.elmag.initSellmeier_r(	B1,
														C1,
														B2,
														C2,
														B3,
														C3,
														B4,
														C4,
														B5,
														C5,
														aux_mu);

				//Sellmeier model
//				work_object.elmag.initSellmeier_r(struct_node.child("object").child("epsilon").child("parameters").attribute("B1").as_double(),
	//					struct_node.child("object").child("epsilon").child("parameters").attribute("C1").as_double(),
		//				struct_node.child("object").child("epsilon").child("parameters").attribute("B2").as_double(),
			//			struct_node.child("object").child("epsilon").child("parameters").attribute("C2").as_double(),
			//			struct_node.child("object").child("epsilon").child("parameters").attribute("B3").as_double(),
			//			struct_node.child("object").child("epsilon").child("parameters").attribute("C3").as_double(),
			//			aux_mu);
			}
		}

		//Center the sphere and push it

		work_object.vR.rrr = 0.0;
		work_object.vR.the = 0.0;
		work_object.vR.phi = 0.0;

		run->geometry.pushObject(work_object);

		//Create vectors for r, theta, x and y, X and Y
		double X[No-1];			// to store x-axis location of particles for all arm
		double Y[No-1];			// to store y-axis location of particles for all arm

		int i=0, j=0;

		// AJ -----------------------------------------------------------------------
		double Theta_rot = 2*consPi/arms;
		double x[Np-1];			// to store x-axis location of particles for each arm
		double y[Np-1];			// to store y-axis location of particles for each arm
		j=0;
		for(i=0; i<Np-1; i++){
			double x_;
			double y_;
			Tools::Pol2Cart(R, i*Theta, x_, y_);
				x[j] = x_ + R;
				y[j] = y_;
				j++;
		}

		j=0;
		for(int arm_i=0; arm_i<arms; arm_i++){
			// Translate points around origin to finally define locations of all particles on arms
			for(i=0; i<Np-1; i++){
				X[j] = x[i]*cos(double(arm_i)*Theta_rot) - y[i]*sin(double(arm_i)*Theta_rot);
				Y[j] = x[i]*sin(double(arm_i)*Theta_rot) + y[i]*cos(double(arm_i)*Theta_rot);
				j++;
			}
		}

		//Determine normal, convert to a spherical object and push

		Cartesian<double> auxCar(0.0, 0.0, 0.0);
		Spherical<double> auxSph(0.0, 0.0, 0.0);

		for(int i=0; i<No-1; i++)
		{
			if(!strcmp(struct_node.child("properties").attribute("normal").value(), "x"))
			{
				//x is normal (conversion is x(pol) -> y; y(pol) -> z
				auxCar.x = 0.0;
				auxCar.y = X[i];
				auxCar.z = Y[i];

				run->geometry.normalToSpiral = 0;
			}

			if(!strcmp(struct_node.child("properties").attribute("normal").value(), "y"))
			{
				//y is normal (conversion is x(pol) -> z; y(pol) -> x
				auxCar.x = Y[i];
				auxCar.y = 0.0;
				auxCar.z = X[i];

				run->geometry.normalToSpiral = 1;
			}

			if(!strcmp(struct_node.child("properties").attribute("normal").value(), "z"))
			{
				//z is normal (conversion is x(pol) -> x; y(pol) -> x
				auxCar.x = X[i];
				auxCar.y = Y[i];
				auxCar.z = 0.0;

				run->geometry.normalToSpiral = 2;
			}

			auxSph = Tools::toSpherical(auxCar);
			work_object.vR = auxSph;
			run->geometry.pushObject(work_object);
		}
	}

	return run->geometry.validate();
}

int Reader::readOutput()
{
    xml_node out_node;  //the main work node

    //Find the source node
    out_node = inputFile.child("output");
    if(!out_node)
    {
        cerr << "Output not defined!" << endl;
        return 1;
    }

    //Determine type
    if(!strcmp(out_node.attribute("type").value(), "coefficients"))
    {
    	run->outputType = 2;
    }

    if(!strcmp(out_node.attribute("type").value(), "field"))
    {
    	run->outputType = 0;	//Field output requested
    	run->singleComponent = 0; //Set this to zero as default

    	run->params[0] = out_node.child("grid").child("x").attribute("min").as_double() * 1e-9;
    	run->params[1] = out_node.child("grid").child("x").attribute("max").as_double()* 1e-9;
    	run->params[2] = out_node.child("grid").child("x").attribute("steps").as_double();
    	run->params[3] = out_node.child("grid").child("y").attribute("min").as_double()* 1e-9;
    	run->params[4] = out_node.child("grid").child("y").attribute("max").as_double()* 1e-9;
    	run->params[5] = out_node.child("grid").child("y").attribute("steps").as_double();
    	run->params[6] = out_node.child("grid").child("z").attribute("min").as_double()* 1e-9;
    	run->params[7] = out_node.child("grid").child("z").attribute("max").as_double()* 1e-9;
    	run->params[8] = out_node.child("grid").child("z").attribute("steps").as_double();

    	if(!strcmp(out_node.child("projection").attribute("spherical").value(), "true"))
    	{
    		run->projection = 1;
    	}
    	else
    	{
    		run->projection = 0;
    	}

    	if(out_node.child("singlemode"))
    	{
    		run->singleMode = true;
    	}
    	else
    	{
    		run->singleMode = false;
    	}

    	if(!strcmp(out_node.child("singlemode").attribute("dominant").value(), "auto"))
    	{
    		run->dominantAuto = true;
    	}
    	else
    	{
    		run->dominantAuto = false;
    		run->singleModeIndex.init(out_node.child("singlemode").attribute("n").as_int(),
    				out_node.child("singlemode").attribute("m").as_int());
    		if(!strcmp(out_node.child("singlemode").attribute("component").value(), "TE"))
    		{
    			run->singleComponent = 1;
    		}
    		if(!strcmp(out_node.child("singlemode").attribute("component").value(), "TM"))
    		{
    			run->singleComponent = 2;
    		}
    	}
    }

    if(!strcmp(out_node.attribute("type").value(), "response"))
    {
    	if(out_node.child("scan").child("wavelength"))
    	{
        	double lam_start(0.), lam_final(0.);
			lam_start = out_node.child("scan").child("wavelength").attribute("initial").as_double();
            lam_final = out_node.child("scan").child("wavelength").attribute("final").as_double();
			run->params[0] = out_node.child("scan").child("wavelength").attribute("initial").as_double() * 1e-9;
            run->params[1] = out_node.child("scan").child("wavelength").attribute("final").as_double()* 1e-9;

			int stepsize(0), steps(0);
			stepsize = out_node.child("scan").child("wavelength").attribute("stepsize").as_double();

			// claculate no of steps
			steps = int(lam_final - lam_start) / stepsize;
            run->params[2] = steps+1;
//			run->params[2] = out_node.child("scan").child("wavelength").attribute("steps").as_double();
            run->outputType = 11;
    	}

    	if(out_node.child("scan").child("radius"))
    	{
        	run->params[3] = out_node.child("scan").child("radius").attribute("initial").as_double() * 1e-9;
            run->params[4] = out_node.child("scan").child("radius").attribute("final").as_double()* 1e-9;
            run->params[5] = out_node.child("scan").child("radius").attribute("steps").as_double();
            run->outputType = 12;
    	}

    	if(out_node.child("scan").child("wavelength") && out_node.child("scan").child("radius"))
    		run->outputType = 112;
    }

   return 0;
}

Reader::Reader(Run *run_)
{
	init(run_);
}

void Reader::init(Run *run_)
{
	run = run_;
	initDone = true;
}

int Reader::readSimulation(char* fileName_)
{
	xml_parse_result fileResult;

	fileResult = inputFile.load_file(fileName_);

	if(!fileResult)
	{
		cerr << "Error reading or parsing input file " << fileName_ << "!" << endl;
		return 1;
	}

	//Read the Geometry
	if(!readGeometry())
	{
		cerr << "Geometry not valid!";
		return 1;
	}

    //Read Excitation
    if(readExcitation())
    {
        cerr << "Source not valid!";
        return 1;
    }

    //Read Excitation
    if(readOutput())
    {
        cerr << "Output not valid!";
        return 1;
    }


	return 0;
}
