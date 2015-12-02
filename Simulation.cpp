#include "Simulation.h"

#include "Reader.h"
#include "CompoundIterator.h"
#include "Excitation.h"
#include "Solver.h"
#include "Result.h"
#include "Run.h"
#include "aliases.h"
#include "Output.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

int Simulation::run()
{

  //Read the case file
  Run run;

  Reader reader;
  reader.init(&run);
  if(reader.readSimulation(caseFile + ".xml"))
    return 1;

  //Initialize the solver

  Solver solver;

  solver.init(&(run.geometry), &(run.excitation), O3DSolverIndirect, run.nMax);
  solver.populate();

  //Determine the simulation type and proceed accordingly
  if(run.outputType == 0) //Field simulation
  {

    Output oFile;
    oFile.init((caseFile + ".h5").c_str());

    Result result;
    result.init(&(run.geometry), &(run.excitation), run.nMax);
    solver.solve(result.scatter_coef, result.internal_coef);

    OutputGrid oEGrid;
    OutputGrid oHGrid;
    oEGrid.init(O3DCartesianRegular, run.params, oFile.getHandle("Field_E"));
    oHGrid.init(O3DCartesianRegular, run.params, oFile.getHandle("Field_H"));

    if(run.singleMode)
    {
      if(run.dominantAuto)
      {
        CompoundIterator p;
        p = result.getDominant();
        result.setFieldsModal(oEGrid, oHGrid, run.projection, p, run.singleComponent);
        std::cout << "Field output finished. Mode given is for n = " << p.first << " and m = " << p. second << "." << std::endl;
      }
      else
      {
        result.setFieldsModal(oEGrid, oHGrid, run.projection, run.singleModeIndex, run.singleComponent);
      }
    }
    else
    {
      result.setFields(oEGrid, oHGrid, run.projection);
    }

    oEGrid.close();
    oHGrid.close();
    oFile.close();
  }

  if(run.outputType == 11) //Wavelength scan
  {
    Result result;

    std::ofstream outASec(caseFile + "_AbsorptionCS.dat");
    std::ofstream outESec(caseFile + "_ExtinctionCS.dat");

    //Now scan over the wavelengths given in params
    double lami = run.params[0];
    double lamf = run.params[1];
    int steps = run.params[2];

    double lam;
    double lams;

    lams = (lamf - lami) / (steps-1);

    for(int i=0; i<steps; i++)
    {
      lam = lami + i * lams;

      std::cout << "Solving for Lambda = " << lam << std::endl;

      run.excitation.updateWavelength(lam);
      run.geometry.update(&(run.excitation));
      solver.update(&(run.geometry), &(run.excitation), run.nMax);

      Result result;
      result.init(&(run.geometry), &(run.excitation), run.nMax);
      solver.solve(result.scatter_coef, result.internal_coef);

      outASec << lam << "\t" << result.getAbsorptionCrossSection() << std::endl;
      outESec << lam << "\t" << result.getExtinctionCrossSection() << std::endl;
    }

    outASec.close();
    outESec.close();
  }

  if(run.outputType == 12) //Radius scan
  {
    Result result;

    std::ofstream outASec(caseFile + "_AbsorptionCS.dat");
    std::ofstream outESec(caseFile + "_ExtinctionCS.dat");

    //Now scan over the wavelengths given in params
    double radi = run.params[3];
    double radf = run.params[4];
    int radsteps = run.params[5];

    double rad;
    double rads;

    rads = (radf - radi) / (radsteps-1);

    for(int i=0; i<radsteps; i++)
    {
      rad = radi + i * rads;

      std::cout << "Solving for R = " << rad << std::endl;

      for(int k=0; k<run.geometry.noObjects; k++)
      {
        run.geometry.updateRadius(rad, k);
      }

      if(run.geometry.structureType == 1)
      {
        run.geometry.rebuildStructure();
      }

      if(!run.geometry.validate())
      {
        std::cerr << "Geometry no longer valid!";
        exit(1);
      }

      solver.update(&(run.geometry), &(run.excitation), run.nMax);

      Result result;
      result.init(&(run.geometry), &(run.excitation), run.nMax);
      solver.solve(result.scatter_coef, result.internal_coef);

      outASec << rad << "\t" << result.getAbsorptionCrossSection() << std::endl;
      outESec << rad << "\t" << result.getExtinctionCrossSection() << std::endl;
    }

    outASec.close();
    outESec.close();
  }

  if(run.outputType == 112) //R vs Wavelength scan
  {
    Result result;

    std::ofstream outASec(caseFile + "_AbsorptionCS.dat");
    std::ofstream outESec(caseFile + "_ExtinctionCS.dat");
    std::ofstream outParams(caseFile + "_RadiusLambda.dat");

    //Now scan over the wavelengths given in params
    double lami = run.params[0];
    double lamf = run.params[1];
    int lamsteps = run.params[2];

    double radi = run.params[3];
    double radf = run.params[4];
    int radsteps = run.params[5];

    double lam, lams, rad, rads;

    lams = (lamf - lami) / (lamsteps-1);
    rads = (radf - radi) / (radsteps-1);

    for(int i=0; i<lamsteps; i++)
    {
      lam = lami + i * lams;

      for(int j=0; j<radsteps; j++)
      {
        rad = radi + j * rads;

        std::cout << "Solving for Lambda = " << lam << " and R =" << rad << std::endl;

        run.excitation.updateWavelength(lam);
        run.geometry.update(&(run.excitation));
        for(int k=0; k<run.geometry.noObjects; k++)
        {
          run.geometry.updateRadius(rad, k);
        }

        if(run.geometry.structureType == 1)
        {
          run.geometry.rebuildStructure();
        }

        if(!run.geometry.validate())
        {
          std::cerr << "Geometry no longer valid!";
          exit(1);
        }

        solver.update(&(run.geometry), &(run.excitation), run.nMax);

        Result result;
        result.init(&(run.geometry), &(run.excitation), run.nMax);
        solver.solve(result.scatter_coef, result.internal_coef);

        outASec << result.getAbsorptionCrossSection() << "\t";
        outESec << result.getExtinctionCrossSection() << "\t";
        outParams << "(" << rad*1e9 << " , " << lam*1e9 << ")" << "\t";
      }

      outASec << std::endl;
      outESec << std::endl;
      outParams << std::endl;
    }

    outASec.close();
    outESec.close();
    outParams.close();
  }

  if(run.outputType == 2)
  {
    //Scattering coefficients requests
    std::ofstream outPCoef(caseFile + "_pCoefficients.dat");
    std::ofstream outQCoef(caseFile + "_qCoefficients.dat");

    Result result;
    result.init(&(run.geometry), &(run.excitation), run.nMax);
    solver.solve(result.scatter_coef, result.internal_coef);

    CompoundIterator p;

    for(p=0; p<p.max(run.nMax); p++)
    {
      outPCoef << p.first << "\t" << p.second << "\t" << abs(result.scatter_coef[p]) << std::endl;
      outQCoef << p.first << "\t" << p.second << "\t" << abs(result.scatter_coef[p.compound + p.max(run.nMax)]) << std::endl;
    }

    outPCoef.close();
    outQCoef.close();
  }

  return 0;
}

int Simulation::done()
{
  //Placeholder method. Not needed at the moment.
  return 0;
}
