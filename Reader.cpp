#include "Reader.h"

#include "Cartesian.h"
#include "Geometry.h"
#include "Scatterer.h"
#include "Spherical.h"
#include "Tools.h"
#include "constants.h"
#include "mpi/Communicator.h"
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

#ifdef OPTIMET_BELOS
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#endif

using namespace pugi;

Reader::Reader() { initDone = false; }

Reader::~Reader() {
  //
}

int Reader::readGeometry() {
  xml_node geo_node; // the main work node

  // Find the simulation node
  geo_node = inputFile.child("simulation");
  if(!geo_node) {
    std::cerr << "Simulation parameters not defined!" << std::endl;
    return 1;
  }

  run->nMax = geo_node.child("harmonics").attribute("nmax").as_int();

  // Find the geometry node
  geo_node = inputFile.child("geometry");
  if(!geo_node) {
    std::cerr << "Geometry not defined!" << std::endl;
    return 1;
  }

  // Check if a structure is defined
  if(geo_node.child("structure")) {
    return readStructure(geo_node);
  }

  run->geometry->structureType = 0;

  // Find all scattering objects
  for(xml_node node = geo_node.child("object"); node; node = node.next_sibling("object"))
    run->geometry->pushObject(readSphericalScatterer(node));

  // Add the background properties
  if(geo_node.child("background")) {
    if(std::strcmp(geo_node.child("background").attribute("type").value(), "relative")) {
      std::complex<double> aux_epsilon(1.0, 0.);
      std::complex<double> aux_mu(1.0, 0.);
      aux_epsilon.real(
          geo_node.child("background").child("epsilon").attribute("value.real").as_double());
      aux_epsilon.imag(
          geo_node.child("background").child("epsilon").attribute("value.imag").as_double());
      aux_mu.real(geo_node.child("background").child("mu").attribute("value.real").as_double());
      aux_mu.imag(geo_node.child("background").child("mu").attribute("value.imag").as_double());
      run->geometry->bground.init_r(aux_epsilon, aux_mu);
    }
  }

  // Validate the geometry in the return
  if(run->geometry->objects.size() == 0)
    throw std::runtime_error("No scatterers defined in input");
  return 1;
}

int Reader::readExcitation() {
  xml_node ext_node; // the main work node

  // Find the source node
  ext_node = inputFile.child("source");
  if(!ext_node) {
    std::cerr << "Source not defined!" << std::endl;
    return 1;
  }

  int source_type;
  double wavelength;
  SphericalP<std::complex<double>> Einc(std::complex<double>(0.0, 0.0),
                                        std::complex<double>(0.0, 0.0),
                                        std::complex<double>(0.0, 0.0));
  Spherical<double> vKinc(0.0, 0.0, 0.0);

  // Determine source type
  if(!std::strcmp(ext_node.attribute("type").value(), "planewave"))
    source_type = 0;
  else // Default is always planewave
    source_type = 0;

  // Determine wavelength
  wavelength = ext_node.child("wavelength").attribute("value").as_double();
  wavelength *= 1e-9;

  // Determine propagation values
  vKinc = Spherical<double>(
      2 * consPi / wavelength,
      ext_node.child("propagation").attribute("theta").as_double() * consPi / 180.0,
      ext_node.child("propagation").attribute("phi").as_double() * consPi / 180.0);

  // Determine polarisation (initial field values)
  SphericalP<std::complex<double>> Eaux;
  Spherical<double> vAux = Spherical<double>(0.0, vKinc.the, vKinc.phi);
  Eaux = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0),
      std::complex<double>(ext_node.child("polarization").attribute("Etheta.real").as_double(),
                           ext_node.child("polarization").attribute("Etheta.imag").as_double()),
      std::complex<double>(ext_node.child("polarization").attribute("Ephi.real").as_double(),
                           ext_node.child("polarization").attribute("Ephi.imag").as_double()));
  Einc = Tools::toProjection(vAux, Eaux);

  // Initialize and populate the excitation
  run->excitation = std::make_shared<optimet::Excitation>(source_type, Einc, vKinc, run->nMax);
  run->excitation->populate();

  // Update the geometry in case we had dynamic models
  run->geometry->update(run->excitation);

  return 0;
}

int Reader::readStructure(xml_node geo_node_) {
  xml_node struct_node = geo_node_.child("structure");

  run->geometry->structureType = 1; // set the spiral structure flag

  if(!std::strcmp(struct_node.attribute("type").value(), "spiral")) {
    // Build a spiral
    double R, d; // radius (length) and distance between spheres
    int Np, No;  // number of points/objects
    double Theta;

    Np = struct_node.child("properties").attribute("points").as_int();
    Theta = consPi / (Np - 1); // Calculate the separation angle

    int arms = struct_node.attribute("arms").as_int();

    if(!std::strcmp(struct_node.child("properties").attribute("length").value(), "")) {
      // Distance based simulation;
      d = 2.0 * struct_node.child("properties").attribute("distance").as_double() * consFrnmTom;
      R = d / (4 * sin(Theta / 2.0)); // Remember that the radius needs to be half the length
    } else {
      // Length based simulation
      R = struct_node.child("properties").attribute("length").as_double() * consFrnmTom;
      d = 2 * R * sin(Theta / 2.0); // Not really used but may be useful
      R = R / 2.0;                  // Remember that the radius needs to be half the length
    }

    No = (Np - 1) * arms + 1; // Number of objects

    // Assign properties to the Scatterer work_object
    if(struct_node.child("object").child("properties").attribute("radius")) {
      auto const radius =
          struct_node.child("object").child("properties").attribute("radius").as_double();
      run->geometry->spiralSeparation = (d / 2) - 2 * radius * consFrnmTom;
    }

    // Create vectors for r, theta, x and y, X and Y
    std::vector<double> X(No - 1); // to store x-axis location of particles for all arm
    std::vector<double> Y(No - 1); // to store y-axis location of particles for all arm

    int i = 0, j = 0;

    // AJ
    // -----------------------------------------------------------------------
    double Theta_rot = 2 * consPi / arms;
    std::vector<double> x(Np - 1); // to store x-axis location of particles for each arm
    std::vector<double> y(Np - 1); // to store y-axis location of particles for each arm
    j = 0;
    for(i = 0; i < Np - 1; i++) {
      double x_;
      double y_;
      Tools::Pol2Cart(R, i * Theta, x_, y_);
      x[j] = x_ + R;
      y[j] = y_;
      j++;
    }

    j = 0;
    for(int arm_i = 0; arm_i < arms; arm_i++) {
      // Translate points around origin to finally define locations of all
      // particles on arms
      for(i = 0; i < Np - 1; i++) {
        X[j] = x[i] * cos(double(arm_i) * Theta_rot) - y[i] * sin(double(arm_i) * Theta_rot);
        Y[j] = x[i] * sin(double(arm_i) * Theta_rot) + y[i] * cos(double(arm_i) * Theta_rot);
        j++;
      }
    }

    // Determine normal, convert to a spherical object and push

    auto const scatterer = readSphericalScatterer(struct_node.child("object"));
    for(int i = 0; i < No - 1; i++) {
      run->geometry->pushObject(scatterer);
      std::string const normal = struct_node.child("properties").attribute("normal").value();
      if(normal == "x") {
        // x is normal (conversion is x(pol) -> y; y(pol) -> z
        run->geometry->normalToSpiral = 0;
        run->geometry->objects.back().vR = Tools::toSpherical({0.0, X[i], Y[i]});
      } else if(normal == "y") {
        // y is normal (conversion is x(pol) -> z; y(pol) -> x
        run->geometry->normalToSpiral = 1;
        run->geometry->objects.back().vR = Tools::toSpherical({Y[i], 0, X[i]});
      } else if(normal == "z") {
        // z is normal (conversion is x(pol) -> x; y(pol) -> x
        run->geometry->normalToSpiral = 2;
        run->geometry->objects.back().vR = Tools::toSpherical({X[i], Y[i], 0});
      } else
        throw std::runtime_error("Unknown normal " + normal);
    }
  }

  if(run->geometry->objects.size() == 0)
    throw std::runtime_error("No scatterers defined in input");
  return 1;
}

int Reader::readOutput() {
  xml_node out_node; // the main work node

  // Find the source node
  out_node = inputFile.child("output");
  if(!out_node) {
    std::cerr << "Output not defined!" << std::endl;
    return 1;
  }

  // Determine type
  if(!std::strcmp(out_node.attribute("type").value(), "coefficients")) {
    run->outputType = 2;
  }

  if(!std::strcmp(out_node.attribute("type").value(), "field")) {
    run->outputType = 0;      // Field output requested
    run->singleComponent = 0; // Set this to zero as default

    run->params[0] = out_node.child("grid").child("x").attribute("min").as_double() * 1e-9;
    run->params[1] = out_node.child("grid").child("x").attribute("max").as_double() * 1e-9;
    run->params[2] = out_node.child("grid").child("x").attribute("steps").as_double();
    run->params[3] = out_node.child("grid").child("y").attribute("min").as_double() * 1e-9;
    run->params[4] = out_node.child("grid").child("y").attribute("max").as_double() * 1e-9;
    run->params[5] = out_node.child("grid").child("y").attribute("steps").as_double();
    run->params[6] = out_node.child("grid").child("z").attribute("min").as_double() * 1e-9;
    run->params[7] = out_node.child("grid").child("z").attribute("max").as_double() * 1e-9;
    run->params[8] = out_node.child("grid").child("z").attribute("steps").as_double();

    if(!std::strcmp(out_node.child("projection").attribute("spherical").value(), "true")) {
      run->projection = 1;
    } else {
      run->projection = 0;
    }

    if(out_node.child("singlemode")) {
      run->singleMode = true;
    } else {
      run->singleMode = false;
    }

    if(!std::strcmp(out_node.child("singlemode").attribute("dominant").value(), "auto")) {
      run->dominantAuto = true;
    } else {
      run->dominantAuto = false;
      run->singleModeIndex.init(out_node.child("singlemode").attribute("n").as_int(),
                                out_node.child("singlemode").attribute("m").as_int());
      if(!std::strcmp(out_node.child("singlemode").attribute("component").value(), "TE")) {
        run->singleComponent = 1;
      }
      if(!std::strcmp(out_node.child("singlemode").attribute("component").value(), "TM")) {
        run->singleComponent = 2;
      }
    }
  }

  if(!std::strcmp(out_node.attribute("type").value(), "response")) {
    if(out_node.child("scan").child("wavelength")) {
      double lam_start(0.), lam_final(0.);
      lam_start = out_node.child("scan").child("wavelength").attribute("initial").as_double();
      lam_final = out_node.child("scan").child("wavelength").attribute("final").as_double();
      run->params[0] =
          out_node.child("scan").child("wavelength").attribute("initial").as_double() * 1e-9;
      run->params[1] =
          out_node.child("scan").child("wavelength").attribute("final").as_double() * 1e-9;

      int stepsize(0), steps(0);
      stepsize = out_node.child("scan").child("wavelength").attribute("stepsize").as_double();

      // claculate no of steps
      steps = int(lam_final - lam_start) / stepsize;
      run->params[2] = steps + 1;
      //      run->params[2] =
      //      out_node.child("scan").child("wavelength").attribute("steps").as_double();
      run->outputType = 11;
    }

    if(out_node.child("scan").child("radius")) {
      run->params[3] =
          out_node.child("scan").child("radius").attribute("initial").as_double() * 1e-9;
      run->params[4] = out_node.child("scan").child("radius").attribute("final").as_double() * 1e-9;
      run->params[5] = out_node.child("scan").child("radius").attribute("steps").as_double();
      run->outputType = 12;
    }

    if(out_node.child("scan").child("wavelength") && out_node.child("scan").child("radius"))
      run->outputType = 112;
  }

  return 0;
}

Reader::Reader(Run *run_) { init(run_); }

void Reader::init(Run *run_) {
  run = run_;
  initDone = true;
}

int Reader::readSimulation(std::string const &fileName_) {
  xml_parse_result fileResult;

  fileResult = inputFile.load_file(fileName_.c_str());

  if(!fileResult) {
    std::cerr << "Error reading or parsing input file " << fileName_ << "!" << std::endl;
    return 1;
  }

  // Read the Geometry
  if(!readGeometry()) {
    std::cerr << "Geometry not valid!";
    return 1;
  }

  // Read Excitation
  if(readExcitation()) {
    std::cerr << "Source not valid!";
    return 1;
  }

  // Read Excitation
  if(readOutput()) {
    std::cerr << "Output not valid!";
    return 1;
  }

  readParallel(inputFile.child("parallel"), run->parallel_params);
  readParameterList(inputFile);

  return 0;
}

void Reader::readParallel(const pugi::xml_node &node,
                          optimet::scalapack::Parameters &parallel_params) {
  parallel_params.block_size = node.attribute("block_size").as_uint(parallel_params.block_size);
  parallel_params.grid.rows =
      node.child("grid").attribute("rows").as_uint(parallel_params.grid.rows);
  parallel_params.grid.cols =
      node.child("grid").attribute("cols").as_uint(parallel_params.grid.cols);
}

Scatterer Reader::readSphericalScatterer(pugi::xml_node const &node) {
  if(node.attribute("type").value() != std::string("sphere"))
    std::runtime_error("Expecting a spherical scatterer");
  Scatterer result(run->nMax);
  // Assign coordinates to the Scatterer work_object
  if(node.child("cartesian")) // Cartesian coordinates
    result.vR = Tools::toSpherical(
        Cartesian<double>{node.child("cartesian").attribute("x").as_double() * consFrnmTom,
                          node.child("cartesian").attribute("y").as_double() * consFrnmTom,
                          node.child("cartesian").attribute("z").as_double() * consFrnmTom});
  else if(node.child("spherical")) // Spherical coordinates
    result.vR = {node.child("spherical").attribute("rrr").as_double() * consFrnmTom,
                 node.child("spherical").attribute("the").as_double(),
                 node.child("spherical").attribute("phi").as_double()};
  else
    result.vR = {0, 0, 0};

  // Assign properties to the Scatterer work_object
  if(node.child("properties").attribute("radius"))
    result.radius = node.child("properties").attribute("radius").as_double() * consFrnmTom;

  // Assign electromagnetic properties to the Scatterer
  if(node.child("epsilon") || node.child("mu")) {
    // only one way to specify mu, AFAIK
    if(node.child("mu").attribute("type").value() != std::string("relative"))
      throw std::runtime_error("The type for mu must be \"relative\"");

    std::complex<double> const aux_mu(node.child("mu").attribute("value.real").as_double(),
                                      node.child("mu").attribute("value.imag").as_double());

    // Now the two epsilon models
    if(node.child("epsilon").attribute("type").value() == std::string("relative")) {
      // Static values
      std::complex<double> epsilon(node.child("epsilon").attribute("value.real").as_double(),
                                   node.child("epsilon").attribute("value.imag").as_double());
      result.elmag.init_r(epsilon, aux_mu);
    } else if(node.child("epsilon").attribute("type").value() == std::string("DrudeModel")) {
      // Drude model
      auto const plasma_freq =
          node.child("epsilon").child("parameters").attribute("plasma_frequency").as_double();
      std::complex<double> const damping_freq(
          0, node.child("epsilon").child("parameters").attribute("damping_frequency").as_double());
      result.elmag.init_r(0, aux_mu);
      result.elmag.initDrudeModel_r(plasma_freq, damping_freq, aux_mu);
    } else if(node.child("epsilon").attribute("type").value() == std::string("sellmeier")) {
      // Sellmeier model
      double B1(0.), C1(0.), B2(0.), C2(0.), B3(0.), C3(0.), B4(0.), C4(0.), B5(0.), C5(0.);
      B1 = node.child("epsilon").child("parameters").attribute("B1").as_double();
      C1 = node.child("epsilon").child("parameters").attribute("C1").as_double();
      B2 = node.child("epsilon").child("parameters").attribute("B2").as_double();
      C2 = node.child("epsilon").child("parameters").attribute("C2").as_double();
      B3 = node.child("epsilon").child("parameters").attribute("B3").as_double();
      C3 = node.child("epsilon").child("parameters").attribute("C3").as_double();
      B4 = node.child("epsilon").child("parameters").attribute("B4").as_double();
      C4 = node.child("epsilon").child("parameters").attribute("C4").as_double();
      B5 = node.child("epsilon").child("parameters").attribute("B5").as_double();
      C5 = node.child("epsilon").child("parameters").attribute("C5").as_double();
      result.elmag.initSellmeier_r(B1, C1, B2, C2, B3, C3, B4, C4, B5, C5, aux_mu.real());
    } else
      throw std::runtime_error("Unknown type for epsilon");
  }
  return result;
};

#ifdef OPTIMET_BELOS
void Reader::readParameterList(pugi::xml_document const &root_node) {
  auto const xml_params = root_node.child("ParameterList");
  std::ostringstream str_params;
  if(not xml_params)
    str_params << "<ParameterList name=\"belos\"></ParameterList>";
  else
    xml_params.print(str_params);
  run->belos_params = Teuchos::getParametersFromXmlString(str_params.str());
  if(not run->belos_params->isParameter("Solver"))
    run->belos_params->set("Solver", "scalapack");
}
#else
void Reader::readParameterList(pugi::xml_document const &) {}
#endif
