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
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

#ifdef OPTIMET_BELOS
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#endif

using namespace pugi;

namespace optimet {
namespace {
std::shared_ptr<Geometry> read_geometry(pugi::xml_document const &node);
Scatterer read_spherical_scatterer(pugi::xml_node const &node, t_int nMax);
std::shared_ptr<Geometry> read_structure(pugi::xml_node const &inputFile, t_int nMax);
std::shared_ptr<Excitation> read_excitation(pugi::xml_document const &inputFile, t_int nMax);
scalapack::Parameters read_parallel(const pugi::xml_node &node);
#ifdef OPTIMET_BELOS
Teuchos::RCP<Teuchos::ParameterList> read_parameter_list(pugi::xml_document const &root_node);
std::tuple<bool, t_int> read_fmm_input(pugi::xml_node const &node);
#endif
Run simulation_input(pugi::xml_document const &inputFile);

std::shared_ptr<Geometry> read_geometry(pugi::xml_document const &inputFile) {
  // Find the simulation node
  auto const simulation_node = inputFile.child("simulation");
  if(!simulation_node)
    throw std::runtime_error("Simulation parameters not defined!");

  auto const nMax = simulation_node.child("harmonics").attribute("nmax").as_int();

  // Find the geometry node
  auto const geo_node = inputFile.child("geometry");
  if(!geo_node)
    throw std::runtime_error("Geometry not defined!");

  // Check if a structure is defined
  if(geo_node.child("structure"))
    return read_structure(geo_node, nMax);

  auto result = std::make_shared<Geometry>();
  result->structureType = 0;

  // Find all scattering objects
  for(xml_node node = geo_node.child("object"); node; node = node.next_sibling("object"))
    result->pushObject(read_spherical_scatterer(node, nMax));

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
      result->bground.init_r(aux_epsilon, aux_mu);
    }
  }

  // Validate the geometry in the return
  if(result->objects.size() == 0)
    throw std::runtime_error("No scatterers defined in input");
  return result;
}

std::shared_ptr<Geometry> read_structure(xml_node const &geo_node_, t_int nMax) {
  auto geometry = std::make_shared<Geometry>();
  xml_node struct_node = geo_node_.child("structure");

  geometry->structureType = 1; // set the spiral structure flag

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
      geometry->spiralSeparation = (d / 2) - 2 * radius * consFrnmTom;
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

    auto const scatterer = read_spherical_scatterer(struct_node.child("object"), nMax);
    for(int i = 0; i < No - 1; i++) {
      geometry->pushObject(scatterer);
      std::string const normal = struct_node.child("properties").attribute("normal").value();
      if(normal == "x") {
        // x is normal (conversion is x(pol) -> y; y(pol) -> z
        geometry->normalToSpiral = 0;
        geometry->objects.back().vR = Tools::toSpherical({0.0, X[i], Y[i]});
      } else if(normal == "y") {
        // y is normal (conversion is x(pol) -> z; y(pol) -> x
        geometry->normalToSpiral = 1;
        geometry->objects.back().vR = Tools::toSpherical({Y[i], 0, X[i]});
      } else if(normal == "z") {
        // z is normal (conversion is x(pol) -> x; y(pol) -> x
        geometry->normalToSpiral = 2;
        geometry->objects.back().vR = Tools::toSpherical({X[i], Y[i], 0});
      } else
        throw std::runtime_error("Unknown normal " + normal);
    }
  }

  if(geometry->objects.size() == 0)
    throw std::runtime_error("No scatterers defined in input");
  return geometry;
}

Scatterer read_spherical_scatterer(pugi::xml_node const &node, t_int nMax) {
  if(node.attribute("type").value() != std::string("sphere"))
    std::runtime_error("Expecting a spherical scatterer");
  Scatterer result(nMax);
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

std::shared_ptr<Excitation> read_excitation(pugi::xml_document const &inputFile, t_int nMax) {
  // Find the source node
  auto const ext_node = inputFile.child("source");
  if(!ext_node)
    std::runtime_error("Source not defined!");

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
  auto result = std::make_shared<optimet::Excitation>(source_type, Einc, vKinc, nMax);
  result->populate();

  return result;
}

void read_output(pugi::xml_document const &inputFile, Run &run) {
  // Find the source node
  auto const out_node = inputFile.child("output");
  if(!out_node)
    throw std::runtime_error("Output not defined!");

  // Determine type
  if(!std::strcmp(out_node.attribute("type").value(), "coefficients"))
    run.outputType = 2;

  if(!std::strcmp(out_node.attribute("type").value(), "field")) {
    run.outputType = 0;      // Field output requested
    run.singleComponent = 0; // Set this to zero as default

    run.params[0] = out_node.child("grid").child("x").attribute("min").as_double() * 1e-9;
    run.params[1] = out_node.child("grid").child("x").attribute("max").as_double() * 1e-9;
    run.params[2] = out_node.child("grid").child("x").attribute("steps").as_double();
    run.params[3] = out_node.child("grid").child("y").attribute("min").as_double() * 1e-9;
    run.params[4] = out_node.child("grid").child("y").attribute("max").as_double() * 1e-9;
    run.params[5] = out_node.child("grid").child("y").attribute("steps").as_double();
    run.params[6] = out_node.child("grid").child("z").attribute("min").as_double() * 1e-9;
    run.params[7] = out_node.child("grid").child("z").attribute("max").as_double() * 1e-9;
    run.params[8] = out_node.child("grid").child("z").attribute("steps").as_double();

    run.projection =
        !std::strcmp(out_node.child("projection").attribute("spherical").value(), "true");

    run.singleMode = out_node.child("singlemode");

    if(!std::strcmp(out_node.child("singlemode").attribute("dominant").value(), "auto"))
      run.dominantAuto = true;
    else {
      run.dominantAuto = false;
      run.singleModeIndex.init(out_node.child("singlemode").attribute("n").as_int(),
                               out_node.child("singlemode").attribute("m").as_int());
      run.singleComponent =
          !std::strcmp(out_node.child("singlemode").attribute("component").value(), "TE");
    }
  }

  if(!std::strcmp(out_node.attribute("type").value(), "response")) {
    if(out_node.child("scan").child("wavelength")) {
      double lam_start(0.), lam_final(0.);
      lam_start = out_node.child("scan").child("wavelength").attribute("initial").as_double();
      lam_final = out_node.child("scan").child("wavelength").attribute("final").as_double();
      run.params[0] =
          out_node.child("scan").child("wavelength").attribute("initial").as_double() * 1e-9;
      run.params[1] =
          out_node.child("scan").child("wavelength").attribute("final").as_double() * 1e-9;

      int stepsize(0), steps(0);
      stepsize = out_node.child("scan").child("wavelength").attribute("stepsize").as_double();

      // claculate no of steps
      steps = int(lam_final - lam_start) / stepsize;
      run.params[2] = steps + 1;
      run.outputType = 11;
    }

    if(out_node.child("scan").child("radius")) {
      run.params[3] =
          out_node.child("scan").child("radius").attribute("initial").as_double() * 1e-9;
      run.params[4] = out_node.child("scan").child("radius").attribute("final").as_double() * 1e-9;
      run.params[5] = out_node.child("scan").child("radius").attribute("steps").as_double();
      run.outputType = 12;
    }

    if(out_node.child("scan").child("wavelength") && out_node.child("scan").child("radius"))
      run.outputType = 112;
  }
}

scalapack::Parameters read_parallel(const pugi::xml_node &node) {
  scalapack::Parameters result;
  result.block_size = node.attribute("block_size").as_uint(result.block_size);
  result.grid.rows = node.child("grid").attribute("rows").as_uint(result.grid.rows);
  result.grid.cols = node.child("grid").attribute("cols").as_uint(result.grid.cols);
  return result;
}

#ifdef OPTIMET_BELOS
Teuchos::RCP<Teuchos::ParameterList> read_parameter_list(pugi::xml_document const &root_node) {
  auto const xml_params = root_node.child("ParameterList");
  std::ostringstream str_params;
  if(not xml_params)
    str_params << "<ParameterList name=\"belos\"></ParameterList>";
  else
    xml_params.print(str_params);
  auto result = Teuchos::getParametersFromXmlString(str_params.str());
  if(not result->isParameter("Solver"))
    result->set("Solver", "scalapack");
  return result;
}

std::tuple<bool, t_int> read_fmm_input(pugi::xml_node const &node) {
  if(not node)
    return std::make_tuple(false, 1);
  if(not node.attribute("subdiagonals"))
    return std::make_tuple(true, std::numeric_limits<t_int>::max());
  return std::make_tuple(true, node.attribute("subdiagonals").as_int());
}
#endif

Run simulation_input(pugi::xml_document const &inputFile) {
  Run result;
  result.geometry = read_geometry(inputFile);
  result.nMax = result.geometry->nMax();

  // Read Excitation
  result.excitation = read_excitation(inputFile, result.nMax);
  // Update the geometry in case we had dynamic models
  result.geometry->update(result.excitation);

  // Read Excitation
  read_output(inputFile, result);

  result.parallel_params = read_parallel(inputFile.child("parallel"));
#ifdef OPTIMET_BELOS
  result.belos_params = read_parameter_list(inputFile);
  std::tie(result.do_fmm, result.fmm_subdiagonals) = read_fmm_input(inputFile.child("FMM"));
#endif

  return result;
}
}

Run simulation_input(std::string const &fileName_) {
  pugi::xml_document inputFile;
  auto const fileResult = inputFile.load_file(fileName_.c_str());
  if(!fileResult) {
    std::ostringstream msg;
    msg << "Error reading or parsing input file " << fileName_ << "!";
    throw std::runtime_error(msg.str());
  }
  return simulation_input(inputFile);
}

Run simulation_input(std::istream &buffer) {
  pugi::xml_document inputFile;
  auto const fileResult = inputFile.load(buffer);
  if(!fileResult)
    throw std::runtime_error("Error reading or parsing istream input");
  return simulation_input(inputFile);
}
}
