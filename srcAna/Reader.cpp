
// (C) University College London 2017
// This file is part of Optimet, licensed under the terms of the GNU Public License
//
// Optimet is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Optimet is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Optimet. If not, see <http://www.gnu.org/licenses/>.

#include "Reader.h"
#include <iostream>
#include <fstream>
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
Scatterer read_scatterer(pugi::xml_node const &node, t_int nMax, t_int nMaxS);
std::shared_ptr<Geometry> read_structure(pugi::xml_node const &inputFile, t_int nMax, t_int nMaxS, bool ACA_cond);
std::shared_ptr<Excitation> read_excitation(pugi::xml_document const &inputFile, t_int nMax, ElectroMagnetic &bground);
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

  bool ACA_cond; //condition for ACA
  if(!std::strcmp(simulation_node.child("ACA").attribute("compression").value(), "yes"))
  ACA_cond = true;
  else
  ACA_cond = false;
  auto const nMax = simulation_node.child("harmonics").attribute("nmax").as_int();
  auto const nMaxS = 1 * nMax; //number of harmonics for SH expansion
  // Find the geometry node
  auto const geo_node = inputFile.child("geometry");
  if(!geo_node)
    throw std::runtime_error("Geometry not defined!");

  // Check if a structure is defined
  if(geo_node.child("structure"))
    return read_structure(geo_node, nMax, nMaxS, ACA_cond);

  auto result = std::make_shared<Geometry>();
  result->ACAcompression(ACA_cond);
  // Find all scattering objects
  for(xml_node node = geo_node.child("object"); node; node = node.next_sibling("object"))
    result->pushObject(read_scatterer(node, nMax, nMaxS));

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
      result->bground.init_r(aux_epsilon, aux_mu, 0.0, 0.0, 0.0, 0.0);
    }
  }
  // Validate the geometry in the return
  if(result->objects.size() == 0)
    throw std::runtime_error("No scatterers defined in input");
  return result;
}

std::shared_ptr<Geometry> read_structure(xml_node const &geo_node_, t_int nMax, t_int nMaxS, bool ACA_cond) {
  auto geometry = std::make_shared<Geometry>();
  xml_node struct_node = geo_node_.child("structure");

  // Add the background properties
  if(geo_node_.child("background")) {
    if(std::strcmp(geo_node_.child("background").attribute("type").value(), "relative")) {
      std::complex<double> aux_epsilon(1.0, 0.);
      std::complex<double> aux_mu(1.0, 0.);
      aux_epsilon.real(
          geo_node_.child("background").child("epsilon").attribute("value.real").as_double());
      aux_epsilon.imag(
          geo_node_.child("background").child("epsilon").attribute("value.imag").as_double());
      aux_mu.real(geo_node_.child("background").child("mu").attribute("value.real").as_double());
      aux_mu.imag(geo_node_.child("background").child("mu").attribute("value.imag").as_double());
      geometry->bground.init_r(aux_epsilon, aux_mu, 0.0, 0.0, 0.0, 0.0);
    }
  }
  
  
  //constructs a cubic cluster of spheres
  if(!std::strcmp(struct_node.attribute("type").value(), "cube")) {
    // Build a cube of spherical scatterers with a corner in origin
    double d; //  distance between spheres
    
    int  No, Ntot;  // number of objects per side of the cube // total number of particles in cube

      No = struct_node.child("properties").attribute("points").as_int();
      
      Ntot = std::pow(No , 3);
 
      d = struct_node.child("properties").attribute("distance").as_double() * consFrnmTom;
  

    // Assign properties to the Scatterer work_object
    if(struct_node.child("object").child("properties").attribute("radius")) {
      auto const radius =
          struct_node.child("object").child("properties").attribute("radius").as_double();
    }

    // Create vectors for X, Y and Z coordinates
    std::vector<double> X(Ntot); // to store x coordinates for all particles in cube
    std::vector<double> Y(Ntot); // to store y coordinates for all particles in cube
    std::vector<double> Z(Ntot); // to store z coordinates for all particles in cube
    
    int i = 0, j = 0, k = 0, br = 0;
   
    for(k = 0; k < No; k++) {
    
     for(j = 0; j < No; j++) {
      
      for(i = 0; i < No; i++) {

        X[br] = d * double(i);
        Y[br] = d * double(j);
        Z[br] = d * double(k);
        
        
        br++;
       }
     }
   }
    

    // convert to a spherical object and push

    auto const scatterer = read_scatterer(struct_node.child("object"), nMax, nMaxS);
    
    geometry->pushObject(scatterer);
   
    for(int i = 1; i < Ntot; i++) {
    
        geometry->objects.back().vR = Tools::toSpherical({X[i], Y[i], Z[i]});
 
        geometry->pushObject(scatterer);
       
  }  
 geometry->ACAcompression(ACA_cond);
}

//assembly of spheres into a zincblende structure (GaAs) with unit cells forming a surface in xy plane

  if(!std::strcmp(struct_node.attribute("type").value(), "GaAssurf")) {
 // Build a noncentrosymmetric surface lattice of GaAs
  double a, R; //  size of cubic unit cell of GaAs
  int  No, Ntot;  // number of cells per side // total number of particles in assembly

  No = struct_node.child("properties").attribute("No_cells").as_int();
  a = struct_node.child("properties").attribute("cell_size").as_double()* consFrnmTom;

   // Assign properties to the Scatterer work_object
     if(struct_node.child("object").child("properties").attribute("radius")) {
     R = struct_node.child("object").child("properties").attribute("radius").as_double()* consFrnmTom ;
     }

   // Create vectors for X1, Y1 and Z1 coordinates
        std::vector<double> X1(18); // to store x coordinates for all particles in first cell
        std::vector<double> Y1(18); // to store y coordinates for all particles in first cell
        std::vector<double> Z1(18); // to store z coordinates for all particles in first cell
        std::vector<double> radii; // radii of atoms in non-centrosymmetric lattice

   // Create vectors for X, Y and Z coordinates
           std::vector<double> X; // to store x coordinates for all particles in structure
           std::vector<double> Y; // to store y coordinates for all particles in structure
           std::vector<double> Z; // to store z coordinates for all particles in structure

    X1 = {a, a, a/2.0, a/2.0, 0.0, 0.0, a/2.0, a, a/2.0, 0.0, 0.0, a, a, 0.0, 3.0*a/4.0, a/4.0, 3.0*a/4.0, a/4.0};
    
    Y1 = {0.0, a/2.0, a/2.0, 0.0, a/2.0, a, a, a, a/2.0, 0.0, 0.0, 0.0, a, a, a/4.0, 3.0*a/4.0, 3.0*a/4.0, a/4.0};
    
    Z1 = {0.0, a/2.0, 0.0, a/2.0, a/2.0, 0.0, a/2.0, a, a, a, 0.0, a, 0.0, a, a/4.0, a/4.0, 3.0*a/4.0, 3.0*a/4.0};

   int i, j, k;
   
      
      for(i = 0; i < No; i++) {
      
      for(j = 0; j < No; j++) {
      
       for(k = 0; k < 18; k++) {
        
        if ((i==0)&&(j==0)) {
        
        
        X.push_back(X1[k]);
        
        Y.push_back(Y1[k]);
        
        Z.push_back(Z1[k]);
        
        if(k<14)
        radii.push_back(R);
        else
        radii.push_back(R/2.0);
             
        }
        
        else if ((j>0)&&(i==0))
        {

        if(Y1[k] > 0.0){
        
        X.push_back(X1[k]);
        Y.push_back(Y1[k] + a * double(j));
        Z.push_back(Z1[k]);
        
        if(k<14)
        radii.push_back(R);
        else
        radii.push_back(R/2.0);

        }
 
       }
       
          else if ((j==0)&&(i>0))
        {

        if(X1[k] > 0.0){
        
        X.push_back(X1[k] + a * double(i));
        Y.push_back(Y1[k]);
        Z.push_back(Z1[k]);
        
        if(k<14)
        radii.push_back(R);
        else
        radii.push_back(R/2.0);

        }
 
       }
       
          else if ((j>0)&&(i>0))
        {

        if((X1[k] > 0.0)&&(Y1[k] > 0.0)){
        
        X.push_back(X1[k] + a * double(i));
        Y.push_back(Y1[k] + a * double(j));
        Z.push_back(Z1[k]);
        
        if(k<14)
        radii.push_back(R);
        else
        radii.push_back(R/2.0);

        }
       } 
    }
  }
 }

 Ntot = X.size();

    // convert to a spherical object and push
 auto const scatterer = read_scatterer(struct_node.child("object"), nMax, nMaxS);
    
    geometry->pushObject(scatterer);
  
    for(int i = 0; i < Ntot; i++) {
        
        if((X[i] == 0.0)&&(Y[i] == 0.0)&&(Z[i] == 0.0))
        continue;
        
        geometry->objects.back().radius = radii[i];
        geometry->objects.back().vR = Tools::toSpherical({X[i], Y[i], Z[i]});
 
        geometry->pushObject(scatterer);
       
  }  
geometry->ACAcompression(ACA_cond);
}
   
//assembly of spheres into a zincblende structure (GaAs) with unit cells forming a cube in first quadrant
if(!std::strcmp(struct_node.attribute("type").value(), "GaAscube")) {
    // Build a noncentrosymmetric cubic lattice of GaAs

   double a, R, ratio; //  size of cubic unit cell of GaAs// radius of Arsenic atom in metaparticle
   int  No, Ntot;  // number of cells per side // total number of particles in assembly

   No = struct_node.child("properties").attribute("No_cells").as_int();
   a = struct_node.child("properties").attribute("cell_size").as_double()* consFrnmTom; 
   
   // Assign properties to the Scatterer work_object
    if(struct_node.child("object").child("properties").attribute("radius")) {
      R =
          struct_node.child("object").child("properties").attribute("radius").as_double()* consFrnmTom;
    }


    // Create vectors for X1, Y1 and Z1 coordinates
    std::vector<double> X1(18); // to store x coordinates for all particles in first cell
    std::vector<double> Y1(18); // to store y coordinates for all particles in first cell
    std::vector<double> Z1(18); // to store z coordinates for all particles in first cell
    std::vector<double> radii; // radii of atoms in non-centrosymmetric lattice

   // Create vectors for X, Y and Z coordinates
     std::vector<double> X; // to store x coordinates for all particles in structure
     std::vector<double> Y; // to store y coordinates for all particles in structure
     std::vector<double> Z; // to store z coordinates for all particles in structure
    
    X1 = {a, a, a/2.0, a/2.0, 0.0, 0.0, a/2.0, a, a/2.0, 0.0, 0.0, a, a, 0.0, 3.0*a/4.0, a/4.0, 3.0*a/4.0, a/4.0};
    
    Y1 = {0.0, a/2.0, a/2.0, 0.0, a/2.0, a, a, a, a/2.0, 0.0, 0.0, 0.0, a, a, a/4.0, 3.0*a/4.0, 3.0*a/4.0, a/4.0};
    
    Z1 = {0.0, a/2.0, 0.0, a/2.0, a/2.0, 0.0, a/2.0, a, a, a, 0.0, a, 0.0, a, a/4.0, a/4.0, 3.0*a/4.0, 3.0*a/4.0};
    
    ratio = 1.7;
    
    int i, j, k, t;

    for(t = 0; t < No; t++) {
      
      for(i = 0; i < No; i++) {
      
      for(j = 0; j < No; j++) {
      
       for(k = 0; k < 18; k++) {
        
        if ((i==0)&&(j==0)&&(t==0)) {
        
        
        X.push_back(X1[k]);
        
        Y.push_back(Y1[k]);
        
        Z.push_back(Z1[k]);

        if(k<14)
        radii.push_back(R);
        else
        radii.push_back(R/ratio);
   
        }
        
        else if ((j>0)&&(i==0)&&(t==0))
        {

        if(Y1[k] > 0.0){
        
        X.push_back(X1[k]);
        Y.push_back(Y1[k] + a * double(j));
        Z.push_back(Z1[k]);

        if(k<14)
        radii.push_back(R);
        else
        radii.push_back(R/ratio);
        
        }
 
       }

     else if ((j==0)&&(i>0)&&(t==0))
        {

        if(X1[k] > 0.0){
        
        X.push_back(X1[k] + a * double(i));
        Y.push_back(Y1[k]);
        Z.push_back(Z1[k]);
        
        if(k<14)
        radii.push_back(R);
        else
        radii.push_back(R/ratio); 
        }
 
       }
       
          else if ((j>0)&&(i>0)&&(t==0))
        {

        if((X1[k] > 0.0)&&(Y1[k] > 0.0)){
        
        X.push_back(X1[k] + a * double(i));
        Y.push_back(Y1[k] + a * double(j));
        Z.push_back(Z1[k]);
        
        if(k<14)
        radii.push_back(R);
        else
        radii.push_back(R/ratio);
        }
       }
       
     else if ((j==0)&&(i==0)&&(t>0))
        {

        if(Z1[k] > 0.0){
        
        X.push_back(X1[k]);
        Y.push_back(Y1[k]);
        Z.push_back(Z1[k]+ a * double(t));
        
        if(k<14)
        radii.push_back(R);
        else
        radii.push_back(R/ratio);
        }
       }  
     else if ((j>0)&&(i==0)&&(t>0))
        {

        if((Z1[k] > 0.0)&&(Y1[k] > 0.0)){
        
        X.push_back(X1[k]);
        Y.push_back(Y1[k]+ a * double(j));
        Z.push_back(Z1[k]+ a * double(t));
        
        if(k<14)
        radii.push_back(R);
        else
        radii.push_back(R/ratio);      
        }
       }   
       
       else if ((j==0)&&(i>0)&&(t>0))
        {

        if((Z1[k] > 0.0)&&(X1[k] > 0.0)){
        
        X.push_back(X1[k]+ a * double(i));
        Y.push_back(Y1[k]);
        Z.push_back(Z1[k]+ a * double(t));
        
        if(k<14)
        radii.push_back(R);
        else
        radii.push_back(R/ratio);
        }
       }   
       
       else if ((j>0)&&(i>0)&&(t>0))
        {

        if((Z1[k] > 0.0)&&(X1[k] > 0.0)&&(Y1[k] > 0.0)){
        
        X.push_back(X1[k]+ a * double(i));
        Y.push_back(Y1[k] + a * double(j));
        Z.push_back(Z1[k]+ a * double(t));
        
        if(k<14)
        radii.push_back(R);
        else
        radii.push_back(R/ratio);      
        }
       }  
       
    }
 }
 }
 }
 
 Ntot = X.size();

 // convert to a spherical object and push
 
 auto const scatterer = read_scatterer(struct_node.child("object"), nMax, nMaxS);
    
    geometry->pushObject(scatterer);
  
    for(int i = 0; i < Ntot; i++) {
        
        if((X[i] == 0.0)&&(Y[i] == 0.0)&&(Z[i] == 0.0))
        continue;
        
        geometry->objects.back().radius = radii[i];       
        geometry->objects.back().vR = Tools::toSpherical({X[i], Y[i], Z[i]});
 
        geometry->pushObject(scatterer);
       
   }  
  geometry->ACAcompression(ACA_cond);
}

 // constructs a cluster of spheres on surface (part of metasurface)
  if(!std::strcmp(struct_node.attribute("type").value(), "surface")) {
    // Build a metasurface of spherical scatterers with a corner in origin
    double d; //  distance between spheres
    
    int  No, Ntot;  // number of objects per side of the metasurface // total number of particles in metasurface

      No = struct_node.child("properties").attribute("points").as_int();
      
      Ntot = std::pow(No , 2);
 
      d = struct_node.child("properties").attribute("distance").as_double() * consFrnmTom;
  

    // Assign properties to the Scatterer work_object
    if(struct_node.child("object").child("properties").attribute("radius")) {
      auto const radius =
          struct_node.child("object").child("properties").attribute("radius").as_double();
    }

    // Create vectors for X , Y and Z coordinates
    std::vector<double> X(Ntot); // to store x coordinates for all particles in surface
    std::vector<double> Y(Ntot); // to store y coordinates for all particles in surface
    std::vector<double> Z(Ntot); // to store z coordinates for all particles in surface, for now 0.0
    
    int i = 0, j = 0, br = 0;
    
     for(j = 0; j < No; j++) {
      
      for(i = 0; i < No; i++) {

        X[br] = d * double(i);
        Y[br] = d * double(j);
        Z[br] = 0.0;

        br++;
       }
     }
    

    // convert to a spherical object and push

    auto const scatterer = read_scatterer(struct_node.child("object"), nMax, nMaxS);
    
    geometry->pushObject(scatterer);
   
    for(int i = 1; i < Ntot; i++) {
    
        geometry->objects.back().vR = Tools::toSpherical({X[i], Y[i], Z[i]});
 
        geometry->pushObject(scatterer);
       
  }  
 geometry->ACAcompression(ACA_cond);
}



  if(geometry->objects.size() == 0)
    throw std::runtime_error("No scatterers defined in input");
  return geometry;


}

Scatterer read_scatterer(pugi::xml_node const &node, t_int nMax, t_int nMaxS) {
  Scatterer result(nMax, nMaxS);
  if(node.attribute("type").value() == std::string("sphere")){
   result.scatterer_type = "sphere";

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
    result.vR = {0.0, 0.0, 0.0};

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
      std::complex<double> epsilon_SH(node.child("epsilon_SH").attribute("value.real").as_double(),
                                   node.child("epsilon_SH").attribute("value.imag").as_double());

   std::complex<double> ksippp(node.child("ksippp").attribute("value.real").as_double(),  
                                   node.child("ksippp").attribute("value.imag").as_double()); 
       
       std::complex<double> ksiparppar(node.child("ksiparppar").attribute("value.real").as_double(), 
                                    node.child("ksiparppar").attribute("value.imag").as_double()); 
       
       std::complex<double> gamma(node.child("gamma").attribute("value.real").as_double(),   
                                  node.child("gamma").attribute("value.imag").as_double()); 
                                                                
      result.elmag.init_r(epsilon, aux_mu, epsilon_SH, ksippp, ksiparppar, gamma);
    }

    else if(node.child("epsilon").attribute("type").value() == std::string("HydroModel")) {
      // Hydrodynamic model (Sipe or Bachelier)  
      std::complex<double> const a_SH(node.child("epsilon").child("parameters").attribute("a.real").as_double(),
          node.child("epsilon").child("parameters").attribute("a.imag").as_double());
      std::complex<double> const b_SH(
          node.child("epsilon").child("parameters").attribute("b.real").as_double(),
           node.child("epsilon").child("parameters").attribute("b.imag").as_double());
      std::complex<double> const d_SH(
          node.child("epsilon").child("parameters").attribute("d.real").as_double(),
           node.child("epsilon").child("parameters").attribute("d.imag").as_double());     
                                                                                       
   
      result.elmag.init_r(0.0, aux_mu, 0.0, 0.0, 0.0, 0.0);
      result.elmag.initHydrodynamicModel_r(a_SH, b_SH, d_SH, aux_mu);
    }

   else if(node.child("epsilon").attribute("type").value() == std::string("SiliconModel")) {
      // Model for Silicon from Schinke valid from 0.5 -1.45um           
      // Surface and bulk tensor values taken as constants

     result.elmag.init_r(0.0, aux_mu, 0.0, 0.0, 0.0, 0.0);
     result.elmag.initSiliconModel_r(aux_mu);  
 }
    
      else
      throw std::runtime_error("Unknown type for epsilon");
  }
  }
// non-spherical particle mesh analysis
else if(node.attribute("type").value() == std::string("arbitrary.shape")){
  result.scatterer_type = "arbitrary.shape";
  
  // Assign coordinates to the center of the object
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
    result.vR = {0.0, 0.0, 0.0};

     
    // Assign properties to the Scatterer 
    if(node.child("properties").attribute("radius"))
    result.radius = node.child("properties").attribute("radius").as_double() * consFrnmTom; //radius of the circumscribed sphere

    // reading the mesh data
    std::ifstream file1("coord.txt"); // coordinates of the vertices
    std::ifstream file2("topol.txt"); // vertices forming a triangle in a counterclockwise manner
 
        std::vector<double> coord;
	std::vector<int> topol;

	double num1 = 0;
	int num2 = 0;
	
	while (file1 >> num1) {
		coord.emplace_back(num1);
		
	}

	while (file2 >> num2) {
		topol.emplace_back(num2);
	}


	
	file1.clear();
	file1.seekg(0, std::ios::beg);
	
	file2.clear();
	file2.seekg(0, std::ios::beg);  // go to the beginning of files
        

	for (int i = 0; i < topol.size(); ++i) {
		topol[i] = topol[i] - 1;
    
	}
        
        result.Mesh(coord, topol);

       if(node.child("epsilon") || node.child("mu")) {
    
    if(node.child("mu").attribute("type").value() != std::string("relative"))
      throw std::runtime_error("The type for mu must be \"relative\"");

    std::complex<double> const aux_mu(node.child("mu").attribute("value.real").as_double(),
                                      node.child("mu").attribute("value.imag").as_double());

    
    if(node.child("epsilon").attribute("type").value() == std::string("relative")) {
      // Static values
      std::complex<double> epsilon(node.child("epsilon").attribute("value.real").as_double(),
                                   node.child("epsilon").attribute("value.imag").as_double());
      std::complex<double> epsilon_SH(node.child("epsilon_SH").attribute("value.real").as_double(),
                                   node.child("epsilon_SH").attribute("value.imag").as_double());
                                   
       std::complex<double> ksippp(node.child("ksippp").attribute("value.real").as_double(), 
                                    node.child("ksippp").attribute("value.imag").as_double()); 
       
       std::complex<double> ksiparppar(node.child("ksiparppar").attribute("value.real").as_double(), 
                                    node.child("ksiparppar").attribute("value.imag").as_double()); 
       
       std::complex<double> gamma(node.child("gamma").attribute("value.real").as_double(),   
                                  node.child("gamma").attribute("value.imag").as_double()); 
                                                                
      result.elmag.init_r(epsilon, aux_mu, epsilon_SH, ksippp, ksiparppar, gamma);
      
      
    } 
    else if(node.child("epsilon").attribute("type").value() == std::string("HydroModel")) {
      // Hydrodynamic model (Sipe or Bachelier)
     std::complex<double> const a_SH(node.child("epsilon").child("parameters").attribute("a.real").as_double(),
          node.child("epsilon").child("parameters").attribute("a.imag").as_double());
      std::complex<double> const b_SH(
          node.child("epsilon").child("parameters").attribute("b.real").as_double(),
           node.child("epsilon").child("parameters").attribute("b.imag").as_double());
      std::complex<double> const d_SH(
          node.child("epsilon").child("parameters").attribute("d.real").as_double(),
           node.child("epsilon").child("parameters").attribute("d.imag").as_double());     
                                                                                       
   
      result.elmag.init_r(0.0, aux_mu, 0.0, 0.0, 0.0, 0.0);
      result.elmag.initHydrodynamicModel_r(a_SH, b_SH, d_SH, aux_mu);
    }
    
    else if(node.child("epsilon").attribute("type").value() == std::string("SiliconModel")) {
      // Model for Silicon from Schinke valid from 0.5 -1.45um
      result.elmag.init_r(0.0, aux_mu, 0.0, 0.0, 0.0, 0.0);
     result.elmag.initSiliconModel_r(aux_mu);

    }
    
      else
      throw std::runtime_error("Unknown type for epsilon");
  }
 
}

  return result;
};

std::shared_ptr<Excitation> read_excitation(pugi::xml_document const &inputFile, t_int nMax, ElectroMagnetic &bground) {
  // Find the source node
  auto const ext_node = inputFile.child("source");
  if(!ext_node)
    std::runtime_error("Source not defined!");

  int source_type;
  double wavelength;
  bool SH_cond;
  SphericalP<std::complex<double>> Einc(std::complex<double>(0.0, 0.0),
                                        std::complex<double>(0.0, 0.0),
                                        std::complex<double>(0.0, 0.0));
  Spherical<double> vKinc(0.0, 0.0, 0.0);
 
  // Determine source type
  if(!std::strcmp(ext_node.attribute("type").value(), "planewave"))
    source_type = 0;
  else // Default is always planewave
    source_type = 0;

   // Check if there are SH sources 
   if(!std::strcmp(ext_node.child("SHsources").attribute("condition").value(), "yes"))
  SH_cond = true;
  else
  SH_cond = false;

  std::complex<double> bgcoeff = std::sqrt(bground.epsilon_r * bground.mu_r);
  // Determine wavelength
  wavelength = ext_node.child("wavelength").attribute("value").as_double();
  wavelength *= 1e-9;
  // Determine propagation values
  vKinc = Spherical<double>(
      (2 * consPi / wavelength),
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
  auto result = std::make_shared<optimet::Excitation>(source_type, Einc, SH_cond, vKinc, nMax, bgcoeff);
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
  result.nMaxS = result.geometry->nMaxS();
  ElectroMagnetic bground =  result.geometry->bground;
  // Read Excitation
  result.excitation = read_excitation(inputFile, result.nMax, bground);
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
