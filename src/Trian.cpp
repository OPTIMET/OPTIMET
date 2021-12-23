// (C) University College London 2017
// This file is part of Optimet, licensed under the terms of the GNU Public License

// Optimet is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Optimet is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Optimet. If not, see <http://www.gnu.org/licenses/>.

#include "Tools.h"
#include "Trian.h"
#include <iostream>

Trian::Trian(const double* c1,  const double* c2, const double* c3): 
	coord1(3), coord2(3), coord3(3), cp(3), nvec(3), lvec(3), dl(3), deter(0.0) 
{
         
	 coord1={ *c1, *(c1 + 1), *(c1 + 2) }; 
	 coord2={ *c2, *(c2 + 1), *(c2 + 2) };  
	 coord3={ *c3, *(c3 + 1), *(c3 + 2) };

	calculate_parameters();
}

// calculating some features of the mesh triangle in Cartesian system
void Trian::calculate_parameters(){
            
      for (unsigned int j = 0; j != 3; ++j) {
		cp[j] = (1.0 / 3.0) * (coord1[j] + coord2[j] + coord3[j]); 
	}
	
	std::vector<double> p21(3), p32(3), p13(3);

        double* p21nvec = new double[3];
	double* p32nvec = new double[3];
	double* p13nvec = new double[3];
	
	for (unsigned int j = 0; j != 3; ++j) {
		p21[j] = coord2[j] - coord1[j];
		p13[j] = coord1[j] - coord3[j];
		p32[j] = coord3[j] - coord2[j];
	}

        for (unsigned int j = 0; j != 3; ++j) {
	
	        p21nvec[j] = p21[j]/(Tools::norm(&p21[0]));
		p32nvec[j] = p32[j]/(Tools::norm(&p32[0]));
		p13nvec[j] = p13[j]/(Tools::norm(&p13[0]));
		
		}
	
	lvec[0]= p21nvec;
	lvec[1]= p32nvec;
	lvec[2]= p13nvec;

        dl[0]= Tools::norm(&p21[0]);
	dl[1]= Tools::norm(&p32[0]);
	dl[2]= Tools::norm(&p13[0]);

	double* nvector = new double[3];
		
    Tools::cross(nvector, &p13[0], &p21[0]);
	nvec = { *nvector, *(nvector + 1), *(nvector + 2) };
	deter = Tools::norm(nvector);

	nvec[0] = nvec[0]/deter;
	nvec[1] = nvec[1]/deter;
	nvec[2] = nvec[2]/deter;
	

}
