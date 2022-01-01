// (C) University College London 2017
// This file is part of Optimet, licensed under the terms of the GNU Public License
//  Optimet is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// Optimet is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with Optimet. If not, see <http://www.gnu.org/licenses/>.


#ifndef OPTIMET_TRIAN_H
#define OPTIMET_TRIAN_H

#include <vector>

class Trian{
private:

	std::vector<double> coord1; 
	std::vector<double> coord2; 
	std::vector<double> coord3; 
	std::vector<double> nvec;
        std::vector<double> cp;
	std::vector<double> dl;
	std::vector<double*> lvec;
	double deter;
	void calculate_parameters();    
	
public:
	Trian(const double* c1, const double* c2, const double* c3);
	const std::vector<double>& getnorm() const { return nvec; }
        const std::vector<double>& getcp() const { return cp; }
	const std::vector<double>& getdl() const { return dl; }
	const std::vector<double*> getlvec() const { return lvec; }
	double getDeter() const { return deter; }

};

#endif	
