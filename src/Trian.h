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
