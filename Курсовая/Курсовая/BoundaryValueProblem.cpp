#include "BoundaryValueProblem.h"
#include "Parameters.h"
#include<vector>

using namespace std;

complex<double> BoundaryValueProblem::IntegrandDenumerator(size_t size) const {
	auto CauchyProblemSolution = CauchyProblemSolutions(1, size);
	if (size == 2) return CauchyProblemSolution[1];
	if (size == 4) return CauchyProblemSolution[1] / CauchyProblemSolution[3];
	if (size == 6) return CauchyProblemSolution[5];
	return 0.0;
}

complex<double> BoundaryValueProblem::waveField() const
{
	auto cauchyProblemSolution = CauchyProblemSolutions(1, 2);
	return cauchyProblemSolution[0] / cauchyProblemSolution[1];
}

void BoundaryValueProblem::SetAlpha(complex<double> _alpha) {
	//alpha = _alpha;
}

double BoundaryValueProblem::getKappa() const
{
	return kappa;
}

double BoundaryValueProblem::getq() const
{
	return 0;
}


complex<double> BoundaryValueProblem::NewtonMethod() {
/*	complex<double> a = IntegrandDenumerator(4);
	while (abs(a) > eps) {
		alpha -= a;
		a = IntegrandDenumerator(4);
	}*/
	return {0,0};//alpha;
}


BoundaryValueProblem::BoundaryValueProblem()
{
}

BoundaryValueProblem::BoundaryValueProblem(double innerRadius, double kappa, double eps) :
	innerRadius(innerRadius), 
	eps(eps), 	
	kappa(kappa)
{
}

BoundaryValueProblem::~BoundaryValueProblem() = default;
