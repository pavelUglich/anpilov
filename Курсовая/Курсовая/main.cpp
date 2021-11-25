#include <iostream>
#include "BoundaryValueProblemSmooth.h"
#include "BoundaryValueProblem.h"
#include "Parameters.h"

using namespace std;

int main()
{
	double innerRadius = 0.001;
	
	BoundaryValueProblemSmooth a([](double x) {return 1.0; }, innerRadius, 0);
	BoundaryValueProblemSmooth b([](double x) {return 0.0; }, innerRadius, 0);

	vector<complex<double>> v1 = { 0,0,0,0 };
	vector<complex<double>> v2 = { 1,0,0,0 };
	vector<complex<double>> v3 = { 0,0,1,0 };

	auto res1 = a.CauchyProblem(innerRadius, 1, v1);
	auto res2 = b.CauchyProblem(innerRadius, 1, v2);
	auto res3 = b.CauchyProblem(innerRadius, 1, v3);

	auto delta = res2[0] * res3[1] - res2[1] * res3[0];
	auto aa = -(res1[0] * res3[1] - res1[1] * res3[0]) / delta;
	auto bb = -(res2[0] * res1[1] - res2[1] * res1[0]) / delta;

	//auto z = aa;

	/*system("Pause");

	delta = res2[0] * res3[2] - res2[2] * res3[0];
	aa = -(res1[0] * res3[2] - res1[2] * res3[0]) / delta;
	bb = -(res2[0] * res1[2] - res2[2] * res1[0]) / delta;*/
	system("Pause");


	return 0;
}