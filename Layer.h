#pragma once
#include <vector>
#include <cstdlib>
#include "OdeSolver.h"
#include "boundary_value_problem.h"
#include "Parameters.h"
#include "plots.h"
#include <iomanip>
#include <functional>
#include <complex>
#include "OdeSolver.h"

const double pi = 3.1415926538;


//template<class T>
class layer {
	//double kappa;
	//std::complex<double> alpha;
	std::vector<std::function<std::complex<double>(double, const std::vector<std::complex<double>>&, std::complex<double>, double)>>
		_equations;
	std::function<double(double)>& mu;
	std::function<double(double)>& rho;


	std::vector<std::function<std::complex<double>(double, const std::vector<std::complex<double>>&)>> get_equation(std::complex<double> alpha, double kappa, size_t size) const;

public:

	layer(/*std::vector<std::function<std::complex<double>(double, const std::vector<std::complex<double>>&)>>
		equations*/) :/*_equations(equations),*/ mu(Parameters::smooth_params[0]), rho(Parameters::smooth_params[1]) {
		_equations = {
	[=](double x, const std::vector<std::complex<double>>& v, std::complex<double> alpha, double kappa) {return v[1] / mu(x); },
	[=](double x, const std::vector<std::complex<double>>& v, std::complex<double> alpha, double kappa) {return (alpha * alpha * mu(x) - kappa * kappa * rho(x)) * v[0]; },
	[=](double x, const std::vector<std::complex<double>>& v, std::complex<double> alpha, double kappa) {return v[3] / mu(x); },
	[=](double x, const std::vector<std::complex<double>>& v, std::complex<double> alpha, double kappa) {return 2.0 * alpha * v[0] + (alpha * alpha * mu(x) - kappa * kappa * rho(x)) * v[2]; },
	[=](double x, const std::vector<std::complex<double>>& v, std::complex<double> alpha, double kappa) {return v[5] / mu(x); },
	[=](double x, const std::vector<std::complex<double>>& v, std::complex<double> alpha, double kappa) {return 2.0 * v[0] + 4.0 * alpha * v[2] + (alpha * alpha * mu(x) - kappa * kappa * rho(x)) * v[4]; }
		};
	}

	double defOnTop(std::complex<double> alpha, double kappa) const;
	double findRoot(double a, double b, const std::function<double(double)>& f, double epsilon = 0.1e-6);
	std::vector<double> find_roots(double a, double b, const std::function<double(double)>& f, size_t n, double epsilon = 0.1e-6);
	std::map<double, std::vector<double>> dispersionalSet(double max_kappa, double step, size_t n, const std::function<double(double, double)>& f);
	std::complex<double> waves(double x1, double kappa, const std::vector<std::complex<double>>& roots, const std::vector<std::complex<double>>& residuals);
	std::map<double, std::complex<double>> waveField(double a, double b, double kappa, double step, const std::vector<std::complex<double>>& roots, const std::vector<std::complex<double>>& residuals);
	std::complex<double> residual(std::complex<double> alpha, double kappa);
	std::vector<std::complex<double>> residualSet(const std::vector<std::complex<double>>& roots, double kappa);
	//std::vector<std::complex<double>> addKer(std::complex<double> alpha, double kappa, double n, double x1);
	std::vector<std::vector<std::complex<double>>> MatrixV(double kappa, double c, double d, size_t rows, size_t columns);
	//std::vector<std::complex<double>> Kernels(std::vector<std::vector<std::complex<double>>> mat);
	std::vector<std::complex<double>> getRoots(double kappa);
	std::vector<std::complex<double>> observed(const std::vector<double>& points, double kappa, std::vector<std::complex<double>> residual_set, std::vector<std::complex<double>> roots);
};

