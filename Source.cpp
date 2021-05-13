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
#include "gen/Galgo.hpp"
//#include "Source.h"

std::vector<std::function<double(double)>> Parameters::smooth_params;
std::vector<double> Parameters::const_params;
std::vector<double> Parameters::points;
std::vector<std::vector<double>> Parameters::piecewise_linear_params;
kind_of_solution Parameters::kind;

const double pi = 3.1415926538;
double kappa;

double defOnTop(std::complex<double> alpha, std::complex<double> kappa) {
	const auto mu = Parameters::smooth_params[0];
	const auto rho = Parameters::smooth_params[1];
	OdeSolver<std::complex<double>> cauchy_problem = { {
	[=](double x, const std::vector<std::complex<double>>& v) {return v[1] / mu(x); },
	[=](double x, const std::vector<std::complex<double>>& v) {return (alpha * alpha * mu(x) - kappa * kappa * rho(x)) * v[0]; }
	}
	, 0.1e-6, RUNGE_KUTTA_FELDBERG };
	auto solution = cauchy_problem.solve(0, 1, { 0, 1 });
	return solution[1].real();
}
std::complex<double> residual(std::complex<double> alpha, std::complex<double> kappa) {
	const auto mu = Parameters::smooth_params[0];
	const auto rho = Parameters::smooth_params[1];
	OdeSolver<std::complex<double>> cauchy_problem = { {
	[=](double x, const std::vector<std::complex<double>>& v) {return v[1] / mu(x); },
	[=](double x, const std::vector<std::complex<double>>& v) {return (alpha * alpha * mu(x) - kappa * kappa * rho(x)) * v[0]; },
	[=](double x, const std::vector<std::complex<double>>& v) {return v[3] / mu(x); },
	[=](double x, const std::vector<std::complex<double>>& v) {return 2.0 * alpha * v[0] + (alpha * alpha * mu(x) - kappa * kappa * rho(x)) * v[2]; }
	}, 0.1e-6, RUNGE_KUTTA_FELDBERG };
	auto solution = cauchy_problem.solve(0, 1, { 0, 1, 0, 0 });
	return solution[0] / solution[3];
}


/*std::vector<double> poly(double kappa, double step, const std::function<std::vector<double>(double, double)>& f) {
	std::vector<double> result;
	double a = 0;
	double b = 1.5 * kappa;
	for (int i = a; i < b; i += step) {
		double c = f(i, kappa)[1];
		double d = f(i, kappa)[0];
		result.push_back(c / d);
	}
	return result;
}*/

std::vector<std::complex<double>> expon(double kappa, double step, double x1) {
	std::vector<std::complex<double>> result;
	int counter = 0;
	double a = 0;
	double b = 1.5 * kappa;
	for (double i = a; i < b; i += step) {
		std::complex<double> im = { 0,1 };
		std::complex<double> c = exp(im * i * x1);
		result.push_back(c);
		counter++;
	}
	return result;
}



double findRoot(double a, double b, const std::function<double(double)>& f, double epsilon = 0.1e-6) {
	while (abs(b - a) > epsilon) {
		const double fa = f(a);
		const double fb = f(b);
		a = b - (b - a) * fb / (fb - fa);
		b = a - (a - b) * fa / (fa - fb);
	}
	return b;
}

std::vector<double> find_roots(double a, double b, const std::function<double(double)>& f, size_t n, double epsilon = 0.1e-6) {
	std::vector<double> result;
	auto step = (b - a) / n;
	for (int i = 0; i < n; i += 1) {
		auto left = a + i * step;
		auto right = left + step;
		if (f(left) * f(right) < 0) {
			double f_r = findRoot(left, right, f);
			result.push_back(f_r);
		}
	}
	return result;
}

std::map<double, std::vector<double>> dispersional_set(double max_kappa, double step, size_t n, const std::function<double(double, double)>& f) {
	std::map<double, std::vector<double>> result;
	double kappa = step;
	while (kappa < max_kappa) {
		auto ff = [kappa, f](double alpha) {return f(alpha, kappa);};
		result[kappa] = find_roots(0, max_kappa, ff, n);
		kappa += step;
	}
	return result;
}


std::map<double, std::vector<double>> extended_sys(double max_kappa, double step, size_t n, const std::function<double(double, double)>& f) {
	std::map<double, std::vector<double>> result;
	double kappa = step;
	auto ff = [kappa, f](double alpha) {return f(alpha, kappa);};
	result[kappa] = find_roots(0, 1.5 * kappa, ff, n);
	kappa += 4;
	while (kappa <= max_kappa) {
		auto ff = [kappa, f](double alpha) {return f(alpha, kappa);};
		result[kappa] = find_roots(0, max_kappa, ff, n);
		kappa += 5;
	}
	return result;
}

void show_vector(std::vector<double>& a)
{
	for (std::vector<double>::iterator it = a.begin(); it != a.end(); ++it)
		std::cout << *it;
}

void plotTheDispersionalCurves(std::map<double, std::vector<double>>& dispersionSet, const std::string& fileName)
{
	std::ofstream os(fileName);
	os << "\\begin{tikzpicture}[scale=1.5]\n";
	os << "\\begin{axis}[grid]\n";
	const auto numberOfCurves = dispersionSet.rbegin()->second.size();
	for (auto& item : dispersionSet) {
		sort(item.second.begin(), item.second.end());//,
			//[](auto x, auto y) { return x < y.real(); });
	}
	const auto maxFreq = dispersionSet.rbegin()->first;
	for (size_t i = 0; i < numberOfCurves; i++) {
		os << "\\addplot[smooth, black] plot coordinates{\n";
		for (auto x : dispersionSet) {
			if (i < x.second.size())
			{
				os << "(" << x.second[i] << ", " << x.first << ") ";
			}
		}
		os << "};\n";
	}
	os << "\\end{axis}\n";
	os << "\\end{tikzpicture}\n";
	os.close();
}

std::map<double, std::vector<double>> split_sets(const std::map<double, std::vector<double>>& left, const std::map<double, std::vector<double>>& right) {
	std::map<double, std::vector<double>> result;
	for (auto x : left)
	{
		result[x.first] = x.second;
		auto value = right.at(x.first);
		std::transform(value.begin(), value.end(), value.begin(), [](auto x) {return -x;});
		result[x.first].insert(result[x.first].end(), value.begin(), value.end());
	}
	return result;
}

std::vector<std::complex<double>> get_roots(double kappa) {
	std::vector<std::complex<double>> result;
	auto fun = [=](double alpha) {
		return defOnTop({ alpha,0 }, kappa);
	};
	auto fun_imag = [=](double alpha) {
		return defOnTop({ 0,alpha }, kappa);
	};
	// вычисляем два вектора корней и сливаем в один
	auto real_roots = find_roots(0, 1.5 * kappa, fun, 20);
	for (int i = 0;i < real_roots.size();i++) {
		std::complex<double> re = {real_roots[i], 0};
		result.push_back(re);
	}
	auto imag_roots = find_roots(0, 50, fun_imag, 20);
	for (int i = 0;i < imag_roots.size();i++) {
		std::complex<double> im = { 0, imag_roots[i] };
		result.push_back(im);
	}
	return result;
}

std::vector<std::complex<double>> reses_set(const std::vector<std::complex<double>> & roots, double kappa) {
	std::vector<std::complex<double>> result;
	for (int i = 0;i < roots.size();i++) {
		result.push_back(residual(roots[i], kappa));
	}
	return result;
}

std::complex<double> waves(double x1, double kappa, const std::vector<std::complex<double>> & roots, const std::vector<std::complex<double>>& residuals) {
	std::complex<double> result = 0;
	std::complex<double> im = { 0,1 };
	for (int i = 0; i < roots.size(); i++) {
		result += im * residuals[i] * exp(im * roots[i] * x1);
	}
	return result;
}

std::map<double, std::complex<double>> wave_field(double a, double b, double kappa, double step, const std::vector<std::complex<double>>& roots, const std::vector<std::complex<double>>& residuals) {
	std::map<double, std::complex<double>> result;
	for (double i = a; i < b; i+=step) {
		result[i] = waves(i, kappa, roots, residuals);
	}
	return result;
}

void plotTheWaveField(const std::map<double, std::complex<double>>& waveField, const std::string& fileName)
{
	std::ofstream os(fileName);
	os << "\\begin{tikzpicture}[scale=1.5]\n";
	os << "\\begin{axis}[grid]\n";
	os << "\\addplot[smooth, red] plot coordinates{\n";
	for (auto x : waveField) {
		os << "(" << x.first << ", " << x.second.real()   << ") ";
	}
	os << "};\n";
	os << "\\addplot[smooth, blue] plot coordinates{\n";
	for (auto x : waveField) {
		os << "(" << x.first << ", " << x.second.imag() << ") ";
	}
	os << "};\n";
	os << "\\end{axis}\n";
	os << "\\end{tikzpicture}\n";
	os.close();
}


// 1. Найти наблюдаемое поле std::vector<std::complex<double>> observed(const std::vector<double> & points, double kappa)
// целевая функция (double objective(const std::vector<double> & parameters, double kappa, const std::vector<std::complex<double>> & observed))
// 
// 2. set_parameters(const std::vector<double> & parameters) // 
// меняем 
// Parameters::smooth_params[0], Parameters::smooth_params[1]
// 3. Находим поле перемещений, соответствующее параметрам parameters в тех же точках points
// 4. Вычисляем целевую функцию

std::vector<std::complex<double>> observed(const std::vector<double>& points, double kappa) {
	std::vector<std::complex<double>> result;
	auto fun = [=](double alpha, double kappa) {
		return defOnTop({ alpha,0 }, kappa);
	};
	auto fun_imag = [=](double alpha, double kappa) {
		return defOnTop({ 0,alpha }, kappa);
	};
	auto u = [=](double alpha, double kappa) {
		return defOnTop({ alpha,0 }, kappa);
	};
	const auto roots = get_roots(kappa);
	const auto reses = reses_set(roots, kappa);
	for (auto point : points)
	{
		result.push_back(waves(point, kappa, roots, reses));
	}
	return result;
}

void set_parameters(const std::vector<double>& parameters) {
	Parameters::smooth_params.clear();
	Parameters::smooth_params.push_back([=](double x) {return 1 + parameters[0] * x + parameters[1] * (1 - x);});
	Parameters::smooth_params.push_back([=](double x) {return 1 + parameters[2] * x + parameters[3] * (1 - x);});
}

double objective(const std::vector<double>& parameters, double kappa, const std::vector<std::complex<double>>& observe, const std::vector<double>& points) {
	double result;
	set_parameters(parameters);
	auto pf = observed(points, kappa);
	for (size_t i = 0; i < pf.size(); i++)
	{
		result += abs(pf[i] - observe[i]) * abs(pf[i] - observe[i]);
	}
	return result;
}


// objective class example
template <typename T>
class MyObjective
{
public:
	// objective function example : Rosenbrock function
	// minimizing f(x,y) = (1 - x)^2 + 100 * (y - x^2)^2

	// ����� ���� ������� �������
	static std::vector<T> Objective(const std::vector<T>& x)
	{
		T obj = -(pow(1 - x[0], 2) + 100 * pow(x[1] - x[0] * x[0], 2));
		return { obj };
	}
	// NB: GALGO maximize by default so we will maximize -f(x,y)
};

// constraints example:
// 1) x * y + x - y + 1.5 <= 0
// 2) 10 - x * y <= 0
template <typename T>
std::vector<T> MyConstraint(const std::vector<T>& x)
{
	return { x[0] * x[1] + x[0] - x[1] + 1.5,10 - x[0] * x[1] };
}

int main()
{
	setlocale(0, "");
	//std::cout << defOnTop({ 0, 1.0 }, 2.1) << std::endl;
	// построить дисперсионные кривые для разгых законов изменения неоднородности (1+sin(pi*x), 1+exp(-x))
	Parameters::smooth_params.push_back([](double x) {return 1 + exp(-x);});//1+0.5*sin(Pi*x)//1+x//1+exp(-x)
	Parameters::smooth_params.push_back([](double x) {return 1;});
	// defining lower bound LB and upper bound UB
	std::vector<double> LB({ 0.0,0.0 });
	std::vector<double> UB({ 1.0,13.0 });

	// initiliazing genetic algorithm
	galgo::GeneticAlgorithm<double> ga(MyObjective<double>::Objective, 100, LB, UB, 50, true);

	// setting constraints
	ga.Constraint = MyConstraint;

	// running genetic algorithm
	ga.run();
	double kappa = 10.0;
	auto fun = [=](double alpha, double kappa) {
		return defOnTop({ alpha,0 }, kappa);
	};
	auto fun_imag = [=](double alpha, double kappa) {
		return defOnTop({ 0,alpha }, kappa);
	};
	auto u = [=](double alpha, double kappa) {
		return defOnTop({ alpha,0 }, kappa);
	};
	auto set_re = dispersional_set(10, 0.2, 20, fun);
	auto set_im = dispersional_set(10, 0.2, 20, fun_imag);
	//auto set_re_1 = extended_sys(10, 1, 20, fun);
	//auto set_im_1 = extended_sys(10, 1, 20, fun_imag);
	auto common_set = split_sets(set_re, set_im);
	//auto common_set_1 = split_sets(set_re, set_im);
	plotTheDispersionalCurves(common_set, "1text.txt");
	const auto roots = get_roots(kappa);
	const auto reses = reses_set(roots, kappa);
	auto wf = wave_field(0, 1, 0, 0.01, roots, reses);
	plotTheWaveField(wf, "5text.txt");
	system("pause");

}
