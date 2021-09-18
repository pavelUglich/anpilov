#include "Layer.h"
#include "gen/Galgo.hpp"
//#include "Source.h"

std::vector<std::function<double(double)>> Parameters::smooth_params;
std::vector<double> Parameters::const_params;
std::vector<double> Parameters::points;
std::vector<std::vector<double>> Parameters::piecewise_linear_params;
kind_of_solution Parameters::kind;

//const double pi = 3.1415926538;
double kappa;

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

/*
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
}*/

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
/*
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
		std::complex<double> re = { real_roots[i], 0 };
		result.push_back(re);
	}
	auto imag_roots = find_roots(0, 50, fun_imag, 20);
	for (int i = 0;i < imag_roots.size();i++) {
		std::complex<double> im = { 0, imag_roots[i] };
		result.push_back(im);
	}
	return result;
}

std::vector<std::complex<double>> getRoots(double kappa, const std::function<double(double, double)>& f) {
	std::vector<std::complex<double>> result;
	auto ff = [kappa, f](double alpha) {return f(alpha, kappa);};
	auto r_roots = find_roots(0, 1.5 * kappa, ff, 20);
	for (int i = 0;i < r_roots.size();i++) {
		std::complex<double> re = { r_roots[i], 0 };
		result.push_back(re);
	}
	return result;
}*/


void plotTheWaveField(const std::map<double, std::complex<double>>& waveField, const std::string& fileName)
{
	std::ofstream os(fileName);
	os << "\\begin{tikzpicture}[scale=1.5]\n";
	os << "\\begin{axis}[grid]\n";
	os << "\\addplot[smooth, red] plot coordinates{\n";
	for (auto x : waveField) {
		os << "(" << x.first << ", " << x.second.real() << ") ";
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
/*
std::vector<std::complex<double>> reconstructed(const std::vector<double>& parameters,
	const std::vector<double>& points, double kappa)
{
	Parameters::smooth_params[0] = [=](auto x) {return 1 + parameters[0] * (1 - x) + parameters[1] * x; };;
	Parameters::smooth_params[1] = [=](auto x) {return 1; };
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
	const auto reses = residualSet(roots, kappa);
	for (auto point : points)
	{
		result.push_back(waves(point, kappa, roots, reses));
	}
	return result;
}

double fitness(const std::vector<double>& parameters,
	const std::vector<double>& points,
	const std::vector<std::complex<double>>& observ, double kappa)
{
	double res = 0;
	if (points.size() != observ.size())
	{
		throw std::invalid_argument("������������� �����");
	}
	auto res1 = reconstructed(parameters, points, kappa);
	for (size_t i = 0; i < points.size(); i++) {
		res += abs(res1[i] - observ[i]) * abs(res1[i] - observ[i]);// +res;
	}
	return res;

}*/

template <typename T>
class MyObjective
{
public:
	static double kappa;
	static std::vector<T> Objective(const std::vector<T>& x)
	{
		std::vector<double> xx;// = rand_vec(64, 0,6.18);
		for (size_t i = 0; i < 10; ++i)
		{
			xx.push_back(0.3 * i);
		}

		const auto observed_data = observed(xx, kappa);
		auto f = [=](const std::vector<double>& x) {return fitness(x, xx, observed_data, kappa); };
		T obj = -f(x);///(pow(1 - x[0], 2) + 100 * pow(x[1] - x[0] * x[0], 2));
		return { obj };
	}
};


// constraints example:
// 1) x * y + x - y + 1.5 <= 0
// 2) 10 - x * y <= 0
template <typename T>
std::vector<T> MyConstraint(const std::vector<T>& x)
{
	return { -x[0], -x[1],x[0] - 1,x[1] - 1 };
}

/*
void set_parameters(const std::vector<double>& parameters) {
	Parameters::smooth_params.clear();
	Parameters::smooth_params.push_back([=](double x) {return 1 + parameters[0] * x + parameters[1] * (1 - x);});
	Parameters::smooth_params.push_back([=](double x) {return 1;});
}
*/

double MyObjective<double>::kappa = 1.0;

/*
// перегрузка observed
std::vector<std::complex<double>> observed(const std::vector<double>& points, double kappa,
	const std::function<double(double)>& mu, const std::function<double(double)>& rho) {
	Parameters::smooth_params[0] = mu;//[=](auto x) {return 1 + x; };
	Parameters::smooth_params[1] = rho;//[=](auto x) {return 1; };
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
	const auto reses = residualSet(roots, kappa);
	for (auto point : points)
	{
		result.push_back(waves(point, kappa, roots, reses));
	}
	return result;
}*/

//========================================================================================================================================================
//========================================================================================================================================================

int main()
{
	setlocale(0, "");
	/*
	//initiliazing genetic algorithm
	galgo::GeneticAlgorithm<double> ga(MyObjective<double>::Objective, 100, LB, UB, 50, true);

	//setting constraints
	ga.Constraint = MyConstraint;

	//running genetic algorithm
	ga.run();

	*/

	//double kappa = 10.0;
	double rho0 = 1;
	const size_t rows = 5;
	double columns = 5;
	double c = 1;
	double d = 2;

	std::vector<double> points;
	for (int k = 0; k < rows;k++) {
		points.push_back(c + k * (d - c) / rows);
	}


	Parameters::smooth_params.push_back([](double x) {return 1 + x; });
	Parameters::smooth_params.push_back([](double x) {return 1; });

	double kappa = 5;
	layer l = {};
	auto roots = l.getRoots(kappa);
	auto residual_set = l.residualSet(roots, kappa);
	auto field = l.observed(points, kappa, residual_set, roots);
	auto mat = l.MatrixV(kappa, c, d, rows, columns);
	//std::vector<std::complex<double>> Vector1 = layer//observed(points, kappa, [=](auto x) {return 1 + x;}, [=](auto x) {return 1;});
	//std::vector<std::vector<std::complex<double>>> Matrix1 = MatrixV(, kappa, columns, points, rows);

	system("pause");

}
