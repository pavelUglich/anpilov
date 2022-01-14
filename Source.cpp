﻿#include "Layer.h"
#include "gen/Galgo.hpp"
#include "VoyevodinMethod.h"
#include<iostream>

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
		std::transform(value.begin(), value.end(), value.begin(), [](auto x) {return -x; });
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

void RecFile(vector<double> x, vector<double> y, string fileName) {
	std::ofstream os(fileName);
	os << "\\begin{tikzpicture}[scale=1.5]\n";
	os << "\\begin{axis}[grid]\n";
	os << "\\addplot[smooth, red] plot coordinates{\n";
	for (size_t i = 0;i < x.size();i++) {
		os << "(" << x[i] << ", " << y[i] << ") ";
	}
	os << "};\n";
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
		throw std::invalid_argument("                   ");
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

//========================================================================================================================================================
//========================================================================================================================================================
template<class T>
void display(const std::vector<T>& vector) {
	for (auto t : vector) {
		cout << t << ", ";
	}
	cout << '\n';
}


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

	/*
		Научиться восстанавливать функцию при непостоянном начальном 
		приближении. Когда появится результат, попытаться сделать следующее:
		1. Строим начальное приближение в виде линейной функции с помощью геретического алгоритма.
		2. Восстанавливаем решение при помощи итерационного процесса
	
	*/


	double kappa = 2.0;
	const size_t rows = 30;
	double columns = 30;
	double c = 0;
	double d = 1;
	std::vector<double> points;
	for (int k = 0; k < rows; k++) {
		points.push_back(c + k * (d - c) / rows);
	}
	
	std::vector<double> rho0(rows, 1.0);

	Parameters::smooth_params.emplace_back([](double x) {return 1 + 0.5*x*x; });
	Parameters::smooth_params.emplace_back([](double x) {return 1; });//!!

	layer l(kappa);
	auto observed_field = l.observed(points);
	Parameters::smooth_params[0] = [](double x) {return 1 + 0.5*x; };
	Parameters::smooth_params[1] = [](double x) {return 1; };
	layer l1(kappa);
	auto field2 = l1.observed(points);
	auto mat = l1.MatrixMu(columns, rows);
	std::vector<double> right_part = observed_field - field2;
	VoyevodinMethod V = { mat, right_part, 1.0 / right_part.size(), Dirichle, Dirichle };
	auto rho1 = V.solution();
	for (size_t i = 0; i < rows; i++)
	{
		rho0.push_back(1.0);
		rho1[i] += rho0[i];
	}
	int it = 0;
	while (norm(right_part) > 0.001 and it < 20)
	{
		// 1. пересчитываем rho1
		// 2. создаём новый слой
		// 3. находим поле и пересчитываем правую часть
		// 4. пересчитывем матрицу
		// 5. новая поправка методом воеводина
		layer l2(kappa, points, rho1, rho0);
		auto field3 = l2.observed(points);
		mat = l2.MatrixMu(columns, rows);
		right_part = observed_field - field3;
		VoyevodinMethod V = { mat, right_part, 1.0 / right_part.size(), Dirichle, Dirichle };
		for (size_t i = 0; i < rows; i++)
		{
			//rho0.push_back(1.0);
			rho1[i] += V.solution()[i];
		}
		it += 1;
	}
	/*
		Функции:
		1 + a*x; // a=0.1, 0.2, 0.5, 1.0 (дирихле, нейман)
		1+a*x*x;
		1+a*sin(pi*x) ; //(дирихле, дирихле)
		1+a*cos(pi*x); // -0.5<=a<=0.5 (нейман, нейман)
		exp(-x)

		1, 2, 5, 10

	*/
	plotTheWaveField({ {"black", rho1} }, "rho.txt", (d - c) / rows);
	cout << it;
	system("pause");
	return 0;
}