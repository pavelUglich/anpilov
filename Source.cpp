#include "Layer.h"
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

double fitness(const std::vector<double>& parameters,
	const std::vector<double>& points, double kappa)
{
	double res = 0;
	Parameters::smooth_params[0] = [=](auto x) {return 1; };
	Parameters::smooth_params[1] = [=](auto x) {return 1 + 0.1 * x; };

	layer l(kappa);
	auto observ = l.observed(points);

	Parameters::smooth_params[0] = [=](auto x) {return 1; };
	Parameters::smooth_params[1] = [=](auto x) {return 1 + parameters[0] * x + parameters[1] * (1 - x); };

	layer l1(kappa);
	auto res1 = l1.observed(points);

	for (size_t i = 0; i < points.size(); i++) {
		res += abs(res1[i] - observ[i]) * abs(res1[i] - observ[i]);// +res;
	}
	return res;

}

template <typename T>
class MyObjective
{
public:
	static double kappa;
	static std::vector<T> Objective(const std::vector<T>& x)
	{
		std::vector<double> xx;// = rand_vec(64, 0,6.18);
		for (size_t i = 0; i < 30; ++i)
		{
			xx.push_back(i / 30);
		}
		auto f = [=](const std::vector<double>& x) {return fitness(x, xx, kappa); };
		T obj = f(x);///(pow(1 - x[0], 2) + 100 * pow(x[1] - x[0] * x[0], 2));
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

//========================================================================================================================================================
//========================================================================================================================================================
template<class T>
void display(const std::vector<T>& vector) {
	for (auto t : vector) {
		cout << t << ", ";
	}
	cout << '\n';
}

double MyObjective<double>::kappa = 2.0;

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

	/*
	double kappa = 1.0;
	const size_t rows = 30;
	const size_t columns = 30;
	double c = 0;
	double d = 1;
	std::vector<double> points;
	for (int k = 0; k < rows; k++) {
		points.push_back(c + k * (d - c) / rows);
	}
	
	std::vector<double> rho0(rows, 1.0);

	Parameters::smooth_params.emplace_back([](double x) {return 1+x; });
	Parameters::smooth_params.emplace_back([](double x) {return 1; });//!!

	layer l(kappa);
	auto observed_field = l.observed(points);
	Parameters::smooth_params[0] = [](double x) {return 1; };
	Parameters::smooth_params[1] = [](double x) {return 1; };
	layer l1(kappa);
	auto field2 = l1.observed(points);
	auto mat = l1.MatrixMu(columns, rows);
	std::vector<double> right_part = observed_field - field2;
	VoyevodinMethod V = { mat, right_part, 1.0 / right_part.size() };
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
		VoyevodinMethod V = { mat, right_part, 1.0 / right_part.size() };
		for (size_t i = 0; i < rows; i++)
		{
			//rho0.push_back(1.0);
			rho1[i] += V.solution()[i];
		}
		it += 1;
	}
	*/
	/*
		Функции:
		1 + a*x; // a=0.1, 0.2, 0.5, 1.0 (дирихле, нейман)
		1+a*x*x;
		1+a*sin(pi*x) ; //(дирихле, дирихле)
		1+a*cos(pi*x); // -0.5<=a<=0.5 (нейман, нейман)
		exp(-x)

		1, 2, 5, 10

	

	plotTheWaveField({ {"black", rho1} }, "rho.txt", (d - c) / rows);
	cout << it;
	*/
	double kappa = 2.0;
	const size_t rows = 30;
	const size_t columns = 30;
	double c = 0;
	double d = 1;
	std::vector<double> points;
	for (int k = 0; k < rows; k++) {
		points.push_back(k / 30.0);
	}

	std::vector<double> rho0(rows, 0.0);

	Parameters::smooth_params.emplace_back([](double x) {return 1; });
	Parameters::smooth_params.emplace_back([](double x) {return 1; });

	std::vector<double> LB(2, 0.0);
	std::vector<double> UB(2, 0.2);

	galgo::GeneticAlgorithm<double> ga(MyObjective<double>::Objective, 50, LB, UB, 10, true);

	ga.Constraint = MyConstraint;

	ga.run();

	auto chrom = ga.result()->getParam();

	Parameters::smooth_params[0] = [=](auto x) {return 1; };
	Parameters::smooth_params[1] = [=](auto x) {return 1 + 0.1 * x; };

	layer ll(kappa);
	auto observed_field = ll.observed(points);

	Parameters::smooth_params[0] = [=](auto x) {return 1; };
	Parameters::smooth_params[1] = [=](auto x) {return 1 + chrom[0] * x + chrom[1] * (1 - x); };
	
	layer lll(kappa);
	auto rec_field = lll.observed(points);
	auto mat = lll.MatrixRho(columns, rows);

	auto right_part = observed_field - rec_field;

	VoyevodinMethod V = { mat, right_part, 1.0 / right_part.size(), Dirichle };
	auto rho1 = V.solution();
	for (size_t i = 0; i < rows; i++)
	{
		rho0[i] += 1 + chrom[0] * points[i] + chrom[1] * (1 - points[i]);
		rho1[i] += rho0[i];
	}
	std::vector<double> mu0(rows, 1.0);
	int it = 0;
	cout << it;
	while (norm(right_part) > 0.001 and it < 20)
	{
		layer l2(kappa, points, mu0, rho1);
		auto field3 = l2.observed(points);
		mat = l2.MatrixRho(columns, rows);
		right_part = observed_field - field3;
		VoyevodinMethod V = { mat, right_part, 1.0 / right_part.size(), Dirichle };
		for (size_t i = 0; i < rows; i++)
		{
			rho1[i] += V.solution()[i];
		}
		it += 1;
		cout << it;
	}
	plotTheWaveField({ {"black", rho1} }, "rho.txt", (d - c) / rows);
	cout << it;
	//метод Бройдена
	//Тихонов-Арсенин
	system("pause");
	return 0;
}