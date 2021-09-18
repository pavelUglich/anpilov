#include "Layer.h"

//std::vector<std::function<std::complex<double>(double, const std::vector<std::complex<double>>&)>> 

std::vector<std::function<std::complex<double>(double, const std::vector<std::complex<double>>&)>> layer::get_equation(std::complex<double> alpha, double kappa, size_t size) const
{
	std::vector<std::function<std::complex<double>(double, const std::vector<std::complex<double>>&)>> result(size);
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i] = [=](double x, const std::vector<std::complex<double>>& v) {return this->_equations[i](x, v, alpha, kappa);};
	}
	return result;
}

double layer::defOnTop(std::complex<double> alpha, double kappa) const
{

	OdeSolver<std::complex<double>> cauchy_problem = { this->get_equation(alpha, kappa, 2)	, 0.1e-6, RUNGE_KUTTA_FELDBERG };
	auto solution = cauchy_problem.solve(0, 1, { 0, 1 });
	return solution[1].real();
}

double layer::findRoot(double a, double b, const std::function<double(double)>& f, double epsilon) // !!
{
	while (abs(b - a) > epsilon) {
		const double fa = f(a);
		const double fb = f(b);
		a = b - (b - a) * fb / (fb - fa);
		b = a - (a - b) * fa / (fa - fb);
	}
	return b;
}

std::vector<double> layer::find_roots(double a, double b, const std::function<double(double)>& f, size_t n, double epsilon)
{
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

std::map<double, std::vector<double>> layer::dispersionalSet(double max_kappa, double step, size_t n, const std::function<double(double, double)>& f)
{
	std::map<double, std::vector<double>> result;
	double kappa = step;
	while (kappa < max_kappa) {
		auto ff = [kappa, f](double alpha) {return f(alpha, kappa);};
		result[kappa] = find_roots(0, max_kappa, ff, n);
		kappa += step;
	}
	return result;
}

std::complex<double> layer::waves(double x1, double kappa, const std::vector<std::complex<double>>& roots, const std::vector<std::complex<double>>& residuals)
{
	std::complex<double> result = 0;
	std::complex<double> im = { 0,1 };
	for (int i = 0; i < roots.size(); i++) {
		result += im * residuals[i] * exp(im * roots[i] * x1);
	}
	return result;
}

std::map<double, std::complex<double>> layer::waveField(double a, double b, double kappa, double step, const std::vector<std::complex<double>>& roots, const std::vector<std::complex<double>>& residuals)
{
	std::map<double, std::complex<double>> result;
	for (double i = a; i < b; i += step) {
		result[i] = waves(i, kappa, roots, residuals);
	}
	return result;
}

std::complex<double> layer::residual(std::complex<double> alpha, double kappa)
{
	OdeSolver<std::complex<double>> cauchy_problem = { this->get_equation(alpha, kappa, 4)	, 0.1e-6, RUNGE_KUTTA_FELDBERG };
	auto solution = cauchy_problem.solve(0, 1, { 0, 1, 0, 0 });
	return solution[0] / solution[3];
}

std::vector<std::complex<double>> layer::residualSet(const std::vector<std::complex<double>>& roots, double kappa)
{
	std::vector<std::complex<double>> result;
	for (int i = 0;i < roots.size();i++) {
		result.push_back(residual(roots[i], kappa));
	}
	return result;
}

/*std::vector<std::complex<double>> layer::addKer(
	//std::complex<double> alpha, double kappa, double n, double x1
	double kappa, double b, double c, size_t rows, size_t columns
)
{
	std::vector<std::complex<double>> result;
	// 3. Делаем расчет по формуле (17)
	/*
	std::vector<double> vec;
	for (int i = 0; i < n;i++) {
		vec.push_back((i + 0.5) / n);
	}
	std::complex<double> im = { 0,1 };
	OdeSolver<std::complex<double>> cauchy_problem = { this->get_equation(alpha, kappa, 6), 0.1e-6, RUNGE_KUTTA_FELDBERG };
	auto solution = cauchy_problem.solve(vec, { 0, 1, 0, 0, 0, 0 });
	std::vector<std::complex<double>> result;
	for (auto x : solution) {
		auto item = (2.0 * x[0] * (x[2] - x[0] * x[5] / x[3]) - im * x[0] * x[0] * x1) / x[3] / x[3];
		result.push_back(item);
	}

	return result;
}*/

std::vector<std::vector<std::complex<double>>> layer::MatrixV(double kappa, double c, double d, size_t rows, size_t columns)
{
	// 1. Определяем корни
	// 2. Определяем вычеты в двухкратных полюсах
	auto roots = getRoots(kappa);
	std::vector<std::vector<std::complex<double>>> mat(roots.size());
	std::complex<double> im = { 0,1 };
	std::vector<double> vec;
	std::vector<std::complex<double>> residuals;
	for (int i = 0; i < rows;i++) {
		vec.push_back((i + 0.5) / rows);
	}

	for (size_t i = 0; i < roots.size(); i++) {
		OdeSolver<std::complex<double>> cauchy_problem = { this->get_equation(roots[i], kappa, 6), 0.1e-6, RUNGE_KUTTA_FELDBERG };
		auto solution = cauchy_problem.solve(vec, { 0, 1, 0, 0, 0, 0 });
		for (size_t ii = 0; ii < solution.size(); ii++)
		{
			auto item = (2.0 * solution[ii][0] * (solution[ii][2] - solution[ii][0] * solution[ii][5] / solution[ii][3]) - im * solution[ii][0] * solution[ii][0] * vec[ii]) / solution[ii][3] / solution[ii][3];
			mat[i].push_back(item);
		}
	}

	std::vector<double> vec1;
	for (size_t i = 0; i < columns; i++)
	{
		vec1.push_back(c + (d - c) / columns * i);
	}

	std::vector<std::vector<std::complex<double>>> result(columns, std::vector<std::complex<double>>(rows, { 0,0 }));
	for (size_t i = 0; i < columns; i++)
	{
		for (size_t j = 0; j < rows; j++) {
			for (size_t k = 0; k < roots.size(); k++)
			{
				result[i][j] += im * exp(im * roots[k] * vec[i]) * mat[k][j];
			}
		}
	}
	return result;
}

/*std::vector<std::complex<double>> layer::Kernels(std::vector<std::vector<std::complex<double>>> mat)
{
	std::vector<std::complex<double>> result;
	for (size_t i = 0; i < mat.size(); i++)
	{
		for (size_t j = 0; j < mat[i].size; i++) {
			result[i] += mat[i][j];
		}
	}
	return result;
}*/

std::vector<std::complex<double>> layer::getRoots(double kappa)
{
	std::vector<std::complex<double>> result;
	auto fun = [=](double alpha) {
		return defOnTop({ alpha,0 }, kappa);
	};
	auto fun_imag = [=](double alpha) {
		return defOnTop({ 0,alpha }, kappa);
	};
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

std::vector<std::complex<double>> layer::observed(const std::vector<double>& points, double kappa, std::vector<std::complex<double>> residual_set, std::vector<std::complex<double>> roots)
{
	std::vector<std::complex<double>> result;
	for (auto point : points)
	{
		result.push_back(waves(point, kappa, roots, residual_set));
	}
	return result;
}
