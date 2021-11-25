#include "Layer.h"

const std::complex<double> im = { 0,1 };

std::vector<std::function<std::complex<double>(double, const std::vector<std::complex<double>>&)>> layer::get_equation(std::complex<double> alpha, size_t size) const
{
	std::vector<std::function<std::complex<double>(double, const std::vector<std::complex<double>>&)>> result(size);
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i] = [=](double x, const std::vector<std::complex<double>>& v) {return this->_equations[i](x, v, alpha, kappa); };
	}
	return result;
}

layer::layer(double kkappa): mu(Parameters::smooth_params[0]), rho(Parameters::smooth_params[1])
{
	kappa = kkappa;
	_equations = {
		[=](double x, const std::vector<std::complex<double>>& v, std::complex<double> alpha, double kappa)
		{
			return v[1] / mu(x);
		},
		[=](double x, const std::vector<std::complex<double>>& v, std::complex<double> alpha, double kappa)
		{
			return (alpha * alpha * mu(x) - kappa * kappa * rho(x)) * v[0];
		},
		[=](double x, const std::vector<std::complex<double>>& v, std::complex<double> alpha, double kappa)
		{
			return v[3] / mu(x);
		},
		[=](double x, const std::vector<std::complex<double>>& v, std::complex<double> alpha, double kappa)
		{
			return 2.0 * alpha * v[0] * mu(x) + (alpha * alpha * mu(x) - kappa * kappa * rho(x)) * v[2];
		},
		[=](double x, const std::vector<std::complex<double>>& v, std::complex<double> alpha, double kappa)
		{
			return v[5] / mu(x);
		},
		[=](double x, const std::vector<std::complex<double>>& v, std::complex<double> alpha, double kappa)
		{
			return 2.0 * v[0] * mu(x) + 4.0 * v[2] * mu(x)*alpha + (alpha * alpha * mu(x) - kappa * kappa * rho(x)) * v[4];//!!!
		}
	};
	roots = getRoots();
	residual_set = residualSet();
}

double layer::defOnTop(std::complex<double> alpha) const
{

	OdeSolver<std::complex<double>> cauchy_problem = { this->get_equation(alpha, 2)	, 0.1e-6, RUNGE_KUTTA_FELDBERG };
	auto solution = cauchy_problem.solve(0, 1, { 0, 1 });
	return solution[1].real();
}

double layer::findRoot(double a, double b, const std::function<double(double)>& f, double epsilon) const
{
	while (abs(b - a) > epsilon) {
		const double fa = f(a);
		const double fb = f(b);
		a = b - (b - a) * fb / (fb - fa);
		b = a - (a - b) * fa / (fa - fb);
	}
	return b;
}

std::vector<double> layer::find_roots(double a, double b, const std::function<double(double)>& f, size_t n, double epsilon) const
{
	std::vector<double> result;
	const auto step = (b - a) / n;
	for (size_t i = 0; i < n; i += 1) {
		const auto left = a + i * step;
		const auto right = left + step;
		if (f(left) * f(right) < 0) {
			double f_r = findRoot(left, right, f, epsilon);
			result.push_back(f_r);
		}
	}
	return result;
}

std::map<double, std::vector<double>> layer::dispersionalSet(double max_kappa, double step, size_t n, const std::function<double(double, double)>& f) const
{
	std::map<double, std::vector<double>> result;
	double kappa = step;
	while (kappa < max_kappa) {
		auto ff = [kappa, f](double alpha) {return f(alpha, kappa); };
		result[kappa] = find_roots(0, max_kappa, ff, n);
		kappa += step;
	}
	return result;
}

std::complex<double> layer::waves(double x1, const std::vector<std::complex<double>>& roots, const std::vector<std::complex<double>>& residuals) const
{
	std::complex<double> result = 0;
	const std::complex<double> im = { 0,1 };
	for (size_t i = 0; i < roots.size(); i++) {
		result += im * residuals[i] * exp(im * roots[i] * x1);
	}
	return result;
}

std::map<double, std::complex<double>> layer::waveField(double a, double b, double step, const std::vector<std::complex<double>>& residuals) const
{
	std::map<double, std::complex<double>> result;
	for (double i = a; i < b; i += step) {
		result[i] = waves(i, roots, residuals);
	}
	return result;
}

std::complex<double> layer::residual(std::complex<double> alpha) const
{
	OdeSolver<std::complex<double>> cauchy_problem = { this->get_equation(alpha, 4)	, 0.1e-6, RUNGE_KUTTA_FELDBERG };
	auto solution = cauchy_problem.solve(0, 1, { 0, 1, 0, 0 });
	return solution[0] / solution[3];
}

std::vector<std::complex<double>> layer::residualSet()
{
	std::vector<std::complex<double>> result(roots.size());
	for (size_t i = 0; i < roots.size(); i++) {
		result[i] = residual(roots[i]);
	}
	return result;
}

std::vector<std::vector<double>> layer::MatrixRho(size_t columns, size_t rows) const
{
	const auto h = 1.0 / columns;
	std::vector<std::vector<std::complex<double>>> mat(roots.size());
	std::vector<double> vec;
	std::vector<std::complex<double>> residuals;
	for (int i = 0; i < rows; i++) {
		vec.push_back((i + 0.5) / rows);
	}
	vec.push_back(1.0);
	for (size_t i = 0; i < roots.size(); i++) {
		OdeSolver<std::complex<double>> cauchy_problem = { this->get_equation(roots[i], 6), 0.1e-6, RUNGE_KUTTA_FELDBERG }; // глобальная точность
		auto solution = cauchy_problem.solve(vec, { 0, 1, 0, 0, 0, 0 });
		auto denum = cauchy_problem.solve(0,1, { 0, 1, 0, 0, 0, 0 });
		//auto denum = solution.back(); // denum[5], denum[3]
		for (size_t ii = 0; ii < solution.size() - 1; ii++)
		{
			auto item = (2.0 * solution[ii][0] * (solution[ii][2] - solution[ii][0] * 0.5 * denum[5] / denum[3]) 
				+ im * solution[ii][0] * solution[ii][0] * vec[ii]) / denum[3] / denum[3];//!!! +
			mat[i].push_back(item);
		}
	}
	std::vector<std::vector<std::complex<double>>> result(columns, std::vector<std::complex<double>>(rows, { 0,0 }));
	for (size_t i = 0; i < columns; i++)
	{
		for (size_t j = 0; j < rows; j++) {
			for (size_t k = 0; k < roots.size(); k++)
			{
				result[i][j] += im * exp(im * roots[k] * vec[i]) * mat[k][j] * kappa * kappa;
			}
		}
	}
	std::vector<std::vector<double>> realresult(columns, std::vector<double>(rows));
	for (size_t i = 0; i < columns; i++)
	{
		for (size_t j = 0; j < rows; j++)
		{
			realresult[i][j] += result[i][j].real();
		}
	}

	return realresult;
}

std::vector<std::vector<double>> layer::MatrixMu(size_t columns, size_t rows)
{
	const auto h = 1.0 / columns;
	std::vector<std::vector<std::complex<double>>> mat(roots.size());
	std::vector<double> vec;
	std::vector<std::complex<double>> residuals;
	for (size_t i = 0; i < rows; i++) {
		vec.push_back((i + 0.5) / rows);
	}
	vec.push_back(1.0);
	for (size_t i = 0; i < roots.size(); i++) {
		OdeSolver<std::complex<double>> cauchy_problem = { this->get_equation(roots[i], 6), 0.1e-6, RUNGE_KUTTA_FELDBERG };
		auto solution = cauchy_problem.solve(vec, { 0, 1, 0, 0, 0, 0 });
		auto denum = solution.back();
		for (size_t ii = 0; ii < solution.size(); ii++)
		{
			auto item = (2.0 * solution[ii][1] * (solution[ii][3] - solution[ii][1] * 0.5 * denum[5] / denum[3]) - im * solution[ii][1] * solution[ii][1] * vec[ii]) / denum[3] / denum[3]
			+ (2.0 * roots[ii] * solution[ii][0] * (solution[ii][0] + roots[ii] * solution[ii][2] - roots[ii] * solution[ii][0] * denum[5] / denum[3]) - im * roots[ii] * solution[ii][0] * roots[ii] * solution[ii][0] * vec[ii]) / denum[3] / denum[3];
			mat[i].push_back(item.real());
		}
	}

	std::vector<std::vector<std::complex<double>>> result(columns, std::vector<std::complex<double>>(rows, { 0,0 }));
	for (size_t i = 0; i < columns; i++)
	{
		for (size_t j = 0; j < rows; j++) {
			for (size_t k = 0; k < roots.size(); k++)
			{
				result[i][j] += im * exp(im * roots[k] * vec[i]) * mat[k][j] * kappa * kappa * h;
			}	
		}
	}
	std::vector<std::vector<double>> realresult(columns, std::vector<double>(rows));
	for (size_t i = 0; i < columns; i++)
	{
		for (size_t j = 0; j < rows; j++)
		{
			realresult[i][j] += result[i][j].real();
		}
	}
	return realresult;
}

std::vector<std::complex<double>> layer::getRoots()
{
	std::vector<std::complex<double>> result;
	auto fun = [=](double alpha) {
		return defOnTop({ alpha,0 });
	};
	auto fun_imag = [=](double alpha) {
		return defOnTop({ 0,alpha });
	};
	auto real_roots = find_roots(0, 1.5 * kappa, fun, 20);
	result.reserve(real_roots.size());
	for (double& real_root : real_roots)
	{
		result.emplace_back(real_root, 0);
	}
	auto imag_roots = find_roots(0, 50, fun_imag, 20);
	for (double& imag_root : imag_roots)
	{
		result.emplace_back(0, imag_root);
	}
	return result;
}

std::vector<double> layer::observed(const std::vector<double>& points) const
{
	std::vector<double> result;
	result.reserve(points.size());
	for (auto point : points)
	{
		result.push_back(waves(point, roots, residual_set).real());
	}
	return result;
}
