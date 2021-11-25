#include "BoundaryValueProblemSmooth.h"
#include"Parameters.h"
#include <cassert>
complex<double> I(0, 1);
const double Pi = 3.1415926535897932384626433832795;

/**
 * \brief Расчет правой части
 * \param t поперечная координата
 * \param vec вектор неизвестных
 * \return значения правой части
 */
vector<complex<double>> BoundaryValueProblemSmooth::Evaluate(double t,
	const vector<complex<double>> & vec) const {
	assert(vec.size() % 2 == 0);
	vector<complex<double>> Result(vec.size());
	for (size_t i = 0; i < vec.size(); i++)
	{
		Result[i] = EquationsSystem[i](kappa, t, vec);
	}
	return Result;
}

/**
 * \brief Расчёт по вложенной формуле Рунге
 * \param x поперечная координата
 * \param step длина отрезка
 * \param InitialValue значения неизвестных функций на левом конце отрезка
 * \return значения неизвестных функций на правом конце отрезка
 */
double BoundaryValueProblemSmooth::EmbeddedFormula(double x, double step,
	vector<complex<double>> & InitialValue) const {
	auto y = Evaluate(x, InitialValue);
	vector<vector<complex<double>>> k;
	vector<complex<double>> kk(y.size());
	for (size_t i = 0; i < y.size(); i++) kk[i] = step * y[i];
	k.push_back(kk);
	for (int i = 0; i < 5; i++) {
		auto yh = InitialValue;
		for (int j = 0; j <= i; j++)
			for (size_t l = 0; l < yh.size(); l++)
				yh[l] += ButcherTableau[i][j + 1] * k[j][l];
		auto y2 = Evaluate(x + ButcherTableau[i][0] * step, yh);
		for (size_t i = 0; i < y.size(); i++) kk[i] = step * y2[i];
		k.push_back(kk);
	}
	vector<complex<double>> y1(y.size());
	vector<complex<double>> y2(y.size());
	for (size_t i = 0; i < y.size(); i++) {
		y1[i] = 0;
		y2[i] = 0;
		for (int j = 0; j <= 5; j++) {
			y1[i] += ButcherTableau[5][j + 1] * k[j][i];
			y2[i] += ButcherTableau[6][j + 1] * k[j][i];
		}
	}
	for (size_t i = 0; i < y1.size(); i++) y2[i] -= y1[i];
	InitialValue = y1;
	return Norma(y2) / Norma(y1);
}

/**
 * \brief Численное решение задачи Коши
 * \param a левый конец отрезка
 * \param b правый конец отрезка
 * \param InitialValue начальные условия
 * \return решение на правом конце отрезка
 */
vector<complex<double>> BoundaryValueProblemSmooth::CauchyProblem(double a,
	double b, vector<complex<double>> & InitialValue) const {
	double x = a;
	double h = b - a;
	while (abs(b - x) > eps) {
		auto yp = InitialValue;
		double n = EmbeddedFormula(x, h, yp);
		while (n > eps) {
			yp = InitialValue;
			h = 0.5*h;
			n = EmbeddedFormula(x, h, yp);
		}
		x += h;
		for (size_t i = 0; i < yp.size(); i++) InitialValue[i] += yp[i];
	}
	return InitialValue;
}


vector<complex<double>> BoundaryValueProblemSmooth::Shoot(double a, double b, vector<complex<double>> & InitialValue) {
	return vector<complex<double>>();
}



double Norma(const vector<complex<double>> & v) {
	double sum = 0;
	for (auto i : v) sum += i.real()* i.real() + i.imag()* i.imag();
	return sqrt(sum);
}

void BoundaryValueProblemSmooth::InitializeTheSystem() {
	D = [](double t) {return 1.0; };
	nu = 0.29;
	EquationsSystem = {
		[this](complex<double> kappa, double t, const vector<complex<double>> & vec) { return -vec[1]; },
		[this](complex<double> kappa, double t, const vector<complex<double>> & vec) { return vec[2] / D(t) - nu / t *vec[1]; },
		[this](complex<double> kappa, double t, const vector<complex<double>> & vec) { return (1 - nu*nu)*D(t) / t / t*vec[1] - (1 - nu)*vec[2] / t + vec[3]; },
		[this](complex<double> kappa, double t, const vector<complex<double>> & vec) { return vec[0] * kappa*kappa - vec[3] / t - q(t); }
	};
}


BoundaryValueProblemSmooth::BoundaryValueProblemSmooth() = default;


BoundaryValueProblemSmooth::~BoundaryValueProblemSmooth() = default;


BoundaryValueProblemSmooth::BoundaryValueProblemSmooth(function<double(double)> q, double innerRadius,
	double kappa, double eps) :
	BoundaryValueProblem(innerRadius, kappa, eps), q(q)
{
	InitializeTheSystem();
}

void BoundaryValueProblemSmooth::CreateTheInnerSolution(size_t points, size_t size)
{
	size_t n = points;//Parameters::rho.size();
	double h = (1.0 - innerRadius) / n;
	vector<complex<double>> initialValue(size, { 0,0 });
	initialValue[0] = { 1,0 };
	InnerCauchyProblemSolution.clear();
	double a = innerRadius;
	CauchyProblem(a, a + h / 2, initialValue);
	InnerCauchyProblemSolution.push_back(initialValue);
	a += h / 2;
	for (size_t i = 0; i < n - 1; i++) {
		CauchyProblem(a, a + h, initialValue);
		InnerCauchyProblemSolution.push_back(initialValue);
		a += h;
	}
	CauchyProblem(1 - h / 2, 1, initialValue);
	InnerCauchyProblemSolution.push_back(initialValue);
}

//                                               
vector<complex<double>> BoundaryValueProblemSmooth::CauchyProblemSolutions(double y, size_t size) const
{
	vector<complex<double>> initialValue(size, complex<double>(0, 0));
	return CauchyProblem(innerRadius, y, initialValue);
}

complex<double> BoundaryValueProblemSmooth::waveField() const
{
	vector<complex<double>> y(2);
	y[1] = { 0,0 };
	y[0] = /*{ exp(-abs(alpha)),0 };//*/{ 1,0 };
	auto cauchyProblemSolution = CauchyProblem(innerRadius, 1, y);
	return cauchyProblemSolution[0] / cauchyProblemSolution[1];
}

