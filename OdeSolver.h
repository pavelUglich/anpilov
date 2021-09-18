#pragma once
//#include <complex>
#include <utility>
#include <vector>
#include <functional>
#include "ButcherTableau.h"
#include <stdexcept>
#include "vector_ops.h"

enum tableau {
	HEUN, BOGACKI_SHAMPINE, CASH_CARP, RUNGE_KUTTA_FELDBERG,
	DORMAND_PRINCE, RKF_78
};

inline std::vector<std::vector<double>> get_butcher_tableau(tableau t)
{
	switch (t) {
	case HEUN: return heun;
	case BOGACKI_SHAMPINE: return bogacki_shampine;
	case CASH_CARP: return cash_carp;
	case RUNGE_KUTTA_FELDBERG: return runge_kutta_feldberg;
	case DORMAND_PRINCE: return dormand_prince;
	case RKF_78: return Rkf78;
	default: throw std::invalid_argument("");
	}
}


template<class T>
class OdeSolver
{
	std::vector<std::function<T(double, const std::vector<T>&)>>
		_equations;
	double _epsilon = 0.1e-6;
	std::vector<std::vector<double>> _butcher_tableau = get_butcher_tableau(RUNGE_KUTTA_FELDBERG);
	std::vector<std::vector<T>> evaluate_k(
		const std::vector<T>& initial_conditions, double a, double h);
	std::pair<std::vector<T>, std::vector<T>> evaluate_uh(
		const std::vector<T>& initial_conditions, double a, double h);
	double calculate_the_residual(double a,
		const std::vector<T>& initial_conditions,
		double h, std::vector<T>& u);


public:
	OdeSolver(
		std::vector<std::function<T(double, const std::vector<T>&)>>
		functions, double epsilon, tableau t);
	std::vector<T> solve(double a, double b,
		const std::vector<T>& initial_conditions);
	std::vector<std::vector<T>> solve(const std::vector<double>& points,
		const std::vector<T>& initial_conditions);
};

template <class T>
std::vector<std::vector<T>> OdeSolver<T>::evaluate_k(
	const std::vector<T>& initial_conditions, double a, double h)
{
	std::vector<std::vector<T>> result;
	std::vector<T> vector(_equations.size());
	for (size_t i = 0; i < _equations.size(); i++)
	{
		vector[i] = h * _equations[i](a, initial_conditions);
	}
	result.push_back(vector);
	for (size_t i = 0; i < _butcher_tableau.size() - 2; i++)
	{
		const double xx = a + _butcher_tableau[i][0] * h;
		std::vector<T> uh(initial_conditions);
		for (size_t ii = 0; ii < uh.size(); ii++)
		{
			for (size_t iii = 0; iii < result.size(); iii++)
			{
				uh[ii] += _butcher_tableau[i][iii + 1] * result[iii][ii];
			}
		}
		for (size_t ii = 0; ii < _equations.size(); ii++)
		{
			vector[ii] = h * _equations[ii](xx, uh);
		}
		result.push_back(vector);
	}
	return result;
}

template<class T>
std::pair<std::vector<T>, std::vector<T>> OdeSolver<T>::evaluate_uh(
	const std::vector<T>& initial_conditions, double a, double h)
{
	std::vector<T> u(initial_conditions);
	std::vector<T> uh(initial_conditions);
	std::vector<std::vector<T>> k = evaluate_k(initial_conditions, a, h);
	const size_t size = _butcher_tableau.size() - 1;
	for (size_t i = 0; i < initial_conditions.size(); i++)
	{
		for (size_t ii = 1; ii < _butcher_tableau[size].size(); ii++)
		{
			u[i] += _butcher_tableau[size - 1][ii] * k[ii - 1][i];
			uh[i] += _butcher_tableau[size][ii] * k[ii - 1][i];
		}
	}
	return { u, uh };
}

template<class T>
double OdeSolver<T>::calculate_the_residual(double a,
	const std::vector<T>& initial_conditions, double h, std::vector<T>& u)
{
	const auto uh = evaluate_uh(initial_conditions, a, h);
	u = uh.first;
	std::vector<T> diff = uh.first - uh.second;
	return norm_(diff);
}


template <class T>
OdeSolver<T>::OdeSolver(
	std::vector<std::function<T(double, const std::vector<T>&)>> functions,
	double epsilon, tableau t) : _equations(std::move(functions))
{
	if (epsilon < 0)
	{
		throw std::invalid_argument("");
	}
	this->_epsilon = epsilon;
	this->_butcher_tableau = get_butcher_tableau(t);
}

template<class T>
std::vector<T> OdeSolver<T>::solve(double a, double b, const std::vector<T>& initial_conditions)
{
	if (a > b)
	{
		throw std::invalid_argument("");
	}
	double h = b - a;
	std::vector<T> initials(initial_conditions);
	while (abs(b - a) > DBL_EPSILON && a < b)
	{
		std::vector<T> u_;
		double r = calculate_the_residual(a, initials, h, u_);
		while (abs(r) > _epsilon)
		{
			h /= 2;
			r = calculate_the_residual(a, initials, h, u_);
		}
		a += h;
		initials = u_;
	}
	return initials;
}

template <class T>
std::vector<std::vector<T>> OdeSolver<T>::solve(
	const std::vector<double>& points,
	const std::vector<T>& initial_conditions)
{
	std::vector<std::vector<T>> result;
	result.push_back(solve(0, points[0], initial_conditions));
	for (size_t i = 1; i < points.size(); i++)
	{
		result.push_back(solve(points[i - 1], points[i], result[i - 1]));
	}
	return result;
}

