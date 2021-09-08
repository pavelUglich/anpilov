#pragma once
#include <vector>
#include <functional>
#include <map>

#include "OdeSolver.h"
#include "Parameters.h"

template<class T>
class boundary_value_problem
{
	std::vector<std::function<T(double, const std::vector<T>&)>>
		_equations;                            // правые части ОДУ
	std::map<size_t, T> left_conditions;  // граничные условия слева
	std::map<size_t, T> right_conditions; // граничные условия справа
	double _epsilon; // погрешность

	std::vector<std::vector<T>> initial_conditions() const;
	std::vector<std::vector<T>> cauchy_problem_solutions() const;
	std::vector<T> get_the_initial_conditions(
		const std::vector<std::vector<T>>& solutions) const;
	T determinant() const;

public:
	boundary_value_problem(
		std::vector<std::function<T(double,
			const std::vector<T>&)>> functions,
		std::map<size_t, T> left_conditions,
		std::map<size_t, T> right_conditions, double epsilon = 0.1e-6);
	std::vector<T> solve() const;
	std::vector<std::vector<T>> solve(
		const std::vector<double>& points) const;
};

template <class T>
std::vector<std::vector<T>>
boundary_value_problem<T>::initial_conditions() const
{
	const auto size = _equations.size();
	const auto cols = size - left_conditions.size();
	std::vector<std::vector<T>> result(size);
	size_t counter = 0;
	for (size_t i = 0; i < result.size(); i++)
	{
		if (left_conditions.find(i) == left_conditions.end())
		{
			result[i] = std::vector<T>(cols, 0);
			result[i][counter++] = 1;
		}
		else
		{
			result[i] = std::vector<T>(cols);
			result[i][0] = left_conditions.at(i);
		}
	}
	return result;
}

template <class T>
std::vector<std::vector<T>>
boundary_value_problem<T>::cauchy_problem_solutions() const
{
	auto conditions = this->initial_conditions();
	const auto size = conditions[0].size();
	std::vector<std::vector<T>> solutions;
	OdeSolver<T> ode_solver(_equations, this->_epsilon, RKF_78);
	for (size_t i = 0; i < size; i++)
	{
		std::vector<T> initials(_equations.size());
		for (size_t ii = 0; ii < initials.size(); ii++)
		{
			initials[ii] = conditions[ii][i];
		}
		auto solution = Parameters::kind == THIRD ?
			ode_solver.solve(Parameters::points, initials).back() :
			ode_solver.solve(0, 1, initials);
		solutions.push_back(solution);
	}
	return solutions;
}

template<class T>
std::vector<T> boundary_value_problem<T>::get_the_initial_conditions(
	const std::vector<std::vector<T>>& solutions) const
{
	std::vector<std::vector<T>> matrix;
	std::vector<T> right_part;
	for (const auto& x : right_conditions)
	{
		right_part.push_back(x.second);
		std::vector<T> row(right_conditions.size());
		for (size_t i = 0; i < row.size(); i++)
		{
			row[i] = solutions[i][x.first];
		}
		matrix.push_back(row);
	}
	const SquareMatrix<T> slau_matrix(right_conditions.size(), matrix);
	return slau_matrix.linSolve(right_part);
}

template <class T>
T boundary_value_problem<T>::determinant() const
{
	const auto solutions = cauchy_problem_solutions();
	std::vector<std::vector<T>> matrix;
	for (const auto& x : right_conditions)
	{
		std::vector<T> row(right_conditions.size());
		for (size_t i = 0; i < row.size(); i++)
		{
			row[i] = solutions[i][x.first];
		}
		matrix.push_back(row);
	}
	const SquareMatrix<T> slau_matrix(right_conditions.size(), matrix);
	return slau_matrix.det();
}

template<class T>
boundary_value_problem<T>::boundary_value_problem(
	std::vector<std::function<T(double, const std::vector<T>&)>> functions, 
	std::map<size_t, T> left_conditions, std::map<size_t, T> right_conditions, 
	double epsilon) :
	_equations(std::move(functions)),
	left_conditions(std::move(left_conditions)),
	right_conditions(std::move(right_conditions)),
	_epsilon(epsilon)
{
}

template<class T>
std::vector<T> boundary_value_problem<T>::solve() const
{
	std::vector<std::vector<T>> solutions = cauchy_problem_solutions();
	std::vector<T> initial_conditions = get_the_initial_conditions(solutions);
	std::vector<T> result(_equations.size(), 0);
	for (size_t i = 0; i < solutions.size(); i++)
	{
		result = result + initial_conditions[i] * solutions[i];
	}
	return result;
}

template<class T>
std::vector<std::vector<T>> boundary_value_problem<T>::solve(
	const std::vector<double>& points) const
{
	const auto solutions = cauchy_problem_solutions();
	const auto initial_conditions = get_the_initial_conditions(solutions);
	std::vector<T> initials;
	size_t counter = 0;
	for (size_t i = 0; i < _equations.size(); i++)
	{
		if (left_conditions.find(i) == left_conditions.end())
		{
			initials.push_back(initial_conditions[counter++]);
		}
		else
		{
			initials.push_back(left_conditions.at(i));
		}
	}
	OdeSolver<T> ode_solver(_equations, this->_epsilon, RKF_78);
	return ode_solver.solve(points, initials);
}