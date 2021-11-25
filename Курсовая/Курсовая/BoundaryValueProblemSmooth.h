#pragma once
#include "BoundaryValueProblem.h"
#include<complex>
#include<functional>
#include"Parameters.h"





/**
 * \brief Таблица Бутчера
 */
const vector<vector<double>> ButcherTableau = {
	{ 0.25,0.25 },
{ 0.375, 0.09375, 0.28125 },
{ 12 / 13.0, 1932 / 2197.0, -7200 / 2197.0, 7296 / 2197.0 },
{ 1, 439 / 216.0, -8, 3680 / 513.0, -845 / 4104.0 },
{ 0.5, -8 / 27.0, 2, -3544 / 2565.0, 1859 / 4104.0, -11 / 40.0 },
{ 0, 16 / 135.0, 0, 6656 / 12825.0, 28561 / 56430.0, -9 / 50.0, 2 / 55.0 },
{ 0, 25 / 216.0, 0, 1408 / 2565.0, 2197 / 4104.0, -0.2, 0 }
};

double Norma(const vector <complex<double>> & v);


/**
 * \brief Краевая задача для однородного упругого слоя
 * (решение через приближенное решение краевой задачи)
 */
class BoundaryValueProblemSmooth :
	public BoundaryValueProblem {
	/**
	 * \brief правые части системы дифференциальных уравнений
	 */
	vector<function<complex<double>(complex<double>, double, const vector<complex<double>>&)>>
		EquationsSystem;

	// коэффициент Пуассона
	double nu;
	// цилиндрическая жёсткость
	function<double(double)> D;
	// внешняя нагрузка
	function<double(double)> q;

	
	/**
	 * \brief Рассчитать правые части дифференциальных уравнений 
	 * \return вектор правых частей
	 */
	vector<complex<double>> Evaluate(double, const vector<complex<double>> &) const;


	/**
	 * \brief Расчёт по вложенным формулам Рунге
	 * \param x поперечна координата
	 * \param step шаг формул Рунге
	 * \param InitialValue начальные уловия
	 * \return результат расчёта
	 */
	double EmbeddedFormula(double x, double step, vector<complex<double>> & InitialValue) const;
protected:

	/**
	 * \brief модуль сдвига в зависимости от координаты
	 */
	//function<complex<double>(double)> Mu;


	/**
	 * \brief плотность в зависмости от координаты
	 */
	//function<double(double)> Rho;


	/**
	 * \brief построение системы
	 */
	void InitializeTheSystem();
	
public:

	vector<complex<double>> Shoot(double a, double b, vector<complex<double>> & InitialValue);

	/**
	 * \brief Решение вспомогательной задачи Коши
	 * \return  Решение вспомогательной задачи Коши
	 */
	virtual vector<complex<double>> CauchyProblem(double, double, vector<complex<double>>&) const;

	vector<complex<double>> CauchyProblemSolutions(double y, size_t size = 2) const override;
	complex<double> waveField() const override;


	BoundaryValueProblemSmooth();
	~BoundaryValueProblemSmooth();
	BoundaryValueProblemSmooth(function<double(double)> q, double innerRadius,
		double kappa, double eps = 0.1e-5);
	/**
 * \brief решение, зависящее от поперечной координаты
 * \param points количество точек
 * \param size количество компонент решения
 */
	void CreateTheInnerSolution(size_t points, size_t size = 2) override;

};

