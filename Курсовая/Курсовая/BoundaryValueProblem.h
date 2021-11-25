#include<vector>
#include<complex>


using namespace std;

#pragma once


/**
 * \brief Краевая задача для однородного упругого слоя 
 * (решение выражается аналитически)
 */
class BoundaryValueProblem {
protected:
	
	/**
	 * \brief Внутренний радиус цилиндра
	 */
	double innerRadius;

	/**
	 * \brief Решение вспомогательной задачи Коши
	 */
	vector<complex<double>> CauchyProblemSolution;

	/**
	 * \brief Погрешность расчетов
	 */
	double eps;

	/**
	 * \brief параметр преобразования Фурье
	 */
	//complex<double> alpha;

	/**
	 * \brief Волновое число
	 */
	double kappa;


	
public:

	//vector<complex<double>> Shoot(double a, double b, vector<complex<double>> & InitialValue) const = 0;

	virtual vector<complex<double>> CauchyProblem(double, double, vector<complex<double>>&) const = 0;

	/**
	 * \brief Решение вспомогательной задачи Коши
	 * \param y 
	 * \param size 
	 * \return 
	 */
	virtual vector<complex<double>> CauchyProblemSolutions(double y,
		size_t size = 4) const = 0;

	vector<vector<complex<double>>> InnerCauchyProblemSolution;	
	BoundaryValueProblem();

	/**
	 * \brief Конструктор 
	 * \param innerRadius внутренний радиус цилиндра
	 * \param kappa волновое число
	 * \param alpha параметр преобразования Фурье
	 * \param eps погрешность вычислений
	 */
	BoundaryValueProblem(double innerRadius, double kappa, double eps = 0.1e-5);
	/**
	 * \brief деструктор
	 */
	virtual ~BoundaryValueProblem();


	/**
	 * \brief Метод Ньютона
	 * \return корень
	 */
	complex<double> NewtonMethod();


	/**
	 * \brief 
	 * \return 
	 */
	vector<complex<double>> GetCauchyProblemSolution() const { return CauchyProblemSolution; };
	
	/**
	 * \brief Знаменатель подынтегрального выражения
	 * \param size размер системы дифференциальных уравнений
	 * \return знаменатель
	 */
	complex<double> IntegrandDenumerator(size_t size = 2) const;


	/**
	 * \brief Волновое поле на внешней поверхности цилиндра
	 * \return Волновое поле на внешней поверхности цилиндра
	 */
	virtual complex<double> waveField() const;

	/**
	 * \brief изменить значение волнового числа.
	 * \param _alpha 
	 */
	void SetAlpha(complex<double> _alpha);
	
	/**
	 * \brief решение, зависящее от поперечной координаты
	 * \param points количество точек
	 * \param size количество компонент решения
	 */
	virtual void CreateTheInnerSolution(size_t points, size_t size = 4) = 0;

	/**
	 * \brief геттер
	 * \return волновое число
	 */
	double getKappa() const;
	double getq() const;
};

