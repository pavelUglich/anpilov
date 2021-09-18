#include "MatrixSystem.h"
#include<vector>
#include<cassert>
#include "vector_ops.h"

double innerprod(const vector<double>& a, const vector<double>& b)
{
	assert(a.size() == b.size());
	double sum = 0;
	for (size_t i = 0; i < a.size(); i++) sum += a[i] * b[i];
	return sum;
}

/**
 * \brief конструктор
 * \param A матрица
 * \param b правая часть
 * \param step шаг
 * \param p параметр стабилизатора
 * \param left граничное условие слева
 * \param right граничное условие справа
 */
matrix_system::matrix_system(vector<vector<double>> A, const vector<double>& b,
	double step, double p, BoundaryCondition left, BoundaryCondition right) :
	Matrix(std::move(A)), right_part(b)
{
	size = Matrix.front().size();
	stabilizer = { size, step, p, left, right };
	multiply_ASinv();
	multiply_transpose_au();
	qpr(); 
	multiply_rx();
}

/**
 * \brief диагональ (геттер)
 * \return диагональ
 */
vector<double> matrix_system::diagonal() const
{
	return p1;
}

/**
 * \brief наддиагональ
 * \return наддиагональ
 */
vector<double> matrix_system::up_diagonal() const
{
	return p2;
}

/**
 * \brief правая часть
 * \return правая часть
 */
vector<double> matrix_system::rightPart() const
{
	return right_part;
}

/**
 * \brief умножение вектора на матрицу, обратную к Q
 * \param v вектор
 * \return новое значение вектора
 */
vector<double> matrix_system::multiply_qtu(const vector<double>& v)
{
	assert(size == v.size());
	auto qtu = v;
	for (size_t i = 0; i < size; i++)
	{
		vector<double> a(size - i);
		for (size_t j = 0; j < size - i; j++) a[j] = Matrix[j + i][i];
		double sc = 0;
		for (size_t k = 0; k < size - i; k++) sc += a[k] * qtu[k + i];
		for (size_t j = i; j < size; j++) qtu[j] -= 2 * a[j - i] * sc;
	}
	return qtu;
}

/**
 * \brief Умножение справа на матрицу, обратную к матрице стабилизатора 
 */
void matrix_system::multiply_ASinv()
{
	auto diagonal = stabilizer.diagonal();
	auto up_diagonal = stabilizer.Updiagonal();
	for (auto& i : Matrix) i[0] /= diagonal[0];
	for (size_t i = 1; i < size; i++)
		for (auto& j : Matrix)
		{
			j[i] -= up_diagonal[i - 1] * j[i - 1];
			j[i] /= diagonal[i];
		}
}

/**
 * \brief обнуление поддиагональных элементов столбца путем умножения слева на 
 * матрицу отражения
 * \param k номер столбца
 * \return значение диагонального элемента матрицы
 */
double matrix_system::del_col(size_t k)
{
	return delete_the_column(this->Matrix, k);
}

/**
 * \brief Умножение справа на матрицу отражения
 * \param k номер строки
 * \return новое значение элемента наддиагонали
 */
double matrix_system::del_row(size_t k)
{
	const auto l = size - k - 1;
	vector<double> av(l);
	for (size_t i = 0; i < l; i++) av[i] = Matrix[k][i + k + 1];
	av[0] -= norm(av);
	normalize(av);
	vector<double> vv(l);
	for (size_t i = 0; i < l; i++) { vv[i] = Matrix[k][i + k + 1]; }
	double sc = innerprod(vv, av);
	const double pp = Matrix[k][k + 1] - 2 * av[0] * sc;
	for (size_t i = k + 1; i < size; i++)
	{
		for (size_t j = 0; j < l; j++) vv[j] = Matrix[i][j + k + 1];
		sc = vv * av;
		for (size_t j = k + 1; j < size; j++)
			Matrix[i][j] -= 2 * av[j - k - 1] * sc;
	}
	for (size_t i = 0; i < l; i++) Matrix[k][i + k + 1] = av[i];
	return pp;
}

/**
 * \brief Преобразование правой части, умножение на транспонированную матрицу
 * СЛАУ;
 */
void matrix_system::multiply_transpose_au()
{
	vector<double> v(size);
	for (size_t i = 0; i < size; i++)
	{
		v[i] = 0;
		for (size_t j = 0; j < size; j++)
			v[i] += Matrix[j][i] * right_part[j];
	}
	right_part = v;
}

/**
 * \brief Сведение матрицы СЛАУ к двухдиагональному виду
 */
void matrix_system::qpr()
{
	p1.resize(size);
	p2.resize(size);
	for (size_t i = 0; i < size - 2; i++)
	{
		p1[i] = del_col(i);
		p2[i] = del_row(i);
	}
	p1[size - 2] = del_col(size - 2);
	p2[size - 2] = Matrix[size - 2][size - 1];
	p1[size - 1] = del_col(size - 1);
	p2[size - 1] = 0;
}

/**
 * \brief Умножение вектора правой части на матрицу, обратную к ортогональной
 * матрице R;
 */
void matrix_system::multiply_rx()
{
	for (size_t i = 0; i < size - 1; i++)
	{
		vector<double> av(size);
		for (size_t j = i + 1; j < size; j++)
			av[j] = Matrix[i][j];
		double sc = 0;
		for (size_t j = i + 1; j < size; j++)
			sc += av[j] * right_part[j];
		for (size_t j = i + 1; j < size; j++) right_part[j] -= 2 * av[j] * sc;
	}
}

/**
 * \brief Умножение вектора u на матрицу, обратную к ортогональной матрице R;
 * \param u вектор правой части
 */
void matrix_system::multiply_rtx(vector<double>& u)
{
	auto v = u;
	for (size_t i = 0; i < size; i++)
	{
		const size_t l = size - i;
		vector<double> a(i);
		for (size_t j = 0; j < i; j++)
			a[j] = Matrix[l - 1][j + l];
		double sc = 0;
		for (size_t j = 0; j < i; j++) sc += a[j] * v[j + l];
		for (auto j = l; j < size; j++) v[j] -= 2 * a[j - l] * sc;
	}
	u = v;
}

/**
 * \brief Умножение вектора u на матрицу, обратную к матрице стабилизатора;
 * \param u вектор правой части
 */
void matrix_system::multiply_Sinv(vector<double>& u) const
{
	auto diagonal = stabilizer.diagonal();
	auto up_diagonal = stabilizer.Updiagonal();
	auto x = u;
	x[size - 1] = u[size - 1] / diagonal[size - 1];
	for (size_t i = 1; i < size; i++)
	{
		const auto j = size - i - 1;
		x[j] = (u[j] - up_diagonal[j] * x[j + 1]) / diagonal[j];
	}
	u = x;
}

double norm(const vector<double>& v)
{
	double result = 0;
	for (auto x : v) result += x * x;
		return sqrt(result);
}
