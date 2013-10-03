#include "stdafx.h"
#include "math_operands.h"

namespace smat{
	IMPEXP_SMALLMAT double norm(t_Complex val){
		return sqrt(std::norm(val));
	}
}

// instantiate all needed functions here and export them
// [-]
template IMPEXP_SMALLMAT void matrix::base::minus<t_Complex, t_Complex>
(const t_Matrix<t_Complex>& l, const t_Matrix<t_Complex>& r,
 t_Matrix<TypeDeduce<t_Complex,t_Complex>::type >& ret);

template IMPEXP_SMALLMAT void matrix::base::minus<double, double>
(const t_Matrix<double>& l, const t_Matrix<double>& r,
 t_Matrix<matrix::TypeDeduce<double,double>::type >& ret);
// [+]
template IMPEXP_SMALLMAT void matrix::base::plus<t_Complex, t_Complex>
(const t_Matrix<t_Complex>& l, const t_Matrix<t_Complex>& r,
 t_Matrix<TypeDeduce<t_Complex,t_Complex>::type >& ret);

template IMPEXP_SMALLMAT void matrix::base::plus<double, double>
(const t_Matrix<double>& l, const t_Matrix<double>& r,
 t_Matrix<matrix::TypeDeduce<double,double>::type >& ret);
// [*]

template IMPEXP_SMALLMAT void matrix::base::mul<t_Complex, t_Complex>
(const t_Complex l, const t_Matrix<t_Complex>& r,
 t_Matrix<matrix::TypeDeduce<t_Complex,t_Complex>::type >& ret);

template IMPEXP_SMALLMAT void matrix::base::mul<double, double>
(const double l, const t_Matrix<double>& r, 
 t_Matrix<TypeDeduce<double, double>::type >& ret);
//[/]

template IMPEXP_SMALLMAT void matrix::base::div<t_Complex, t_Complex>
(const t_Matrix<t_Complex>& r, const t_Complex l, 
 t_Matrix<matrix::TypeDeduce<t_Complex,t_Complex>::type >& ret);

template IMPEXP_SMALLMAT void matrix::base::div<double, double>
(const t_Matrix<double>& l, const double r,
 t_Matrix<TypeDeduce<double, double>::type >& ret);

//[A*B]
template IMPEXP_SMALLMAT void matrix::base::mat_mul<t_Complex, t_Complex>
(const t_Matrix<t_Complex>& l, const t_Matrix<t_Complex>& r,
 t_Matrix<TypeDeduce<t_Complex,t_Complex>::type >& ret);

template IMPEXP_SMALLMAT void matrix::base::mat_mul<double, double>
(const t_Matrix<double>& l, const t_Matrix<double>& r,
 t_Matrix<matrix::TypeDeduce<double,double>::type >& ret);
