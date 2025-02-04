#include "stdafx.h"
#include "small_mat.h"

// try to instantiate needed classes

template class IMPEXP_SMALLMAT t_Matrix<t_Complex>;
template class IMPEXP_SMALLMAT t_Matrix<double>;

template class IMPEXP_SMALLMAT t_Vec<t_Complex>;
template class IMPEXP_SMALLMAT t_Vec<double>;

template class IMPEXP_SMALLMAT t_Vec3<t_Complex>;
template class IMPEXP_SMALLMAT t_Vec3<double>;

template class IMPEXP_SMALLMAT t_SqMatrix<t_Complex>;
template class IMPEXP_SMALLMAT t_SqMatrix<double>;

template class IMPEXP_SMALLMAT t_SqMat3<t_Complex>;
template class IMPEXP_SMALLMAT t_SqMat3<double>;

// matrix operators 
template IMPEXP_SMALLMAT void matrix::base::minus<t_Complex, t_Complex>
(const t_Matrix<t_Complex>& l, const t_Matrix<t_Complex>& r,
	t_Matrix<TypeDeduce<t_Complex, t_Complex>::type >& ret);

template IMPEXP_SMALLMAT void matrix::base::minus<double, double>
(const t_Matrix<double>& l, const t_Matrix<double>& r,
	t_Matrix<matrix::TypeDeduce<double, double>::type >& ret);
// [+]
template IMPEXP_SMALLMAT void matrix::base::plus<t_Complex, t_Complex>
(const t_Matrix<t_Complex>& l, const t_Matrix<t_Complex>& r,
	t_Matrix<TypeDeduce<t_Complex, t_Complex>::type >& ret);

template IMPEXP_SMALLMAT void matrix::base::plus<double, double>
(const t_Matrix<double>& l, const t_Matrix<double>& r,
	t_Matrix<matrix::TypeDeduce<double, double>::type >& ret);
// [*]

template IMPEXP_SMALLMAT void matrix::base::mul<t_Complex, t_Complex>
(const t_Complex l, const t_Matrix<t_Complex>& r,
	t_Matrix<matrix::TypeDeduce<t_Complex, t_Complex>::type >& ret);

template IMPEXP_SMALLMAT void matrix::base::mul<double, double>
(const double l, const t_Matrix<double>& r,
	t_Matrix<TypeDeduce<double, double>::type >& ret);
//[/]

template IMPEXP_SMALLMAT void matrix::base::div<t_Complex, t_Complex>
(const t_Matrix<t_Complex>& r, const t_Complex l,
	t_Matrix<matrix::TypeDeduce<t_Complex, t_Complex>::type >& ret);

template IMPEXP_SMALLMAT void matrix::base::div<double, double>
(const t_Matrix<double>& l, const double r,
	t_Matrix<TypeDeduce<double, double>::type >& ret);

//[A*B]
template<> IMPEXP_SMALLMAT void matrix::base::mat_mul<t_Complex, t_Complex>
	(const t_Matrix<t_Complex>& l, const t_Matrix<t_Complex>& r,
		t_Matrix<TypeDeduce<t_Complex, t_Complex>::type >& ret) {
	__mat_mul<t_Complex, t_Complex, TypeDeduce<t_Complex, t_Complex>::type>(l, r, ret);
};

	template<> IMPEXP_SMALLMAT void matrix::base::mat_mul<double, double>
		(const t_Matrix<double>& l, const t_Matrix<double>& r,
			t_Matrix<matrix::TypeDeduce<double, double>::type >& ret) {
		__mat_mul<double, double, TypeDeduce<double, double>::type>(l, r, ret);
	};

// vector operands


template<> IMPEXP_SMALLMAT t_Vec<double> operator*(const t_Vec<double>& vec, const double num) {
	return __vecbynum_mul<double, double, matrix::TypeDeduce<double, double>::type>(vec, num);
};

template<> IMPEXP_SMALLMAT t_Vec<double> operator*(const double num, const t_Vec<double>& vec) {
	return __numbyvec_mul<double, double, matrix::TypeDeduce<double, double>::type>(num, vec);
};

template<> IMPEXP_SMALLMAT t_Vec<double> operator/(const t_Vec<double>& vec, const double num) {
	return __vecbynum_div<double, double, matrix::TypeDeduce<double, double>::type>(vec, num);
};

template<> IMPEXP_SMALLMAT t_Vec<double> operator+(const t_Vec<double>& l, const t_Vec<double>& r) {
	return __vec_plus<double, double, matrix::TypeDeduce<double, double>::type>(l, r);
};

template<> IMPEXP_SMALLMAT t_Vec<double> operator-(const t_Vec<double>& l, const t_Vec<double>& r) {
	return __vec_minus<double, double, matrix::TypeDeduce<double, double>::type>(l, r);
};

template<> IMPEXP_SMALLMAT t_Vec<double> operator*(const t_Matrix<double>& l, const t_Vec<double>& r) {
	return __matbyvec_mul<double, double, matrix::TypeDeduce<double, double>::type>(l, r);
};


template<> IMPEXP_SMALLMAT t_Vec<t_Complex> operator*(const t_Matrix<double>& l, const t_Vec<t_Complex>& r) {
	return __matbyvec_mul<double, t_Complex, matrix::TypeDeduce<double, t_Complex>::type>(l, r);
};

// 3d vector operands

template<> IMPEXP_SMALLMAT t_Vec3<double> operator*(const double num, const t_Vec3<double>& vec) {
	return __numbyvec3_mul<double, double, matrix::TypeDeduce<double, double>::type>(num, vec);
};

template<> IMPEXP_SMALLMAT  t_Vec3<double> operator+(const t_Vec3<double>& l, const t_Vec3<double>& r) {
	return __vec3_plus<double, double, matrix::TypeDeduce<double, double>::type>(l, r);
};

template<> IMPEXP_SMALLMAT t_Vec3<double> operator-(const t_Vec3<double>& l, const t_Vec3<double>& r) {
	return __vec3_minus<double, double, matrix::TypeDeduce<double, double>::type>(l, r);
};

template<> IMPEXP_SMALLMAT void t_SqMatrix<t_CompVal>::setToHermConj() {

	t_CompVal swp_val;
	for (int i = 0; i<_ncols; i++)
		for (int j = 0; j <= i; j++) {
			swp_val = _cont[i][j];
			_cont[i][j] = std::conj(_cont[j][i]);
			_cont[j][i] = std::conj(swp_val);

		}
};

template<> IMPEXP_SMALLMAT void t_SqMatrix<double>::setToHermConj() {
	wxString msg(_T("Smat error: setToHermConj not implemented \
					for specified template parameter"));
	ssuTHROW(t_GenException, msg);
};

namespace smat{
	IMPEXP_SMALLMAT double norm(t_Complex val){
		return sqrt(std::norm(val));
	}

	IMPEXP_SMALLMAT void vec_cart_to_cyl(t_Vec3<double>& vec){
		wxLogError(_T("Vec3D: Transition to cyl rf not implemented!"));
	}

	IMPEXP_SMALLMAT void vec_cart_to_cyl(t_Vec3<t_Complex>& vec){
		wxLogError(_T("Vec3D: Transition to cyl rf not implemented!"));
	}

	IMPEXP_SMALLMAT void vec_cart_to_cone(t_Vec3<double>& vec, double teta){

		const double PI = acos(-1.0);

		double x = vec[0];
		double y = vec[1];
		double z = vec[2];

		double r_yz = sqrt(y*y + z*z);

		double l = cos(teta)*x + sin(teta)*r_yz;

		double phi = 0.0;

		if (r_yz>0.0) phi = (z>0.0) ? acos(y/r_yz) : 2*PI - acos(y/r_yz);

		double h = -sin(teta)*x + cos(teta)*r_yz;

		vec.set(l, phi, h);

	}

	IMPEXP_SMALLMAT void vec_cart_to_cone(t_Vec3<t_Complex>& vec, double teta){

		t_Vec3Dbl tmp_vec;
		t_Complex val0, val1, val2;
		t_Complex iun(0.0, 1.0);

		tmp_vec.set(vec[0].real(), vec[1].real(), vec[2].real());
		vec_cart_to_cone(tmp_vec, teta);

		val0 = tmp_vec[0];
		val1 = tmp_vec[1];
		val2 = tmp_vec[2];

		tmp_vec.set(vec[0].imag(), vec[1].imag(), vec[2].imag());
		vec_cart_to_cone(tmp_vec, teta);

		val0+=tmp_vec[0]*iun;
		val1+=tmp_vec[1]*iun;
		val2+=tmp_vec[2]*iun;

		vec.set(val0, val1, val2);

	}
}

namespace vector{

    template<> IMPEXP_SMALLMAT t_Vec3Dbl cross(const t_Vec3Dbl& l, const t_Vec3Dbl& r){
	return __cross<double, double, double>(l,r);
    };

}

// interpolation functions

#include "fun_interpolate.h"

template IMPEXP_SMALLMAT double smat::interpolate_parab<double>
(double x1, double f1, double x2, double f2, double x3, double f3, double x);

template IMPEXP_SMALLMAT t_Complex smat::interpolate_parab<t_Complex>
(double x1, t_Complex f1, double x2, t_Complex f2, double x3, t_Complex f3, double x);