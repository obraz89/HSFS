///////////////////////////////////////////////////////////////////////////////
// Project:	SmallMat
// Purpose:	Small matrix linear algebra
///////////////////////////////////////////////////////////////////////////////
// File:        math_operands.h
// Purpose:     Interface of the linear algebra classes:
//					vectors, matrices, ... 
//					supports int, double and Complex values
// Author:      A.Obraz
///////////////////////////////////////////////////////////////////////////////

#ifndef __SMALL_MAT_OPERANDS
#define __SMALL_MAT_OPERANDS

#include <complex>
#include <cmath>
#include <vector>

#include "math_impexp.h"

#include "gen_exception.h"
#include "io_helpers.h"


typedef std::complex<double> t_Complex;
typedef std::complex<double> t_CompVal;
// TODO: rename and finally replace!!!
typedef std::vector<double> t_DblVec;

template<typename T> class t_Matrix;
typedef t_Matrix<t_Complex> t_MatCmplx;
typedef t_Matrix<double> t_MatDbl;

template<typename T> class t_Vec;
typedef t_Vec<t_Complex> t_VecCmplx;
typedef t_Vec<double> t_VecDbl;

template<typename T> class t_Vec3;
typedef t_Vec3<t_Complex> t_Vec3Cmplx;
typedef t_Vec3<double> t_Vec3Dbl;

template<typename T> class t_SqMatrix;
typedef t_SqMatrix<t_Complex> t_SqMatCmplx;
typedef t_SqMatrix<double> t_SqMatDbl;

template<typename T> class t_SqMat3;
typedef t_SqMat3<t_Complex> t_SqMat3Cmplx;
typedef t_SqMat3<double> t_SqMat3Dbl;

/************************************************************************/
/* Basic linear algebra operators
** used to instantiate vectors of complex or double or smth constructable
** from 0.0
*/
/************************************************************************/

// sad but true, but in std std::norm(10)=100
// TODO: wrap everything with smat
namespace smat{
	IMPEXP_SMALLMAT double norm(t_Complex val);
}


namespace matrix{
	inline double herm_prod(double l, double r){ return l*r;};
	inline t_Complex herm_prod(t_Complex l, t_Complex r){return l*std::conj(r);};

	// type guess
	template<typename t1, typename t2>class TypeDeduce{};

	template<> class TypeDeduce<int, int>{public:typedef int type;};
	template<> class TypeDeduce<int, double>{public:typedef double type;};
	template<> class TypeDeduce<int, t_Complex>{public:typedef t_Complex type;};

	template<> class TypeDeduce<double, int>{public:typedef double type;};
	template<> class TypeDeduce<double, double>{public:typedef double type;};
	template<> class TypeDeduce<double, t_Complex>{public:typedef t_Complex type;};

	template<> class TypeDeduce<t_Complex, int>{public:typedef t_Complex type;};
	template<> class TypeDeduce<t_Complex, double>{public:typedef t_Complex type;};
	template<> class TypeDeduce<t_Complex, t_Complex>{public:typedef t_Complex type;};

};

/************************************************************************/
/* matrix class
** operations and storage are column-based
** _cont is a set of vectors===columns
** multiplication implies ordinary right mult: A.mul(B) <=> 
** C=A*B set again by columns: C[i][j] - column i, row j
** if matrix A is mul by set of vectors  h={h1, h2,}
** this means h*A ( like "basis" transition)
*/
/************************************************************************/

template<typename T> class  t_Matrix;

/*
template<typename T> t_Matrix<T> operator*
(const t_Matrix<T> & l, const t_Matrix<T>& r);
*/

template<typename T> class t_Vec;

template<typename T> class  t_Matrix{
public:
	struct t_Col{
		std::vector<T> _cont;

		t_Col(int dim, T val=0.0):_cont(dim, val){};
		std::vector<T>& std(){return _cont;};
		const std::vector<T>& std()const{return _cont;};
		//operators
		T& operator[](int n){ return _cont[n];};
		const T& operator[](int n) const{ return _cont[n];};
		t_Col& operator=(const t_Vec<T>& r){
			if (this->size()!=r.size())
				ssuTHROW(matrix::t_SizeMismatch, 
				_T("Col assign error: col vec size mismatch"));
			for (int i=0; i<this->size(); i++){
				_cont[i] = r[i];
			};
			return *this;	
		}; ;
		int size()const{return _cont.size();};
	};
protected:
	std::vector<t_Col> _cont;
	void _chk_col_ind(const int n) const;
	void _chk_row_ind(const int n) const; 
public:

	t_Matrix(const int a_nVecs=0, const int a_nElemInVec=0, T val = 0.0);
	t_Matrix(const int a_nVecs, const t_Col& col);

	int nCols() const{return _cont.size();};
	int nRows() const{if (nCols()>0) {return _cont[0].size();} else {return 0;}; };
	t_Col& operator[](const int n);
	const t_Col& operator [](const int n) const;

	friend std::ostream& operator<<(std::ostream& str, const t_Matrix& m){
		for (int i=0; i<m.nRows(); i++){
			str<<"\n[";
			for (int j=0; j<m.nCols(); j++){
				str<<std_manip::std_format_fixed<T>(m[j][i])<<"  ";
			}
			str<<"]";
		}
		return str;
	}
};

template<typename T> t_Matrix<T>::t_Matrix<T>
(const int a_nVecs, const int a_nElemInVec, T val /*=0.0*/)
:_cont(a_nVecs, t_Col(a_nElemInVec, val)){};

template<typename T>t_Matrix<T>::t_Matrix<T>
(const int a_nVecs, const t_Col& col)
:_cont(a_nVecs, col){};


template<typename T> void t_Matrix<T>::_chk_col_ind(const int n) const{
	if ((n<0)||(n>=nCols())) 
		ssuTHROW(matrix::t_ColOutOFRange, _T("Matrix Error: Column index out of range"));
};

template<typename T> void t_Matrix<T>::_chk_row_ind(const int m) const{
	if ((m<0)||(m>=nRows())) 
		ssuTHROW(matrix::t_RowOutOFRange, _T("Matrix Error: Row index out of range"));
};

template<typename T> 
const t_Matrix<T>::t_Col& t_Matrix<T>::operator [](const int n) const{
	_chk_col_ind(n);
	return _cont[n];
};

template<typename T>
t_Matrix<T>::t_Col& t_Matrix<T>::operator [](const int n){
	_chk_col_ind(n);
	return _cont[n];
};

/************************************************************************/
/* Matrix operators
** with vectors, square matrix and everything else
** which is derived from matrix
*/
/************************************************************************/

// this is in fact static t_Matrix methods
// I just want symmetric notation matrix::call_smth<t1,t2>(l,r);
// instead of t_Matrix<t1>::call_smth<t2>(l,r);
namespace matrix{
// exceptions
	class t_ColOutOFRange: public t_GenException{
	public:
		t_ColOutOFRange(const wxString& what, const wxChar* szFile,  const int line):
		  t_GenException(what, szFile, line){};
	};
	class t_RowOutOFRange: public t_GenException{
	public:
		t_RowOutOFRange(const wxString& what, const wxChar* szFile,  const int line):
		  t_GenException(what, szFile, line){};
	};
	class t_SizeMismatch: public t_GenException{
	public:
		t_SizeMismatch(const wxString& what, const wxChar* szFile,  const int line):
		  t_GenException(what, szFile, line){};
	};

	class t_BadIndex: public t_GenException{
	public:
		t_BadIndex(const wxString& what, const wxChar* szFile,  const int line):
		  t_GenException(what, szFile, line){};
	};

	class t_NotSquare: public t_GenException{
	public:
		t_NotSquare(const wxString& what, const wxChar* szFile,  const int line):
		  t_GenException(what, szFile, line){};
	};

	class t_Not3D: public t_GenException{
	public:
		t_Not3D(const wxString& what, const wxChar* szFile,  const int line):
		  t_GenException(what, szFile, line){};
	};

// base
	namespace base{
		template<typename t1, typename t2>
		t_Matrix<TypeDeduce<t1,t2>::type > mul
			(const t1 num, const t_Matrix<t2>& mat){
				t_Matrix<TypeDeduce<t1,t2>::type > ret(mat.nCols(), mat.nRows());
				for (int i=0; i<mat.nCols(); i++)
					for (int j=0; j<mat.nRows(); j++)
						ret[i][j]=mat[i][j]*num;
				return ret;
		};

		template<typename t1, typename t2>
		t_Matrix<TypeDeduce<t1,t2>::type > div
			(const t_Matrix<t1>& mat, const t2 num){
				return mul(1.0/num, mat);
		};

		template<typename t1, typename t2> 
		void chk_size_match_add
			(const t_Matrix<t1>& l, const t_Matrix<t2>& r){
				if ((l.nCols()!=r.nCols()) || (l.nRows()!=r.nRows()))
					ssuTHROW(matrix::t_SizeMismatch, 
						_T("Matrix Error: Add size mismatch"));
		};

		template<typename t1, typename t2> 
		void chk_size_match_mul
			(const t_Matrix<t1>& l, const t_Matrix<t2>& r){
				if (l.nCols()!=r.nRows())
					ssuTHROW(matrix::t_SizeMismatch, 
						_T("Matrix Error: Mul size mismatch"));
		};

		template<typename t1, typename t2>
		t_Matrix<TypeDeduce<t1,t2>::type > plus
			(const t_Matrix<t1>& l, const t_Matrix<t1>& r){

				chk_size_match_add<t1,t2>(l,r);
				t_Matrix<TypeDeduce<t1,t2>::type > ret(l.nCols(), l.nRows());

				for (int i=0; i<l.nCols(); i++)
					for (int j=0; j<l.nRows(); j++)
						ret[i][j]=l[i][j] + r[i][j];
				return ret;
		};

		template<typename t1, typename t2>
		t_Matrix<TypeDeduce<t1,t2>::type > minus
			(const t_Matrix<t1>& l, const t_Matrix<t1>& r){

				chk_size_match_add<t1,t2>(l,r);
				t_Matrix<TypeDeduce<t1,t2>::type > ret(l.nCols(), l.nRows());

				for (int i=0; i<l.nCols(); i++)
					for (int j=0; j<l.nRows(); j++)
						ret[i][j]=l[i][j] - r[i][j];
				return ret;
		};


		template<typename t1, typename t2>
		t_Matrix<TypeDeduce<t1,t2>::type > mat_mul
			(const t_Matrix<t1>& l, const t_Matrix<t2>& r){

				typedef TypeDeduce<t1,t2>::type type;

				chk_size_match_mul<t1,t2>(l,r);
				t_Matrix<type> ret(r.nCols(), l.nRows());

				for (int i=0; i<r.nCols(); i++){
					for (int j=0; j<l.nCols(); j++){
						t_Matrix<type>::t_Col tmp(l[j].size());
						type val = r[i][j];

						for (int k=0; k<tmp.size(); k++){
							tmp[k]=l[j][k];
							tmp[k]*=val;
						};
						for (int k=0; k<ret[i].size(); k++){
							ret[i][k]+=tmp[k]; 
						};
					};
				};
				return ret;
		};
	};
};

//matrix-specific operators

template<typename t1, typename t2>
t_Matrix<matrix::TypeDeduce<t1,t2>::type > operator*
	(const t_Matrix<t1>& mat, const t2 num){
		return matrix::base::mul<t1,t2>(num, mat);
}

template<typename t1, typename t2>
t_Matrix<matrix::TypeDeduce<t1,t2>::type > operator*
(const t2 num, const t_Matrix<t1>& mat){
	return matrix::base::mul<t1,t2>(num, mat);
}

template<typename t1, typename t2>
t_Matrix<matrix::TypeDeduce<t1,t2>::type > operator*
(const t_Matrix<t1>& l, const t_Matrix<t2>& r){
	return matrix::base::mat_mul<t1,t2>(l,r);
}

template<typename t1, typename t2>
t_Matrix<matrix::TypeDeduce<t1,t2>::type > operator/
(const t_Matrix<t1>& l, const t2 r){
	return matrix::base::div<t1,t2>(l,r);
}

template<typename t1, typename t2>
t_Matrix<matrix::TypeDeduce<t1,t2>::type > operator+
(const t_Matrix<t1>& l, const t_Matrix<t2>& r){
	return matrix::base::plus<t1,t2>(l,r);
}

template<typename t1, typename t2>
t_Matrix<matrix::TypeDeduce<t1,t2>::type > operator-
(const t_Matrix<t1>& l, const t_Matrix<t2>& r){
	return matrix::base::minus<t1,t2>(l,r);
}

/************************************************************************/
/* Linear Algebra Vector class
** with a second norm as norm
*/
/************************************************************************/

namespace vector{
	template<typename t1, typename t2>
	matrix::TypeDeduce<t1,t2>::type dot
		(const t_Vec<t1>&, const t_Vec<t2>&);
}

template<typename T> class  t_Vec: public t_Matrix<T>{
	t_Matrix<T>::t_Col& _get_col(){return _cont[0];};
	const t_Matrix<T>::t_Col& _get_col() const{return _cont[0];};
public:
	t_Vec(const int dim=0, const T val=T(0.0)):t_Matrix(1, dim, val){};
	t_Vec(const t_Matrix<T>::t_Col& col):t_Matrix(1, col){};

	int size() const{
		const t_Matrix<T>::t_Col& col = _cont[0];
		return col.size();
	};
	// TODO: maybe some check that imag is small ?
	T norm() const{return sqrt(vector::dot(*this, *this));};
	void normalize();
	// TODO: check index? hmm... if something is wrong fail anyway)
	inline T& operator[](const int ind);
	inline const T& operator[](const int ind) const;

};

template<typename T>T& t_Vec<T>::operator[](const int ind){
	t_Matrix<T>::t_Col& col = _cont[0];
	if ((ind>=0)&&(ind<size())) return col[ind];
	ssuTHROW(matrix::t_BadIndex, _T("t_Vec error: Index out of range"));
};
template<typename T>const T& t_Vec<T>::operator[](const int ind) const{
	const t_Matrix<T>::t_Col& col = _cont[0];
	if ((ind>=0)&&(ind<size())) return col[ind];
	ssuTHROW(matrix::t_BadIndex, _T("t_Vec error: Index out of range"));
};

template<typename T>void t_Vec<T>::normalize(){
	double norm = this->norm();
	std::vector<T>& col = _cont[0].std();
	for (std::vector<T>::iterator it = col.begin(); it<col.end(); it++){
		*it/=norm;
	};
};

namespace vector{
	template<typename t1, typename t2> 
		matrix::TypeDeduce<t1,t2>::type dot
		(const t_Vec<t1>& l, const t_Vec<t2>& r) {
			typedef matrix::TypeDeduce<t1,t2>::type type;
			type ret = 0;
			matrix::base::chk_size_match_add(l, r);
			for (int i=0; i<l.size(); i++){
				ret+=matrix::herm_prod(l[i], r[i]); 
			}
			return ret;
	}; 

}

//Vec-specific operators

template<typename t1, typename t2>
t_Vec<matrix::TypeDeduce<t1,t2>::type > operator*
(const t_Vec<t1>& vec, const t2 num){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_Vec<type> ret(matrix::base::mul<t1,t2>(num, vec)[0]);
	return ret;
}

template<typename t1, typename t2>
t_Vec<matrix::TypeDeduce<t1,t2>::type > operator*
(const t1 num, const t_Vec<t2>& vec){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_Vec<type> ret(matrix::base::mul<t1,t2>(num, vec)[0]);
	return ret;
}

template<typename t1, typename t2>
t_Vec<matrix::TypeDeduce<t1,t2>::type > operator/
(const t_Vec<t1>& vec, const t2 num){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_Vec<type> ret(matrix::base::div<t1,t2>(vec, num)[0]);
	return ret;
}

template<typename t1, typename t2>
t_Vec<matrix::TypeDeduce<t1,t2>::type > operator+
(const t_Vec<t1>& l, const t_Vec<t2>& r){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_Vec<type> ret(matrix::base::plus<t1,t2>(l,r)[0]);
	return ret;
}

template<typename t1, typename t2>
t_Vec<matrix::TypeDeduce<t1,t2>::type > operator-
(const t_Vec<t1>& l, const t_Vec<t2>& r){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_Vec<type> ret(matrix::base::minus<t1,t2>(l,r)[0]);
	return ret;
}


// A*h=b
template<typename t1, typename t2>
t_Vec<matrix::TypeDeduce<t1,t2>::type > operator*
(const t_Matrix<t1>& l, const t_Vec<t2>& r){
	typedef matrix::TypeDeduce<t1,t2>::type type;
 	t_Vec<type> ret(matrix::base::mat_mul<t1,t2>(l,r)[0]);
	return ret;
}

/************************************************************************/
/* 3D vector                                                            */
/************************************************************************/


template<typename T> class t_Vec3 : public t_Vec<T>{
public:
	t_Vec3(T x=0.0, T y=0.0, T z=0.0):t_Vec(3){
		this->operator[](0) = x;
		this->operator[](1) = y;
		this->operator[](2) = z;
	};
	t_Vec3(const t_Vec<T>& r): t_Vec(r){
		if (r.size()!=3) 
			ssuTHROW(matrix::t_Not3D, 
				_T("3D vector error: init by non-3D general vector"));
	};
	t_Vec3& operator=(const t_Vec<T>& r){
		if (r.size()!=3) 
			ssuTHROW(matrix::t_Not3D, 
				_T("3D vector error: init by non-3D general vector"));
		(t_Vec&)(*this) = r;
		return *this;
	};
};

namespace vector{
	template<typename t1, typename t2> 
	t_Vec3<matrix::TypeDeduce<t1,t2>::type> cross
		(const t_Vec3<t1>& l, const t_Vec3<t2>& r){
			typedef matrix::TypeDeduce<t1,t2>::type type;
			type v0 = l[1]*r[2] - l[2]*r[1];
			type v1 = l[2]*r[0] - l[0]*r[2];
			type v2 = l[0]*r[1] - l[1]*r[2];
			return t_Vec3<type>(v0,v1,v2);
	};
}

//3d-vector-specific operators

template<typename t1, typename t2>
t_Vec3<matrix::TypeDeduce<t1,t2>::type > operator*
(const t_Vec3<t1>& vec, const t2 num){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_Vec3<type> ret(matrix::base::mul<t1,t2>(num, vec)[0]);
	return ret;
}

template<typename t1, typename t2>
t_Vec3<matrix::TypeDeduce<t1,t2>::type > operator*
(const t2 num, const t_Vec3<t1>& vec){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_Vec3<type> ret(matrix::base::mul<t1,t2>(num, vec)[0]);
	return ret;
}

template<typename t1, typename t2>
t_Vec3<matrix::TypeDeduce<t1,t2>::type > operator/
(const t_Vec3<t1>& vec, const t2 num){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_Vec3<type> ret(matrix::base::div<t1,t2>(vec, num)[0]);
	return ret;
}

template<typename t1, typename t2>
t_Vec3<matrix::TypeDeduce<t1,t2>::type > operator+
(const t_Vec3<t1>& l, const t_Vec3<t2>& r){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_Vec3<type> ret(matrix::base::plus<t1,t2>(l,r)[0]);
	return ret;
}

template<typename t1, typename t2>
t_Vec3<matrix::TypeDeduce<t1,t2>::type > operator-
(const t_Vec3<t1>& l, const t_Vec3<t2>& r){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_Vec3<type> ret(matrix::base::minus<t1,t2>(l,r)[0]);
	return ret;
}



/************************************************************************/
/* Square Matrix                                                        */
/************************************************************************/
template<typename T>class  t_SqMatrix: public t_Matrix<T>{
public:
	t_SqMatrix():t_Matrix(0,0){};
	t_SqMatrix(const int dim):t_Matrix(dim, dim){};;
	t_SqMatrix(const t_Matrix<T>&);
	int size() const{return nCols();};
	void resize(const int n){_cont.resize(n, t_Col(n, 0.0));};
	// construct a minor from matrix
	//t_SqMatrix(const t_SqMatrix& mat, const int raw, const int col);
	//instead use this
	// TODO: better col first
	t_SqMatrix<T> get_minor(const int raw, const int col) const;
	void setToUnity();
	//void resize(int a_size);
	T det() const;
	t_SqMatrix<T> inverse() const;
};

template<typename T> t_SqMatrix<T>::t_SqMatrix<T>(const t_Matrix<T>& m){
	if (m.nCols()==m.nRows()){
		(t_Matrix&)*this = m;
	}else{
		ssuTHROW(matrix::t_NotSquare, 
			_T("Square matrix error: init by non-square matrix"));
	}
}

template<typename T> t_SqMatrix<T> t_SqMatrix<T>::get_minor
(const int row, const int col) const{
	_chk_col_ind(col);
	_chk_row_ind(row);
	t_SqMatrix<T> ret(size()-1);
	for (int i=0; i<col; i++){
		for (int j=0; j<row; j++){
			ret[i][j] = _cont[i][j];
		};
		for (int j=row+1; j<nRows(); j++){
			ret[i][j-1] = _cont[i][j];
		};
	};
	for (int i=col+1; i<nCols(); i++){
		for (int j=0; j<row; j++){
			ret[i-1][j] = _cont[i][j];
		};
		for (int j=row+1; j<nRows(); j++){
			ret[i-1][j-1] = _cont[i][j];
		};
	};
	return ret;
};

// calculate determinant by first row (0):
// det = SUM(j)((-1)^(0+j)*det(d(0,j)))
// d(0,j) - algebraic minor
// TODO: rewrite with gauss (QR)!!!
template<typename T> T t_SqMatrix<T>::det() const{
	if (this->nCols()==1){
		return _cont[0][0];
	};
	// optimize for small matrices 2x2, 3x3
	if (this->nCols()==2){
		return _cont[0][0]*_cont[1][1] - 
			_cont[0][1]*_cont[1][0];
	};
	if (this->nCols()==3){
		return 
			_cont[0][0]*
			(
			_cont[1][1]*_cont[2][2] - 
			_cont[2][1]*_cont[1][2]
		) - 
			_cont[1][0]*
			(
			_cont[0][1]*_cont[2][2] - 
			_cont[2][1]*_cont[0][2]
		) +
			_cont[2][0]*
			(
			_cont[0][1]*_cont[1][2]-
			_cont[1][1]*_cont[0][2]
		);
	};
	T res(0.0);
	for (int i=0; i<nCols(); i++){
		t_SqMatrix minor = get_minor(0, i);
		T val = _cont[i][0];
		res+= pow(-1.0, i)*val*minor.det();
	}
	return res;
}

template<typename T> void t_SqMatrix<T>::setToUnity(){
	for (int i=0; i<this->size(); i++){
		for (int j=0; j<this->size(); j++){
			if (i==j) 
				_cont[i][j]=1.0;
			else 
				_cont[i][j]=0.0;
		};
	};
};
// TODO: rewrite with robust procedures!!!
template<typename T> t_SqMatrix<T> t_SqMatrix<T>::inverse() const{
	int dim = this->nCols();
	t_SqMatrix<T> ret(dim);
	t_SqMatrix<T> unity(dim);
	unity.setToUnity();
	T coef = 1.0/this->det();
	for (int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			t_SqMatrix temp = *this;
			temp[j] = unity[i];
			ret[i][j] = coef*temp.det();
		}
	}
	return ret;
};

//sq-matrix-specific operators

template<typename t1, typename t2>
t_SqMatrix<matrix::TypeDeduce<t1,t2>::type > operator*
(const t_SqMatrix<t1>& mat, const t2 num){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMatrix<type> ret(matrix::base::mul<t1,t2>(num, mat));
	return ret;
}

template<typename t1, typename t2>
t_SqMatrix<matrix::TypeDeduce<t1,t2>::type > operator*
(const t2 num, const t_SqMatrix<t1>& mat){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMatrix<type> ret(matrix::base::mul<t1,t2>(num, mat));
	return ret;
}

template<typename t1, typename t2>
t_SqMatrix<matrix::TypeDeduce<t1,t2>::type > operator*
(const t_SqMatrix<t1>& l, const t_SqMatrix<t2>& r){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMatrix<type> ret(matrix::base::mat_mul<t1,t2>(l,r));
	return ret;
}

template<typename t1, typename t2>
t_SqMatrix<matrix::TypeDeduce<t1,t2>::type > operator/
(const t_SqMatrix<t1>& l, const t2 r){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMatrix<type> ret(matrix::base::div<t1,t2>(l,r));
	return ret;
}

template<typename t1, typename t2>
t_SqMatrix<matrix::TypeDeduce<t1,t2>::type > operator+
(const t_SqMatrix<t1>& l, const t_SqMatrix<t2>& r){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMatrix<type> ret(matrix::base::plus<t1,t2>(l,r));
	return ret;
}

template<typename t1, typename t2>
t_SqMatrix<matrix::TypeDeduce<t1,t2>::type > operator-
(const t_SqMatrix<t1>& l, const t_SqMatrix<t2>& r){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMatrix<type> ret(matrix::base::minus<t1,t2>(l,r));
	return ret;
}

/************************************************************************/
/* 3D matrix (orthogonal transformation matrix)                         */
/************************************************************************/

template<typename T> class t_SqMat3 : public t_SqMatrix<T>{
public:
	enum t_Init{Zero=0, Unity};
	t_SqMat3(t_Init init_type = Zero):t_SqMatrix(3){
		if (init_type==Unity){
			_cont[0][0]=1.0;
			_cont[1][1]=1.0;
			_cont[2][2]=1.0;
		};
	};
	
	t_SqMat3 inverse() const{
		t_SqMat3 ret;
		const t_SqMatrix& rM = *this;
		(t_SqMatrix&)ret = rM.inverse();
		return ret;
	};

};

// 3d-matrix-specific operators

template<typename t1, typename t2>
t_SqMat3<matrix::TypeDeduce<t1,t2>::type > operator*
(const t_SqMat3<t1>& mat, const t2 num){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMat3<type> ret(matrix::base::mul<t1,t2>(num, mat));
	return ret;
}

template<typename t1, typename t2>
t_SqMat3<matrix::TypeDeduce<t1,t2>::type > operator*
(const t1 l, const t_SqMat3<t2>& r){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMat3<type> ret(matrix::base::mul<t1,t2>(l,r));
	return ret;
}

template<typename t1, typename t2>
t_Vec3<matrix::TypeDeduce<t1,t2>::type > operator*
(const t_SqMat3<t1>& l, const t_Vec3<t2>& r){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_Vec3<type> ret(matrix::base::mat_mul<t1,t2>(l,r)[0]);
	return ret;
}

template<typename t1, typename t2>
t_SqMat3<matrix::TypeDeduce<t1,t2>::type > operator*
(const t_SqMat3<t1>& l, const t_SqMat3<t1>& r){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMat3<type> ret(matrix::base::mat_mul<t1,t2>(l,r));
	return ret;
}

template<typename t1, typename t2>
t_SqMat3<matrix::TypeDeduce<t1,t2>::type > operator/
(const t_SqMat3<t1>& l, const t2 r){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMat3<type> ret(matrix::base::div<t1,t2>(l,r));
	return ret;
}

template<typename t1, typename t2>
t_SqMat3<matrix::TypeDeduce<t1,t2>::type > operator+
(const t_SqMat3<t1>& l, const t_SqMat3<t2>& r){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMat3<type> ret(matrix::base::plus<t1,t2>(l,r));
	return ret;
}

template<typename t1, typename t2>
t_SqMat3<matrix::TypeDeduce<t1,t2>::type > operator-
(const t_SqMat3<t1>& l, const t_SqMat3<t2>& r){
	typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMat3<type> ret(matrix::base::minus<t1,t2>(l,r));
	return ret;
}

#endif // __SMALL_MAT_OPERANDS
