///////////////////////////////////////////////////////////////////////////////
// Project:	SmallMat
// Purpose:	various math routines
///////////////////////////////////////////////////////////////////////////////
// File:        math_operands.h
// Purpose:     Interface of the linear algebra classes:
//					vectors, matrices, ... 
//					supports int, double and Complex values
// Author:      A.Obraz
///////////////////////////////////////////////////////////////////////////////

#ifndef __MATH_OPERANDS
#define __MATH_OPERANDS

#include <complex>
#include <cmath>
#include <vector>

#include "dll_impexp_smat.h"

#include "gen_exception.h"
#include "io_helpers.h"

#include "type_deduce.h"

typedef std::complex<double>  t_Complex;
typedef std::complex<double> t_CompVal;

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

template<typename T>  class t_Vec;

template<typename T> class t_Matrix{
protected:
	T** _cont;
	const int _ncols, _nrows;
	void _chk_col_ind(const int n) const;
	void _chk_row_ind(const int n) const; 
	t_Matrix& _set(const t_Matrix& a_mat);
public:

	t_Matrix(const int a_nVecs=0, const int a_nElemInVec=0, T val = 0.0);
	t_Matrix(const t_Matrix& a_mat);

	~t_Matrix();

	inline int nCols() const{return _ncols;};
	inline int nRows() const{return _nrows;};

	inline T* operator[](int n);
	void col_to_vec(int col, t_Vec<T>& vec ) const;
	void set_col(int col, const t_Vec<T>& vec);
	void setToZero();
	t_Matrix& mul_by_factor(T val);
	// set is alias for operator= 
	// use when operator= can be ambiguous
	t_Matrix& set(const t_Matrix& a_mat);
	t_Matrix& operator=(const t_Matrix& a_mat);
	const T* operator[](int n) const;

	friend std::wostream& operator<<(std::wostream& str, const t_Matrix& m){
		for (int i=0; i<m.nRows(); i++){
			str<<_T("\n[");
			for (int j=0; j<m.nCols(); j++){
				str<<std_manip::std_format_fixed<T>(m[j][i])<<_T("  ");
			}
			str<<_T("]");
		}
		return str;
	}

	wxString to_wxstr() const {
		std::wostringstream ostr;
		ostr << *this;
		return wxString(&ostr.str()[0]);
	}
};

//template<typename T> t_Matrix<T>::t_Matrix<T>
//(const int a_nVecs, const int a_nElemInVec, T val /*=0.0*/)
//:_cont(a_nVecs, t_Col(a_nElemInVec, val)){};

template<typename T> t_Matrix<T>::t_Matrix(const int a_nVecs, const int a_nElemInVec, T val /*=0.0*/)
:_ncols(a_nVecs), _nrows(a_nElemInVec){
	_cont = new T*[_ncols];
	for (int i=0; i<_ncols; i++) _cont[i] = new T[_nrows];
}

template<typename T> t_Matrix<T>::t_Matrix(const t_Matrix& a_mat)
:_ncols(a_mat.nCols()), _nrows(a_mat.nRows()){
	_cont = new T*[_ncols];
	for (int i=0; i<_ncols; i++) _cont[i] = new T[_nrows];
	for (int i=0; i<_ncols; i++)
		for (int j=0; j<_nrows; j++)
			_cont[i][j] = a_mat[i][j];
}

template<typename T>t_Matrix<T>::~t_Matrix(){
	for (int i=0; i<_ncols; i++) delete[] _cont[i];
	delete[] _cont;
}

// assign column to a destination vector
// 
template<typename T> void t_Matrix<T>::col_to_vec(int col, t_Vec<T>& vec ) const{
#ifdef _DEBUG
	if (vec.size()!=_nrows) ssuTHROW(matrix::t_SizeMismatch, 
		_T("Matrix col to Vec error: wrong size of dest vec"));
#endif
	for (int j=0; j<_nrows; j++) vec[j] = _cont[col][j];
}

// set column values
template<typename T> void t_Matrix<T>::set_col(int col, const t_Vec<T>& vec){
#ifdef _DEBUG
	if (vec.size()!=_nrows) ssuTHROW(matrix::t_SizeMismatch, 
		_T("Matrix col to Vec error: wrong size of dest vec"));
#endif
	for (int j=0; j<_nrows; j++) _cont[col][j] = vec[j];
}

template<typename T> void t_Matrix<T>::setToZero(){
	for (int i=0; i<_ncols; i++) 
		for (int j=0; j<_nrows; j++) 
			_cont[i][j]=0;
}

template<typename T> t_Matrix<T>& t_Matrix<T>::mul_by_factor(T val){
	for (int i=0; i<_ncols; i++) 
		for (int j=0; j<_nrows; j++) 
			_cont[i][j]*=val;

	return *this;
}

template<typename T> inline t_Matrix<T>& t_Matrix<T>::_set(const t_Matrix<T>& a_mat){
#ifdef _DEBUG
	matrix::base::chk_size_match_add(*this, a_mat);
#endif
	for (int i=0; i<_ncols; i++)
		for (int j=0; j<_nrows; j++)
			_cont[i][j] = a_mat[i][j];
	return *this;
}

template<typename T> t_Matrix<T>& t_Matrix<T>::set(const t_Matrix<T>& a_mat){
	return this->_set(a_mat);
};
template<typename T> t_Matrix<T>& t_Matrix<T>::operator=(const t_Matrix<T>& a_mat){
	return this->_set(a_mat);
};
/************************************************************************/
/* Exceptions
*/
/************************************************************************/

namespace matrix{

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


}

template<typename T>inline void t_Matrix<T>::_chk_col_ind(const int n) const{
	if ((n<0)||(n>=nCols())) 
		ssuTHROW(matrix::t_ColOutOFRange, _T("Matrix Error: Column index out of range"));
};

template<typename T>inline void t_Matrix<T>::_chk_row_ind(const int m) const{
	if ((m<0)||(m>=nRows())) 
		ssuTHROW(matrix::t_RowOutOFRange, _T("Matrix Error: Row index out of range"));
};

template<typename T> inline const T* t_Matrix<T>::operator [](int n) const{
#ifdef _DEBUG
	_chk_col_ind(n);
#endif
	return _cont[n];
};

template<typename T> inline T* t_Matrix<T>::operator [](int n){
#ifdef _DEBUG
	_chk_col_ind(n);
#endif
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

	// base
	namespace base{
		template<typename t1, typename t2> void mul
			(const t1 num, const t_Matrix<t2>& mat, 
			t_Matrix<typename TypeDeduce<t1,t2>::type >& ret){

				int ncols = mat.nCols();
				int nrows = mat.nRows();
#ifdef _DEBUG
				if(ret.nCols()!=ncols||ret.nRows()!=nrows)
					ssuTHROW(matrix::t_SizeMismatch, _T("Mat mul scalar error : sizes do not match"));
#endif
				for (int i=0; i<ncols; i++)
					for (int j=0; j<nrows; j++)
						ret[i][j]=mat[i][j]*num;
		};

		template<typename t1, typename t2> void div
			(const t_Matrix<t1>& mat, const t2 num,
			t_Matrix<typename TypeDeduce<t1,t2>::type >& ret){
				mul(1.0/num, mat, ret);
		};

		template<typename t1, typename t2> inline
			void chk_size_match_add
			(const t_Matrix<t1>& l, const t_Matrix<t2>& r){
				if ((l.nCols()!=r.nCols()) || (l.nRows()!=r.nRows()))
					ssuTHROW(matrix::t_SizeMismatch, 
					_T("Matrix Error: Add size mismatch"));
		};

		template<typename t1, typename t2> inline
			void chk_size_match_mul
			(const t_Matrix<t1>& l, const t_Matrix<t2>& r){
				if (l.nCols()!=r.nRows())
					ssuTHROW(matrix::t_SizeMismatch, 
					_T("Matrix Error: Mul size mismatch"));
		};

		template<typename t1, typename t2> void plus
			(const t_Matrix<t1>& l, const t_Matrix<t1>& r,
			t_Matrix<typename TypeDeduce<t1,t2>::type>& ret){
#ifdef _DEBUG
				typedef TypeDeduce<t1,t2>::type type; 
				chk_size_match_add<t1,t2>(l,r);
				chk_size_match_add<t1,type>(l, ret);
#endif
				int ncols, nrows;
				ncols = l.nCols();
				nrows = l.nRows();

				for (int i=0; i<ncols; i++)
					for (int j=0; j<nrows; j++)
						ret[i][j]=l[i][j] + r[i][j];
				return;
		};

		template<typename t1, typename t2> void minus
			(const t_Matrix<t1>& l, const t_Matrix<t1>& r,
			t_Matrix<typename TypeDeduce<t1,t2>::type>& ret){
#ifdef _DEBUG
				typedef TypeDeduce<t1,t2>::type type; 
				chk_size_match_add<t1,t2>(l,r);
				chk_size_match_add<t1,type>(l, ret);
#endif
				int ncols, nrows;
				ncols = l.nCols();
				nrows = l.nRows();

				for (int i=0; i<ncols; i++)
					for (int j=0; j<nrows; j++)
						ret[i][j]=l[i][j] - r[i][j];
				return;
		};
		
		

		template<typename t1, typename t2, typename tt> void __mat_mul
			(const t_Matrix<t1>& l, const t_Matrix<t2>& r, 
			t_Matrix<tt >& ret){
			
				typedef tt type;
				
				int ncols, nrows;
				ncols = r.nCols();
				nrows = l.nRows();
#ifdef _DEBUG
				chk_size_match_mul<t1,t2>(l,r);
				if (ret.nCols()!=ncols||ret.nRows()!=nrows){
				
					wxLogMessage(_T("l:c=%d;r=%d"), l.nRows(), l.nCols());
					wxLogMessage(_T("r:c=%d;r=%d"), r.nRows(), r.nCols());
					wxLogMessage(_T("ret:c=%d;r=%d"), ret.nRows(), ret.nCols());
					ssuTHROW(matrix::t_SizeMismatch, _T("Mat mul error: destination matrix has wrong size"));
				}
#endif
				//type* tmp = new type[nrows];
				//type val; 
				ret.setToZero();
				/*
				for (int i=0; i<r.nCols(); i++){
					for (int j=0; j<l.nCols(); j++){
						val = r[i][j];

						for (int k=0; k<nrows; k++){
							tmp[k]=l[j][k];
							tmp[k]*=val;
						};
						for (int k=0; k<nrows; k++){
							ret[i][k]+=tmp[k]; 
						};
					};
				};
				delete[] tmp;
				return;
				*/
				for (int i = 0; i < nrows; i++)
					for (int j = 0; j < ncols; j++)
						for (int k = 0; k < l.nCols(); k++)
							ret[j][i] += l[k][i] * r[j][k];
		};
		
		template<typename t1, typename t2> void mat_mul
		    (const t_Matrix<t1>& l, const t_Matrix<t2>& r, 
		    t_Matrix<typename matrix::TypeDeduce<t1,t2>::type>& ret);
	};
};

//matrix-specific operators

template<typename t1, typename t2> 
t_Matrix<typename matrix::TypeDeduce<t1, t2>::type > operator*
(const t_Matrix<t1>& mat, const t2 num){
	t_Matrix<typename matrix::TypeDeduce<t1,t2>::type > ret(mat.nCols(), mat.nRows());
	return matrix::base::mul<t1,t2>(num, mat, ret);
	return ret;
}

template<typename t1, typename t2> 
t_Matrix<typename matrix::TypeDeduce<t1,t2>::type > operator*
(const t2 num, const t_Matrix<t1>& mat){
	t_Matrix<typename matrix::TypeDeduce<t1,t2>::type > ret(mat.nCols(), mat.nRows());
	matrix::base::mul<t1,t2>(num, mat, ret);
	return ret;
}

template<typename t1, typename t2> 
t_Matrix<typename matrix::TypeDeduce<t1,t2>::type > operator*
(const t_Matrix<t1>& l, const t_Matrix<t2>& r){
	t_Matrix<typename matrix::TypeDeduce<t1,t2>::type > ret(r.nCols(), l.nRows());
	matrix::base::mat_mul<t1,t2>(l,r, ret);
	return ret;
}

template<typename t1, typename t2> 
t_Matrix<typename matrix::TypeDeduce<t1,t2>::type > operator/
(const t_Matrix<t1>& l, const t2 r){
	t_Matrix<typename matrix::TypeDeduce<t1,t2>::type > ret(l.nCols(), l.nRows());
	matrix::base::div<t1,t2>(l, r, ret);
	return ret;
}

template<typename t1, typename t2>  
t_Matrix<typename matrix::TypeDeduce<t1,t2>::type > operator+
(const t_Matrix<t1>& l, const t_Matrix<t2>& r){
	t_Matrix<typename matrix::TypeDeduce<t1,t2>::type > ret(l.nCols(), r.nRows());
	matrix::base::plus<t1,t2>(l,r,ret);
	return ret;
}

template<typename t1, typename t2>  
t_Matrix<typename matrix::TypeDeduce<t1,t2>::type > operator-
(const t_Matrix<t1>& l, const t_Matrix<t2>& r){
	t_Matrix<typename matrix::TypeDeduce<t1,t2>::type > ret(l.nCols(), l.nRows());
	matrix::base::minus<t1,t2>(l,r,ret);
	return ret;
}

/************************************************************************/
/* Linear Algebra Vector class
** with a second norm as norm
*/
/************************************************************************/

namespace vector{
	template<typename t1, typename t2, typename TT> inline
		TT __dot(const t_Vec<t1>& l, const t_Vec<t2>& r) {
			//typedef matrix::TypeDeduce<t1,t2>::type type;
			TT ret = 0;
#ifdef _DEBUG
			matrix::base::chk_size_match_add(l, r);
#endif
			for (int i=0; i<l.size(); i++){
				ret+=matrix::herm_prod(l[i], r[i]); 
			}
			return ret;
	}; 

	template<typename t1, typename t2>
		typename matrix::TypeDeduce<t1,t2>::type 
			dot(const t_Vec<t1>& l, const t_Vec<t2>& r){
			    return __dot<t1, t2, matrix::TypeDeduce<t1,t2>::type>(l,r);
			};
	template<typename t1, typename t2, typename TT> inline
			TT __plain_prod(const t_Vec<t1>& l, const t_Vec<t2>& r) {
				//typedef matrix::TypeDeduce<t1,t2>::type type;
				TT ret = 0;
#ifdef _DEBUG
				matrix::base::chk_size_match_add(l, r);
#endif
				for (int i=0; i<l.size(); i++){
					ret+=l[i]*r[i]; 
				}
				return ret;
		}; 

	template<typename t1, typename t2>
	typename matrix::TypeDeduce<t1,t2>::type 
		plain_prod(const t_Vec<t1>& l, const t_Vec<t2>& r){
			return __plain_prod<t1, t2, matrix::TypeDeduce<t1,t2>::type>(l,r);
	};
}

template<typename T> class  t_Vec: public t_Matrix<T>{
	T* _get_col(){return this->_cont[0];};
	const T* _get_col() const{return this->_cont[0];};
public:
	t_Vec<T>(const int dim=0, const T val=0.0):t_Matrix<T>(1, dim, val){};

	int size() const{return t_Matrix<T>::_nrows;};
	void set_vals(T val=0.0);
	// TODO: maybe some check that imag is small ?
	T norm() const{return sqrt(vector::dot(*this, *this));};
	t_Vec& normalize();
	// TODO: check index? hmm... if something is wrong fail anyway)
	inline T& operator[](int ind);
	inline T operator[](int ind) const;

};

template<typename T>inline T& t_Vec<T>::operator[](int ind){
	//if ((ind>=0)&&(ind<size())) return _cont[0][ind];
	return t_Matrix<T>::_cont[0][ind];
	//ssuTHROW(matrix::t_BadIndex, _T("t_Vec error: Index out of range"));
};
template<typename T>inline T t_Vec<T>::operator[](int ind) const{
	//if ((ind>=0)&&(ind<size())) return _cont[0][ind];
	return this->_cont[0][ind];
	//ssuTHROW(matrix::t_BadIndex, _T("t_Vec error: Index out of range"));
};

template<typename T>void t_Vec<T>::set_vals(T val/*=0.0*/){
	for(int i=0; i<size(); i++) this->_cont[0][i] = val;
}

template<typename T>t_Vec<T>& t_Vec<T>::normalize(){
	T norm = this->norm();
	for (int i=0; i<this->_nrows; i++) this->_cont[0][i]/=norm;
	return *this;
};


//Vec-specific operators
// inline funcs are to avoid typedef type guesses by compiler 
// (not all compilers can easily do that)

template<typename t1, typename t2, typename TT> inline
t_Vec<TT> __vecbynum_mul(const t_Vec<t1>& vec, const t2 num){
//	typedef matrix::TypeDeduce<t1,t2>::type type;
	//t_Vec<type> ret(matrix::base::mul<t1,t2>(num, vec)[0]);
	t_Vec<TT> ret(vec);
	for (int i=0; i<ret.size(); i++) ret[i]*=num;
	return ret;
};

template<typename t1, typename t2> 
    t_Vec<typename matrix::TypeDeduce<t1,t2>::type> 
	operator*(const t_Vec<t1>& vec, const t2 num); 

template<typename t1, typename t2>
t_Vec<typename matrix::TypeDeduce<t1, t2>::type>
	operator*(const t1 num, const t_Vec<t2>& vec);


template<typename t1, typename t2, typename TT> inline 
t_Vec<TT> __numbyvec_mul(const t1 num, const t_Vec<t2>& vec){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	//t_Vec<type> ret(matrix::base::mul<t1,t2>(num, vec)[0]);
	//return ret;
	return __vecbynum_mul<t2,t1,TT>(vec, num);
};

template<typename t1, typename t2, typename TT> inline 
t_Vec<TT> __vecbynum_div(const t_Vec<t1>& vec, const t2 num){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	//t_Vec<type> ret(matrix::base::div<t1,t2>(vec, num)[0]);
	t_Vec<TT> ret(vec);
	for (int i=0; i<ret.size(); i++) ret[i]/=num;
	return ret;
};

template<typename t1, typename t2> 
    t_Vec<typename matrix::TypeDeduce<t1,t2>::type> 
	operator/(const t_Vec<t1>& vec, const t2 num); 


template<typename t1, typename t2, typename TT> inline
t_Vec<TT > __vec_plus(const t_Vec<t1>& l, const t_Vec<t2>& r){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	//t_Vec<type> ret(matrix::base::plus<t1,t2>(l,r)[0]);
#ifdef _DEBUG
	matrix::base::chk_size_match_add(l, r);
#endif
	t_Vec<TT> ret(l);
	for (int i=0; i<ret.size(); i++) ret[i]= l[i] + r[i];
	return ret;
}

template<typename t1, typename t2> 
    t_Vec<typename matrix::TypeDeduce<t1,t2>::type> 
	operator+(const t_Vec<t1>& l, const t_Vec<t2>& r); 

template<typename t1, typename t2, typename TT> inline 
t_Vec<TT> __vec_minus
(const t_Vec<t1>& l, const t_Vec<t2>& r){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	//t_Vec<type> ret(matrix::base::minus<t1,t2>(l,r)[0]);
#ifdef _DEBUG
	matrix::base::chk_size_match_add(l, r);
#endif
	t_Vec<TT> ret(l);
	for (int i=0; i<ret.size(); i++) ret[i]= l[i] - r[i];
	return ret;
};

template<typename t1, typename t2> 
    t_Vec<typename matrix::TypeDeduce<t1,t2>::type> 
	operator-(const t_Vec<t1>& l,const t_Vec<t2>& r); 


// A*h=b
template<typename t1, typename t2, typename TT> inline
t_Vec<TT> __matbyvec_mul(const t_Matrix<t1>& l, const t_Vec<t2>& r){
	// TODO: check sizes
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	t_Vec<TT> ret(r.size());
	//matrix::base::mat_mul<t1,t2>(l,r, ret);
	matrix::base::__mat_mul<t1,t2,TT>(l,r,ret);
	return ret;
}

template<typename t1, typename t2> 
    t_Vec<typename matrix::TypeDeduce<t1,t2>::type> 
	operator*(const t_Matrix<t1>& l, const t_Vec<t2>& r); 

/************************************************************************/
/* 3D vector                                                            */
/************************************************************************/


template<typename T> class t_Vec3 : public t_Vec<T>{
public:
	t_Vec3& set(T x, T y, T z){		
		this->operator[](0) = x;
		this->operator[](1) = y;
		this->operator[](2) = z;
		return *this;
	}
	t_Vec3<T>(T x=0.0, T y=0.0, T z=0.0):t_Vec<T>(3){set(x,y,z);};
	t_Vec3<T>(const t_Vec<T>& r): t_Vec<T>(r){
#ifdef _DEBUG
		if (r.size()!=3) 
			ssuTHROW(matrix::t_Not3D, 
			_T("3D vector error: init by non-3D general vector"));
#endif
	};
	t_Vec3& operator=(const t_Vec<T>& r){
#ifdef _DEBUG
		if (r.size()!=3) 
			ssuTHROW(matrix::t_Not3D, 
			_T("3D vector error: init by non-3D general vector"));
#endif
		(t_Vec<T>&)(*this) = r;
		return *this;
	};
};

namespace vector{
	template<typename t1, typename t2, typename TT> inline
	t_Vec3<TT> __cross(const t_Vec3<t1>& l, const t_Vec3<t2>& r){
			//typedef matrix::TypeDeduce<t1,t2>::type type;
			TT v0 = l[1]*r[2] - l[2]*r[1];
			TT v1 = l[2]*r[0] - l[0]*r[2];
			TT v2 = l[0]*r[1] - l[1]*r[0];
			return t_Vec3<TT>(v0,v1,v2);
	};
	
	template<typename t1, typename t2> IMPEXP_SMALLMAT
	    t_Vec3<typename matrix::TypeDeduce<t1,t2>::type> 
		cross(const t_Vec3<t1>& l, const t_Vec3<t2>& r);
}

//3d-vector-specific operators
template<typename t1, typename t2, typename TT> inline
t_Vec3<TT> __numbyvec3_mul(const t1 num, const t_Vec3<t2>& vec){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	t_Vec3<TT> ret;
	matrix::base::mul<t1,t2>(num, vec, ret);
	return ret;
}

template<typename t1, typename t2>
    t_Vec3<typename matrix::TypeDeduce<t1,t2>::type> 
	operator*(const t1 num, const t_Vec3<t2>& vec);

template<typename t1, typename t2, typename TT> inline
t_Vec3<TT> __vec3bynum_div(const t_Vec3<t1>& vec, const t2 num){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	t_Vec3<TT> ret;
	matrix::base::div<t1,t2>(vec, num, ret);
	return ret;
}

template<typename t1, typename t2, typename TT> inline
t_Vec3<TT> __vec3_plus(const t_Vec3<t1>& l, const t_Vec3<t2>& r){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	t_Vec3<TT> ret;
	matrix::base::plus<t1,t2>(l,r, ret);
	return ret;
}

template<typename t1, typename t2>
    t_Vec3<typename matrix::TypeDeduce<t1,t2>::type> 
	operator+(const t_Vec3<t1>& l, const t_Vec3<t2>& r);

template<typename t1, typename t2, typename TT> inline
t_Vec3<TT> __vec3_minus(const t_Vec3<t1>& l, const t_Vec3<t2>& r){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	t_Vec3<TT> ret;
	matrix::base::minus<t1,t2>(l,r, ret);
	return ret;
}

template<typename t1, typename t2>
    t_Vec3<typename matrix::TypeDeduce<t1,t2>::type> 
	operator-(const t_Vec3<t1>& l, const t_Vec3<t2>& r);

// useful geometry transformations of vectors

// get cartesian vector coords in cylindrical ref frame (x, phi, r)
namespace smat{

	IMPEXP_SMALLMAT void vec_cart_to_cyl(t_Vec3<double>& vec);
	IMPEXP_SMALLMAT void vec_cart_to_cyl(t_Vec3<t_Complex>& vec);

	// get cartesian vector coords in conical ref frame  (l, phi, h)
	IMPEXP_SMALLMAT void vec_cart_to_cone(t_Vec3<double>& vec, double semi_cone_angle);
	IMPEXP_SMALLMAT void vec_cart_to_cone(t_Vec3<t_Complex>& vec, double semi_cone_angle);

}
/************************************************************************/
/* Square Matrix                                                        */
/************************************************************************/
template<typename T>class  t_SqMatrix: public t_Matrix<T>{
	void _set_minor(t_SqMatrix& minor, int row, int col) const;
public:
	//t_SqMatrix():t_Matrix(0,0){};
	t_SqMatrix(const int dim):t_Matrix<T>(dim, dim){};;
	t_SqMatrix(const t_Matrix<T>&);
	int size() const{return this->_ncols;};
	//void resize(const int n){_cont.resize(n, t_Col(n, 0.0));};
	// construct a minor from matrix
	//t_SqMatrix(const t_SqMatrix& mat, const int raw, const int col);
	//instead use this
	// TODO: better col first
	t_SqMatrix<T> get_minor(const int raw, const int col) const;
	void setToUnity();
	void transpose();
	void setToHermConj();
	//void resize(int a_size);
	T det() const;
	t_SqMatrix<T> inverse() const;
};

template<typename T> t_SqMatrix<T>::t_SqMatrix(const t_Matrix<T>& m):t_Matrix<T>(m){
#ifdef _DEBUG
	if (m.nCols()!=m.nRows())		
		ssuTHROW(matrix::t_NotSquare, 
		_T("Square matrix error: init by non-square matrix"));
#endif
	//(t_Matrix<T>&)*this = m;
}

template<typename T> void t_SqMatrix<T>::_set_minor(t_SqMatrix<T>& ret, int row, int col) const{
#ifdef _DEBUG
	if (ret.size()!=(size()-1)) 
		ssuTHROW(matrix::t_SizeMismatch, 
		_T("Matrix set minor error: wrong minor size"));
#endif
	for (int i=0; i<col; i++){
		for (int j=0; j<row; j++){
			ret[i][j] = this->_cont[i][j];
		};
		for (int j=row+1; j<this->nRows(); j++){
			ret[i][j-1] = this->_cont[i][j];
		};
	};
	for (int i=col+1; i<this->nCols(); i++){
		for (int j=0; j<row; j++){
			ret[i-1][j] =this->_cont[i][j];
		};
		for (int j=row+1; j<this->nRows(); j++){
			ret[i-1][j-1] = this->_cont[i][j];
		};
	};
}

template<typename T> t_SqMatrix<T> t_SqMatrix<T>::get_minor
(const int row, const int col) const{
#ifdef _DEBUG
	_chk_col_ind(col);
	_chk_row_ind(row);
#endif
	t_SqMatrix<T> ret(size()-1);
	for (int i=0; i<col; i++){
		for (int j=0; j<row; j++){
			ret[i][j] = this->_cont[i][j];
		};
		for (int j=row+1; j<this->nRows(); j++){
			ret[i][j-1] = this->_cont[i][j];
		};
	};
	for (int i=col+1; i<this->nCols(); i++){
		for (int j=0; j<row; j++){
			ret[i-1][j] = this->_cont[i][j];
		};
		for (int j=row+1; j<this->nRows(); j++){
			ret[i-1][j-1] = this->_cont[i][j];
		};
	};
	return ret;
};

// calculate determinant by first row (0):
// det = SUM(j)((-1)^(0+j)*det(d(0,j)))
// d(0,j) - algebraic minor
// TODO: rewrite with gauss (QR)!!!
template<typename T> T t_SqMatrix<T>::det() const{

	if (this->_ncols==1){
		return this->_cont[0][0];
	};
	// optimize for small matrices 2x2, 3x3
	if (this->_ncols==2){
		return this->_cont[0][0]*this->_cont[1][1] - 
			this->_cont[0][1]*this->_cont[1][0];
	};
	if (this->_ncols==3){
		return 
			this->_cont[0][0]*
			(
			this->_cont[1][1]*this->_cont[2][2] - 
			this->_cont[2][1]*this->_cont[1][2]
		) - 
			this->_cont[1][0]*
			(
			this->_cont[0][1]*this->_cont[2][2] - 
			this->_cont[2][1]*this->_cont[0][2]
		) +
			this->_cont[2][0]*
			(
			this->_cont[0][1]*this->_cont[1][2]-
			this->_cont[1][1]*this->_cont[0][2]
		);
	};
	T res(0.0);
	
	t_SqMatrix<T> minor_d(this->_ncols-1); 
	for (int i=0; i<this->nCols(); i++){
		this->_set_minor(minor_d, 0, i);
		//minor = get_minor(0, i);
		T val = this->_cont[i][0];
		res+= pow(-1.0, i)*val*minor_d.det();
	}
	return res;
}

template<typename T> void t_SqMatrix<T>::setToUnity(){
	for (int i=0; i<this->size(); i++){
		for (int j=0; j<this->size(); j++){
			if (i==j) 
				this->_cont[i][j]=1.0;
			else 
				this->_cont[i][j]=0.0;
		};
	};
};

template<typename T> void t_SqMatrix<T>::transpose(){
	T swp_val;
	for (int i=0; i<this->_ncols; i++)
		for (int j=0; j<=i; j++){
			swp_val = this->_cont[i][j];
			this->_cont[i][j] = this->_cont[j][i];
			this->_cont[j][i] = swp_val;
		}
};

// TODO: rewrite with robust procedures!!!
template<typename T> t_SqMatrix<T> t_SqMatrix<T>::inverse() const{
	int dim = this->nCols();
	t_SqMatrix<T> ret(dim);
	t_SqMatrix<T> unity(dim);
	t_SqMatrix<T> temp(dim);
	t_Vec<T> buf_vec(dim);
	unity.setToUnity();
	T coef = 1.0/this->det();
	for (int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			temp = *this;
			//temp[j] = unity[i];
			unity.col_to_vec(i, buf_vec);
			temp.set_col(j, buf_vec);
			ret[i][j] = coef*temp.det();
		}
	}
	return ret;
};

//sq-matrix-specific operators

template<typename t1, typename t2, typename TT> inline
t_SqMatrix<TT> __sqmatbynum_mul(const t_SqMatrix<t1>& mat, const t2 num){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMatrix<TT> ret(mat.size());
	matrix::base::mul<t1,t2>(num, mat, ret);
	return ret;
}

template<typename t1, typename t2, typename TT> inline
t_SqMatrix<TT> __numbysqmat_mul(const t2 num, const t_SqMatrix<t1>& mat){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMatrix<TT> ret(mat.size());
	matrix::base::mul<t1,t2>(num, mat, ret);
	return ret;
}

template<typename t1, typename t2, typename TT> inline
t_SqMatrix<TT> __sqmat_mul(const t_SqMatrix<t1>& l, const t_SqMatrix<t2>& r){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMatrix<TT> ret(l.size());
	matrix::base::mat_mul<t1,t2>(l, r, ret);
	return ret;
}

template<typename t1, typename t2, typename TT> inline
t_SqMatrix<TT> _sqmatbynum_div(const t_SqMatrix<t1>& l, const t2 r){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMatrix<TT> ret(matrix::base::div<t1,t2>(l,r));
	return ret;
}

template<typename t1, typename t2, typename TT> inline
t_SqMatrix<TT> __sqmat_plus(const t_SqMatrix<t1>& l, const t_SqMatrix<t2>& r){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMatrix<TT> ret(l.size());
	matrix::base::plus<t1,t2>(l,r, ret)
	return ret;
}

template<typename t1, typename t2, typename TT> inline
t_SqMatrix<TT > __sqmat_minus(const t_SqMatrix<t1>& l, const t_SqMatrix<t2>& r){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMatrix<TT> ret(l.size());
	matrix::base::minus<t1,t2>(l,r, ret);
	return ret;
}
/************************************************************************/
/* 3D matrix (orthogonal transformation matrix)                         */
/************************************************************************/

template<typename T> class t_SqMat3 : public t_SqMatrix<T>{
public:
	enum t_Init{Zero=0, Unity};
	t_SqMat3(t_Init init_type = Zero):t_SqMatrix<T>(3){
		if (init_type==Unity){
			this->_cont[0][0]=1.0;
			this->_cont[1][1]=1.0;
			this->_cont[2][2]=1.0;
		};
	};

	t_SqMat3 inverse() const{
		t_SqMat3 ret;
		const t_SqMatrix<T>& rM = *this;
		(t_SqMatrix<T>&)ret = rM.inverse();
		return ret;
	};

};

// 3d-matrix-specific operators

template<typename t1, typename t2, typename TT> inline
t_SqMat3<TT> __sqmat3bynum_mul(const t_SqMat3<t1>& mat, const t2 num){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMat3<TT> ret;
	matrix::base::mul<t1,t2>(num, mat, ret);
	return ret;
}

template<typename t1, typename t2, typename TT> inline
t_SqMat3<TT > __numbysqmat3_mul(const t1 l, const t_SqMat3<t2>& r){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMat3<TT> ret;
	matrix::base::mul<t1,t2>(l,r, ret);
	return ret;
}

template<typename t1, typename t2, typename TT> inline
t_Vec3<TT> __sqmat3byvec3_mul(const t_SqMat3<t1>& l, const t_Vec3<t2>& r){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	t_Vec3<TT> ret;
	matrix::base::mat_mul<t1,t2>(l,r, ret);
	return ret;
}

template<typename t1, typename t2, typename TT> inline
t_SqMat3<TT> __sqmat3_mul(const t_SqMat3<t1>& l, const t_SqMat3<t1>& r){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMat3<type> ret;
	matrix::base::mat_mul<t1,t2>(l,r, ret);
	return ret;
}

template<typename t1, typename t2, typename TT> inline
t_SqMat3<TT> __sqmat3bynum_div(const t_SqMat3<t1>& l, const t2 r){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMat3<TT> ret;
	matrix::base::div<t1,t2>(l,r, ret);
	return ret;
}

template<typename t1, typename t2, typename TT> inline
t_SqMat3<TT> __sqmat3_plus(const t_SqMat3<t1>& l, const t_SqMat3<t2>& r){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMat3<TT> ret;
	matrix::base::plus<t1,t2>(l,r, ret);
	return ret;
}

template<typename t1, typename t2, typename TT> inline
t_SqMat3<TT> __sqmat3_minus(const t_SqMat3<t1>& l, const t_SqMat3<t2>& r){
	//typedef matrix::TypeDeduce<t1,t2>::type type;
	t_SqMat3<TT> ret;
	matrix::base::minus<t1,t2>(l,r,ret);
	return ret;
}

typedef t_Matrix<t_Complex> t_MatCmplx;
typedef t_Matrix<double> t_MatDbl;

typedef t_Vec<t_Complex> t_VecCmplx;
typedef t_Vec<double> t_VecDbl;

typedef t_Vec3<t_Complex> t_Vec3Cmplx;
typedef t_Vec3<double> t_Vec3Dbl;

typedef t_SqMatrix<t_Complex> t_SqMatCmplx;
typedef t_SqMatrix<double> t_SqMatDbl;

typedef t_SqMat3<t_Complex> t_SqMat3Cmplx;
typedef t_SqMat3<double> t_SqMat3Dbl;

#endif // __MATH_OPERANDS
