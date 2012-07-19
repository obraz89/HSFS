#ifndef __MATH_OPERANDS
#define __MATH_OPERANDS
#include <complex>
#include <vector>
#include <cmath>
#include "stdafx.h"
#include "gen_exception.h"

// TO REMOVE
#include <boost/numeric/mtl/mtl.hpp>
// ~TO REMOVE

typedef std::complex<double> t_Complex;
typedef std::complex<double> t_CompVal;
typedef std::vector<double> t_DblVec;


// TO REMOVE
using namespace mtl;
typedef vector::parameters<tag::col_major, vector::fixed::dimension<3> > _param_fixed_vec3;
typedef matrix::parameters<tag::col_major, mtl::index::c_index, mtl::fixed::dimensions<3, 3> > _param_fixed_mat33;
typedef  mtl::dense2D<double, _param_fixed_mat33> t_SqMat3;
typedef  mtl::dense_vector<double, _param_fixed_vec3> t_Vec3;
typedef mtl::dense_vector<t_CompVal, _param_fixed_vec3> t_CompVec3;

// ~TO REMOVE

/************************************************************************/
/* Linear Algebra Vector class
** with a second norm as norm
** the Vector also is used to form
** a matrix by columns (see t_Matrix)
*/
/************************************************************************/

class  t_Vec{
protected:
	std::vector< t_Complex > _cont;
	typedef std::vector<double>::iterator Iter;
	typedef std::vector<double>::const_iterator Citer;

	inline static void _chk_size_match(const t_Vec& l, const t_Vec& r);
	// no default constructors 
public:
	t_Vec(const int dim, const t_Complex val=0.0):_cont(dim,val){};
	int size() const{return _cont.size();};
	// TODO: maybe some check that imag is small ?
	double norm() const{return dot(*this, *this).real();}; 
	// TODO: check index? hmm... if something is wrong fail anyway)
	t_Complex& operator[](const int ind){return _cont[ind];};
	const t_Complex& operator[](const int ind) const{return _cont[ind];};

	friend t_Vec operator*(const t_Vec& v, const t_Complex c){
		t_Vec ret(v);
		for (int i=0; i<v.size(); i++){
			ret[i]=v[i]*c;
		}
		return ret;
	};
	friend t_Vec operator*(const t_Complex c, const t_Vec& v){
		return v*c;
	};
	friend t_Vec operator/(const t_Vec& v, const t_Complex c){
		return v*(1.0/c);
	};
	friend t_Vec operator+(const t_Vec& l, const t_Vec& r){
		_chk_size_match(l,r);
		t_Vec ret(l);
		for (int i=0; i<l.size(); i++){
			ret[i]+=r[i];
		}
		return ret;
	};
	friend t_Vec operator-(const t_Vec& l, const t_Vec& r){
		// to avoid copy overhead we can't do
		// l+r*(-1)
		_chk_size_match(l,r);
		t_Vec ret(l);
		for (int i=0; i<l.size(); i++){
			ret[i]-=r[i];
		}
		return ret;
	};

	static t_Complex dot(const t_Vec& l, const t_Vec& r) {
		t_Complex ret = 0.0;
		_chk_size_match(l, r);
		for (int i=0; i<l.size(); i++){
			ret+=l[i]*std::conj(r[i]);
		}
		return ret;
	};

	//exceptions
	class t_SizeMismatch: public t_GenException{
	public:
		t_SizeMismatch(const wxString& what, const wxChar* szFile,  const int line):
			t_GenException(what, szFile, line){};
	};
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
class  t_Matrix{
protected:
	std::vector<t_Vec> _cont;
	static void _chk_size_match(const t_Matrix& l, const t_Matrix& r); 
	void _chk_col_ind(const int n) const;
	void _chk_row_ind(const int n) const; 
public:
	int nCols() const{return _cont.size();};
	int nRows() const{if (nCols()>0) {return _cont[0].size();} else {return 0;}; };
	t_Vec& operator[](const int n);
	const t_Vec& t_Matrix::operator [](const int n) const;
	t_Matrix(const int a_nVecs, const int a_nElemInVec);
	friend t_Matrix operator*(const t_Matrix& l, const t_Matrix& r);

	//exceptions
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
};

/************************************************************************/
/* Square Matrix                                                        */
/************************************************************************/
class  t_SqMatrix: public t_Matrix{
public:
	t_SqMatrix():t_Matrix(0,0){};
	t_SqMatrix(const int dim);
	int size() const{return nCols();};
	void resize(const int n){_cont.resize(n, t_Vec(n, 0.0));};
	// construct a minor from matrix
	//t_SqMatrix(const t_SqMatrix& mat, const int raw, const int col);
	//instead use this
	// TODO: better col first
	t_SqMatrix get_minor(const int raw, const int col) const;
	void setToUnity();
	//void resize(int a_size);
	t_Complex det() const;
	t_SqMatrix inverse() const;

	// exceptions
	class t_NotSquare: public t_GenException{
	public:
		t_NotSquare(const wxString& what, const wxChar* szFile,  const int line):
		  t_GenException(what, szFile, line){};
	};
};

t_SqMatrix operator*(const t_SqMatrix& l, const t_SqMatrix& r);


// Cartesian vector class
// {i,j,k} <-> {x, y, z}
// To replace MTL code
/*
class t_Vec3: public t_Vec{
public:
	t_Vec3(const t_Complex val=0.0):t_Vec(3, val){};
	t_Vec3 cross(const t_Vec3& r) const{
		t_Vec3 ret;
		// use sq matrix ?
		//ret[0] = 
		return ret;
	};

};
*/

#endif // __MATH_OPERANDS