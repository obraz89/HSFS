#include <vector>
#include <complex>
// all computations with double precision
#ifndef __ODES_OPERANDS
#define __ODES_OPERANDS
typedef std::complex<double> t_Complex;
class t_Vec;
//typedef t_Vec (*pFunRHS)(const double&, const t_Vec&);
class t_Vec{
	std::vector< t_Complex > _cont;
// no default constructors 
public:
	t_Vec(const int& dim, const t_Complex val);
	~t_Vec();
	double norm() const;
	int size() const;
	t_Complex& operator[](const int& ind);
	const t_Complex& operator[](const int& ind) const;
	t_Vec operator*(const t_Complex&) const;
	t_Vec operator/(const t_Complex&) const;
	t_Vec operator+(const t_Vec&) const;
	t_Vec operator-(const t_Vec&) const;
	t_Complex scalProd(const t_Vec&) const;	
};

// arrrgghhh;
t_Vec operator*(const t_Complex&, const t_Vec&);

// matrix class
// operations are column-based
// e.g. matrix of solution contains at a given point
// _dim columns - vectors of solution
// inverse matrix is set by columns too
// multiplication implies ordinary right mult: A.mul(B) <=> C=A*B set again by columns: C[i][j] - column j, row i
// if matrix A is mul by set of vectors  h={h1, h2,} this means h*A ( like "basis" transition)

class t_Matrix{
protected:
	std::vector<t_Vec> _cont;
public:
	int nCols, nRows;
	t_Vec& operator[](const int& n);
	const t_Vec& t_Matrix::operator [](const int& n) const;
	t_Matrix(const int& a_nVecs, const int& a_nElemInVec);
	
	t_Matrix mul(const t_Matrix&);
	~t_Matrix();
};

class t_SqMatrix: public t_Matrix{
public:
	t_SqMatrix(const int& dim);
	// construct a minor from matrix
	t_SqMatrix(const t_SqMatrix& mat, const int& raw, const int& col);
	t_SqMatrix& operator=(const t_Matrix& rval);
	~t_SqMatrix();
	void setToUnity();
	void resize(int a_size);
	t_Complex det() const;
	t_SqMatrix inverse() const;
};
#endif //__ODES_OPERANDS