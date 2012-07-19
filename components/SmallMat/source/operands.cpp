#include "math_operands.h"

/************************************************************************/
/*                                                                t_Vec */
/************************************************************************/

void t_Vec::_chk_size_match(const t_Vec& l, const t_Vec& r){
	if (l.size()==r.size()) return;
	ssuTHROW(t_SizeMismatch, _T("Vector binary operand error: size mismatch"));
};

/************************************************************************/
/*                                                             t_Matrix */
/************************************************************************/

// zero matrix with correct sizes
t_Matrix::t_Matrix(const int a_nVecs, const int a_nElemInVec)
:_cont(a_nVecs, t_Vec(a_nElemInVec, 0.0)){};

void t_Matrix::_chk_col_ind(const int n) const{
	if ((n<0)||(n>=nCols())) 
		ssuTHROW(t_ColOutOFRange, _T("Matrix Error: Column index out of range"));
};

void t_Matrix::_chk_row_ind(const int m) const{
	if ((m<0)||(m>=nRows())) 
		ssuTHROW(t_RowOutOFRange, _T("Matrix Error: Row index out of range"));
};

void t_Matrix::_chk_size_match(const t_Matrix& l, const t_Matrix& r){
	if (l.nCols()!=r.nRows())
		ssuTHROW(t_SizeMismatch, _T("Matrix Error: Binary operation size mismatch"));
};

const t_Vec& t_Matrix::operator [](const int n) const{
	_chk_col_ind(n);
	return _cont[n];
};

t_Vec& t_Matrix::operator [](const int n){
	_chk_col_ind(n);
	return _cont[n];
};

t_Matrix operator*(const t_Matrix & l, const t_Matrix& r){
	t_Matrix::_chk_size_match(l, r);
	t_Matrix ret(r.nCols(), l.nRows());
	for (int i=0; i<r.nCols(); i++){
		for (int j=0; j<l.nCols(); j++){
			ret[i] = ret[i] + r[i][j]*l[j];
		};
	};
	return ret;
};


/************************************************************************/
/*                                                           t_SqMatrix */
/************************************************************************/

t_SqMatrix::t_SqMatrix(const int dim):t_Matrix(dim, dim){};
t_SqMatrix t_SqMatrix::get_minor(const int row, const int col) const{
	_chk_col_ind(col);
	_chk_row_ind(row);
	t_SqMatrix ret(size()-1);
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

t_Complex t_SqMatrix::det() const{
	// recursive procedure
	// calculate determinant by first row (0):
	// det = SUM(j)((-1)^(0+j)*det(d(0,j)))
	// d(0,j) - algebraic minor
	// TODO: rewrite with gauss (QR)
	if (this->nCols()==1){
		return _cont[0][0];
	}; 
	t_Complex res=0.0;
	for (int i=0; i<nCols(); i++){
		t_SqMatrix minor = get_minor(0, i);
		t_Complex val = _cont[i][0];
		res+= pow(-1.0, i)*val*minor.det();
	}
	return res;
}

void t_SqMatrix::setToUnity(){
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
t_SqMatrix t_SqMatrix::inverse() const{
	int dim = this->nCols();
	t_SqMatrix ret(dim);
	t_SqMatrix unity(dim);
	unity.setToUnity();
	t_Complex coef = 1.0/this->det();
	for (int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			t_SqMatrix temp = *this;
			temp[j] = unity[i];
			ret[i][j] = coef*temp.det();
		}
	}
	return ret;
}

t_SqMatrix operator*(const t_SqMatrix & l, const t_SqMatrix& r){
	t_SqMatrix ret;
	const t_Matrix &rml = l;
	const t_Matrix &rmr = r;
	t_Matrix &rret = ret;
	rret = rml*rmr;
	return ret;
};
