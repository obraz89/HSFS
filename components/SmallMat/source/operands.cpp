#include "math_operands.h"

t_Vec::t_Vec(const int& dim, const t_Complex val=0.0):_cont(dim,val){};
t_Vec::~t_Vec(){};

t_Complex& t_Vec::operator[](const int &ind){
	// TODO: check bounds
	return this->_cont[ind];
};

int t_Vec::size() const{
	return (this->_cont.size());
};

t_Complex t_Vec::scalProd(const t_Vec& r) const{
	// check size match, throw
	t_Complex ret(0.0);
	for (int i=0; i<this->size(); i++){
		ret+=this->_cont[i]*std::conj(r[i]);
	}
	return ret;
};

double t_Vec::norm() const{
	return sqrt((this->scalProd(*this)).real());
};

const t_Complex& t_Vec::operator [](const int& ind) const{
	return this->_cont[ind];
};

t_Vec t_Vec::operator*(const t_Complex& val) const{
	t_Vec ret(*this);
	for (int i=0; i<ret.size(); i++){
		ret[i]=(*this)[i]*val;
	}
	return ret;
};

t_Vec t_Vec::operator/(const t_Complex& val) const{
	t_Vec ret(*this);
	t_Complex rev = 1.0/val;
	for (int i=0; i<ret.size(); i++){
		ret[i]=(*this)[i]*rev;
	}
	return ret;
};


t_Vec t_Vec::operator+(const t_Vec& r) const{
	// TODO: check sizes and throw size exceptions
	t_Vec ret(*this);
	for (int i=0; i<ret.size(); i++){
		ret[i] = this->_cont[i] + r[i];
	}
	return ret;
}

t_Vec t_Vec::operator-(const t_Vec& r) const{
	// TODO: check sizes and throw size exceptions
	t_Vec ret(*this);
	for (int i=0; i<ret.size(); i++){
		ret[i] = this->_cont[i] - r[i];
	}
	return ret;
}

// global operators
t_Vec operator*(const t_Complex& l, const t_Vec& r){
	t_Vec ret(r);
	for (int i=0; i<r.size(); i++){
		ret[i] = r[i]*l;
	}
	return ret;
};

// zero matrix with correct sizes
// no default constructor allowed to avoid size errors
t_Matrix::t_Matrix(const int& a_nVecs, const int& a_nElemInVec)
	:nCols(a_nVecs), nRows(a_nElemInVec),_cont(a_nVecs, t_Vec(a_nElemInVec, 0.0)){}
const t_Vec& t_Matrix::operator [](const int& n) const{
	return _cont[n];
};

t_Vec& t_Matrix::operator [](const int& n){
	return _cont[n];
};

t_Matrix t_Matrix::mul(const t_Matrix & r){
	// TODO: exceptions
	if (this->nCols!=r.nRows){
		std::cerr<<"Matrix multiplication error: sizes do not match\n";
	}; 
	t_Matrix ret(r.nCols, this->nRows);
	for (int i=0; i<r.nCols; i++){
		for (int j=0; j<this->nCols; j++){
					ret[i] = ret[i] + r[i][j]*this->_cont[j];
		};
	};
	return ret;
}
t_Matrix::~t_Matrix(){};

t_SqMatrix::t_SqMatrix():t_Matrix(0,0){};
t_SqMatrix::t_SqMatrix(const int& dim):t_Matrix(dim, dim){};
t_SqMatrix::t_SqMatrix(const t_SqMatrix& mat, const int& row, const int& col):t_Matrix(mat.nCols-1, mat.nCols-1){
	t_Matrix& self = *this;
	for (int i=0; i<col; i++){
		for (int j=0; j<row; j++){
			self[i][j] = mat[i][j];
		};
		for (int j=row+1; j<mat.nRows; j++){
			self[i][j-1] = mat[i][j];
		};
	};
	for (int i=col+1; i<mat.nCols; i++){
		for (int j=0; j<row; j++){
			self[i-1][j] = mat[i][j];
		};
		for (int j=row+1; j<mat.nRows; j++){
			self[i-1][j-1] = mat[i][j];
		};
	}
};
t_SqMatrix::~t_SqMatrix(){};


t_Complex t_SqMatrix::det() const{
	// recursive procedure
	// calculate determinant by first row (0):
	// det = SUM(j)((-1)^(0+j)*det(d(0,j)))
	// d(0,j) - algebraic minor
	// TODO: rewrite with gauss (QR)
	if (this->nCols==1){
		//const t_Matrix& self = *this;
		//return self[0][0]*self[1][1] - self[0][1]*self[1][0];
		return (*this)[0][0];
	} 
	t_Complex res=0.0;
	for (int i=0; i<nCols; i++){
		t_SqMatrix minor(*this, 0, i);
		t_Complex val = (*this)[i][0];
		res+= pow(-1.0, i)*val*minor.det();
	}
	return res;
}

void t_SqMatrix::setToUnity(){
	for (int i=0; i<this->nCols; i++){
		for (int j=0; j<this->nCols; j++){
			if (i==j) 
				(*this)[i][j]=1.0;
			else 
				(*this)[i][j]=0.0;
		};
	};
};
// TODO: rewrite with robust procedures
// e.g. Gauss
t_SqMatrix t_SqMatrix::inverse() const{
	int dim = this->nCols;
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

// operators
t_SqMatrix& t_SqMatrix::operator=(const t_Matrix& rval){
	//TODO: size exception
	if (rval.nCols!=rval.nRows) std::cerr<<"Matrix error: assigning non-square matrix to square one"<<std::endl;
	t_Matrix& base = *this;
	base = rval;
	return *this;
}

void t_SqMatrix::resize(int a_size){
	_cont.resize(a_size, t_Vec(a_size, 0.0));
};
