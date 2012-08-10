#include "ODES.h"
#include <iostream>
// default 3D context
t_ODES::t_ODES():
_dim(4), _nnodes(), varRange(), 
solution(0, t_MatCmplx(_dim, _dim)){};

void t_ODES::set3DContext(){
		_dim = 4;
		solution.clear();
		varRange.clear();
		_orthStack.clear();
};

void t_ODES::resizeGrid(const int& a_newNnodes){
	_nnodes = a_newNnodes;
	varRange.resize(a_newNnodes);
	solution.resize(a_newNnodes, t_MatCmplx(_dim, 2*_dim));
};

// orthoStack
	t_ODES::t_OrthPoint::t_OrthPoint(const int &a_ind, const int &a_dim)
		:ind(a_ind), orthMatrix(a_dim){};

t_VecCmplx t_ODES::stepRK3D(const int& ind, const t_VecCmplx& h){
	double x = this->varRange[ind];
	double step = this->varRange[ind+1] - x;
	t_VecCmplx k1 = this->formRHS3D(x, h);
	t_VecCmplx k2 = this->formRHS3D(x + 0.5*step, h+0.5*step*k1);
	t_VecCmplx k3 = this->formRHS3D(x + 0.5*step, h+0.5*step*k2);
	t_VecCmplx k4 = this->formRHS3D(x + step, h+step*k3);
	// at x+dx
	return (h+1.0/6.0*step*(k1+2.0*k2+2.0*k3+k4));
};

void t_ODES::solve(){
	for (int i=0; i<this->_nnodes-1; i++){
		// ortho ...
		if (this->needOrtho(solution[i])){
			this->ortho(i);
			// ortho solution
			this->solution[i] = this->solution[i]*(this->_orthStack.back().orthMatrix);
		};
		// get solution
		for (int j=0; j<this->_dim;j++){
			solution[i+1][j] = stepRK3D(i, solution[i][j]);
		}

	}
}

t_Complex t_ODES::detGS(const t_MatCmplx& sol, const int& rank) const{
	// exceptions 
	if (rank>sol.nCols()) std::cerr<<"GS determinant error: Rank exceeds task dimension\n";
	if (rank<0) std::cerr<<"GS determinant error: Rank < 0 \n";
	if (rank==0) return 1.0;
	t_SqMatCmplx mat(rank);
	for (int i=0; i<rank; i++){
		for (int j=0; j<=i; j++){
			t_Complex val = 
				vector::dot(t_VecCmplx(sol[i]), t_VecCmplx(sol[j]));
			// Ermith Matrix
			mat[i][j]= val;
			mat[j][i]= std::conj(val);	
		}
	}
	t_Complex det = mat.det();
	if (abs(det.imag()/det.real())>1.0e-6) std::cerr<<"GS determinant error: determinant is complex; imag:"
										 <<abs(det.imag())<<std::endl;
	return mat.det();
}

t_Complex t_ODES::minorGS(const t_MatCmplx& vectors, const int& dim, const int& nExcludeCol) const{
	if (dim==0) return 1.0;
	t_SqMatCmplx ret(dim);
	for (int i=0; i<nExcludeCol; i++){
		for(int j=0; j<dim; j++){
			ret[i][j] = 
				vector::dot(t_VecCmplx(vectors[i]), t_VecCmplx(vectors[j]));
		}
	}
	for (int i=nExcludeCol+1; i<=dim; i++){
		for(int j=0; j<dim; j++){
			ret[i-1][j] =
				vector::dot(t_VecCmplx(vectors[i]), t_VecCmplx(vectors[j]));
		}
	}
	return ret.det();
}

void t_ODES::ortho(const int& nnode){
	// matrix S  - transition from basis 
	// e(A,B,C...)  to a new one orthonormal e'(A~, B~, C~,...)
	// e'=e*S
	// via Gramm-Shmidt orthonormization
	// stored by columns
	t_SqMatCmplx trans_matrix(_dim);
	for (int i=0; i<_dim; i++){
		t_VecCmplx column(_dim, 0.0);
		// fill non-zero elements
		for (int j=0; j<=i; j++){
			double dj0 = detGS(solution[nnode], i).real();
			double dj1 = detGS(solution[nnode], i+1).real();
			double coef = 1.0/sqrt(dj0*dj1);
			column[j] = pow(-1.0, i+j)*coef*minorGS(solution[nnode], i, j);
		}
		trans_matrix[i] = column;
	}
	t_ODES::t_OrthPoint orth_point(nnode, _dim);
	orth_point.orthMatrix = trans_matrix;
	orth_point.ind = nnode;
	this->_orthStack.push_back(orth_point);

};

std::vector<t_MatCmplx> t_ODES::reconstruct(){
	std::vector<t_MatCmplx> ret(this->solution);
	t_SqMatCmplx trans(_dim);
	t_SqMatCmplx direct(_dim);
	int orthStackInd = 0;
	direct.setToUnity();
	trans.setToUnity();
	for (int i=0; i<_nnodes-1; i++){
		if (i==this->_orthStack[orthStackInd].ind){
			direct = direct*(_orthStack[i].orthMatrix);
			trans = direct.inverse();		
			orthStackInd++;
		}
			ret[i] = solution[i]*(trans);
	}

	// reconstruct last record (for which we don't do ortho)
	ret[_nnodes-1] = solution[_nnodes-1]*(trans);
	return ret;
}