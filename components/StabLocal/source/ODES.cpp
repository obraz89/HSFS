#include "ODES.h"
#include "log.h"
#include <iostream>
// default 3D context

const int t_ODES::_dim = 4;

t_ODES::t_ODES():
_nnodes(), varRange(), 
solution(0, t_MatCmplx(_dim, 2*_dim)){};

void t_ODES::clean(){
		_orthStack.clear();
};

void t_ODES::setContext(const int& a_newNnodes){
	clean();
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

	return h+1.0/6.0*step*(k1+2.0*k2+2.0*k3+k4);
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
		for (int j=0; j<_dim;j++){
			solution[i+1][j] = stepRK3D(i, solution[i][j]);
		}

	}
}

t_Complex t_ODES::detGS(const t_MatCmplx& sol, const int& rank) const{

	if (rank>sol.nCols()) Log<<_T("GS determinant error: Rank exceeds task dimension\n");
	if (rank<0) Log<<_T("GS determinant error: Rank < 0 \n");
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
	if (abs(det.imag()/det.real())>1.0e-6){
		Log<<_T("GS determinant error: determinant is complex; imag:")
			<<abs(det.imag())<<_T("\n");
	}
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
	// ugly reconstruction from upper boundary down to wall
	// this doesn't work for sure
	/*
	std::vector<t_MatCmplx> ret(this->solution);
	t_SqMatCmplx trans(_dim);
	t_SqMatCmplx direct(_dim);
	int orthStackInd = 0;
	direct.setToUnity();
	trans.setToUnity();
	int orthStackSize = _orthStack.size();
	for (int i=0; i<_nnodes; i++){
		if (orthStackInd<orthStackSize){
			if (i==this->_orthStack[orthStackInd].ind){
				direct = direct*(_orthStack[orthStackInd].orthMatrix);
				trans = direct.inverse();		
				orthStackInd++;
			}
		}
		ret[i] = solution[i]*(trans);
	}
	return ret;
	*/

	// reconstruct from wall to outer region
	std::vector<t_MatCmplx> ret(_nnodes, t_MatCmplx(_dim, 2*_dim));
	int orthStackInd = _orthStack.size()-1;

	t_SqMatCmplx direct(_dim);
	direct.setToUnity();

	for (int i=_nnodes-1; i>=0; i--){
		ret[i] = solution[i]*direct;
		if (orthStackInd>=0){
			if (i==_orthStack[orthStackInd].ind){
				// TODO: WTF??? What's the right order?
				direct = (_orthStack[orthStackInd].orthMatrix)*direct;
				orthStackInd--;
			}
		}
	}
	return ret;
}