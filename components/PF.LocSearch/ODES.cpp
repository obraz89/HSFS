#include "stdafx.h"

#include "ODES.h"
#include "log.h"
#include <iostream>
// default 3D context


static const int NNODES_BIG = 5001;

t_ODES::t_ODES(const int NSolVecs, const int SolVecDim):
_NSolVecs(NSolVecs), _SolVecDim(SolVecDim),
_nnodes(0), varRange(NNODES_BIG), 
solution(NNODES_BIG, t_MatCmplx(NSolVecs, SolVecDim)), 
_orthStack(NNODES_BIG, t_OrthPoint(-1, NSolVecs)),
_orth_stack_len(0),
_rk1(SolVecDim), _rk2(SolVecDim), _rk3(SolVecDim), _rk4(SolVecDim),
_hl(SolVecDim), _hr(SolVecDim), _h_res(SolVecDim),
_sol_vec_tmp1(SolVecDim), _sol_vec_tmp2(SolVecDim),
_sol_tmp(NSolVecs, SolVecDim),
_trans_mat_vec_tmp1(NSolVecs), _trans_mat_vec_tmp2(NSolVecs){};

void t_ODES::clean(){
		//_nnodes = 0;
		_orth_stack_len = 0;
		// do not clear orth stack for speed up
		//_orthStack.clear();
};

void t_ODES::setContext(const int& a_newNnodes){
	clean();
	if (a_newNnodes>=NNODES_BIG){
		wxString msg(_T("Too many points in ODES discretization, aborting"));
		wxLogError(msg);
		ssuGENTHROW(msg);
	}
	_nnodes = a_newNnodes;
	//varRange.resize(a_newNnodes);
	//solution.resize(a_newNnodes, t_MatCmplx(_dim, 2*_dim));
};

// orthoStack
	t_ODES::t_OrthPoint::t_OrthPoint(const int &a_ind, const int &a_dim)
		:ind(a_ind), orthMatrix(a_dim){};

// stepRK - Runge-Kutta step for a particular vector, returns solution vector
void t_ODES::stepRK3D(const int& ind, const t_VecCmplx& h, t_VecCmplx& dest){
	double x = varRange[ind];
	double step = varRange[ind+1] - x;
	//t_VecCmplx k1 = this->formRHS3D(x, h);
	//t_VecCmplx k2 = this->formRHS3D(x + 0.5*step, h+0.5*step*k1);
	//t_VecCmplx k3 = this->formRHS3D(x + 0.5*step, h+0.5*step*k2);
	//t_VecCmplx k4 = this->formRHS3D(x + step, h+step*k3);
	//return h+1.0/6.0*step*(k1+2.0*k2+2.0*k3+k4);

	formRHS3D(x, h, _rk1);

	for (int i = 0; i < _SolVecDim; i++)
		_h_res[i] = h[i] + 0.5*step*_rk1[i];	

	formRHS3D(x + 0.5*step, _h_res, _rk2);

	for (int i = 0; i < _SolVecDim; i++)
		_h_res[i] = h[i] + 0.5*step*_rk2[i];

	formRHS3D(x + 0.5*step, _h_res, _rk3);

	for (int i = 0; i < _SolVecDim; i++)
		_h_res[i] = h[i] + step*_rk3[i];

	formRHS3D(x + step, _h_res, _rk4);

	double c = 1.0 / 6.0*step;

	for (int i = 0; i < _SolVecDim; i++)
		dest[i] = h[i] + c*(_rk1[i] + 2.0*_rk2[i] + 2.0*_rk3[i] + _rk4[i]);

};

void t_ODES::solve(){

	for (int i=0; i<_nnodes-1; i++){
		// ortho ...
		if (this->needOrtho(solution[i])){
			this->ortho(i);
			// ortho solution
			//this->solution[i] = this->solution[i]*(this->_orthStack.back().orthMatrix);

			const t_SqMatCmplx& rOrthM = _orthStack[_orth_stack_len-1].orthMatrix;
			matrix::base::mat_mul(solution[i], rOrthM, _sol_tmp);
			solution[i] = _sol_tmp;
			
		};
		// get solution
		for (int j=0; j<_NSolVecs;j++){

			//solution[i].col_to_vec(j, sol_vec);
			//solution[i+1].set_col(j, stepRK3D(i, sol_vec));

			solution[i].col_to_vec(j, _sol_vec_tmp1);
			stepRK3D(i, _sol_vec_tmp1, _sol_vec_tmp2);
			solution[i+1].set_col(j, _sol_vec_tmp2);
		}

	}
}

// ortho - makes transition to new basis of solution vectors

// matrix S  - transition from basis 
// (v1, v2, v3, v4)  to a new one orthonormal (e1, e2, e3, e4)
// (e1, e2, e3, e4) = (v1, v2, v3, v4)*S
// make Gramm-Shmidt orthonormization
// S<=>trans_matrix
void t_ODES::ortho(const int& nnode){

	_orth_stack_len++;
	t_OrthPoint& orth_point = _orthStack[_orth_stack_len-1];
	orth_point.ind = nnode;

	t_SqMatCmplx& trans_matrix = orth_point.orthMatrix;

	// variant 1 - calculate orthMatrix via determinant formulas
	/*{	
		double dj0, dj1, coef;
		for (int i = 0; i<_NSolVecs; i++) {
			_trans_mat_vec_tmp1.set_vals(0.0);
			// fill non-zero elements
			for (int j = 0; j <= i; j++) {
				dj0 = detGS(solution[nnode], i).real();
				dj1 = detGS(solution[nnode], i + 1).real();
				coef = 1.0 / sqrt(dj0*dj1);
				_trans_mat_vec_tmp1[j] = pow(-1.0, i + j)*coef*minorGS(solution[nnode], i, j);
			}

			//trans_matrix[i] = column;
			trans_matrix.set_col(i, _trans_mat_vec_tmp1);
		}
		//wxLogMessage(_T("============\nTrans mat old:%s"), trans_matrix.to_wxstr());
	}*/

	// variant 2 - calculate orthMatrix via vector orthogonalization formulas
	{
		// make orthogonalization (v1,v2,v3,v4)-> (e1,e2,e3,e4)
		const t_MatCmplx& sol = solution[nnode];
		t_MatCmplx& sol_orth = _sol_tmp; 
		t_VecCmplx& vl = _sol_vec_tmp1;
		t_VecCmplx& vr = _sol_vec_tmp2;
		t_VecCmplx& vtl = _trans_mat_vec_tmp1;
		t_VecCmplx& vtr = _trans_mat_vec_tmp2;
		// orthMatrix coefficients calculated iteratively
		// column by column: 1->2->3->4
		t_Complex dotpr, coef, norm;
		for (int iCol = 0; iCol < _NSolVecs; iCol++) {
			for (int iRow = 0; iRow < _NSolVecs; iRow++) {
				coef = (iCol == iRow) ? 1.0 : 0.0;
				// ortho
				for (int j = 1; j <= iCol; j++) {
					sol.col_to_vec(iCol, vl);
					sol_orth.col_to_vec(j - 1, vr);
					dotpr = vector::dot(vl, vr);
					coef += -1.0*dotpr*trans_matrix[j-1][iRow];
				}
				trans_matrix[iCol][iRow] = coef;
			}

			// first, compute u[iCol] - orth vector that is to be normalized
			// e[iCol] = u[iCol]/norm(u[iCol])
			trans_matrix.col_to_vec(iCol, vtl);
			matrix::base::mat_mul(sol, vtl, vl);		
			norm = vl.norm();
			// multiply coefs by inverse norm of u, this is final values 
			// of the column iCol of transition matrix
			matrix::base::mul(1.0 / norm, vtl, vtr);
			trans_matrix.set_col(iCol, vtr);

			// now the correct vector e[iCol] can be computed
			// e[iCol] is stored in sol_orth matrix
			trans_matrix.col_to_vec(iCol, vtl);
			matrix::base::mat_mul(sol, vtl, vl);
			sol_orth.set_col(iCol, vl);
		}
		//wxLogMessage(_T("Trans mat new:%s\n===========\n"), trans_matrix.to_wxstr());
	}

};

// reconstruction is t-consuming - to obtain reconstructed solution
// reconstruct should be invoked - this is only needed when 
// eigenfunctions are of interest (in addition to eigenvalues)

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
	std::vector<t_MatCmplx> ret(_nnodes, t_MatCmplx(_NSolVecs, 2*_NSolVecs));

	//int orthStackInd = _orthStack.size()-1;
	int orthStackInd = _orth_stack_len - 1;

	t_SqMatCmplx direct(_NSolVecs);
	direct.setToUnity();

	t_SqMatCmplx tmp(_NSolVecs);

	for (int i=_nnodes-1; i>=0; i--){
		ret[i] = solution[i]*direct;
		if (orthStackInd>=0){
			if (i==_orthStack[orthStackInd].ind){
				const t_SqMatCmplx& l = _orthStack[orthStackInd].orthMatrix;
				// bug fixed in smallmat lib
				//direct=l*direct;
				// better make the same without implicit matrix operations 
				matrix::base::mat_mul<t_Complex, t_Complex>(_orthStack[orthStackInd].orthMatrix,direct,tmp);
				direct.set(tmp);
				orthStackInd--;
			}
		}
	}
	return ret;
}

// Gramm - Shmidt determinant of rank k, basis v0, v1, ...
// Ermith matrix
//	( (v0, v0)		(v1, v0)*		 ... (v[k-1], v0)*		)
//	(  ...													)
//	(  (v[k-1], v0) (v[k-1], v1)	 ... (v[k-1], v[k-1])*	)

t_Complex t_ODES::detGS(const t_MatCmplx& sol, const int& rank) const {

#ifdef _DEBUG
	if (rank>sol.nCols()) wxLogMessage(_T("GS determinant error: Rank exceeds task dimension\n"));
	if (rank<0) wxLogMessage(_T("GS determinant error: Rank < 0 \n"));
#endif

	if (rank == 0) return 1.0;
	t_SqMatCmplx mat(rank);
	int sol_nrows = sol.nRows();
	t_VecCmplx sol_i(sol_nrows), sol_j(sol_nrows);
	for (int i = 0; i<rank; i++) {
		for (int j = 0; j <= i; j++) {
			sol.col_to_vec(i, sol_i);
			sol.col_to_vec(j, sol_j);
			//t_Complex val = vector::dot(t_VecCmplx(sol[i]), t_VecCmplx(sol[j]));
			t_Complex val = vector::dot(sol_i, sol_j);
			// Ermith Matrix
			mat[i][j] = val;
			mat[j][i] = std::conj(val);
		}
	}

	//t_Complex det = mat.det();
	// IMPORTANT TODO: EXPLORE WHY!
	//if (abs(det.imag()/det.real())>1.0e-6){
	//std::wostringstream ostr;
	//ostr<<_T("GS determinant error: determinant is complex; imag:")
	//<<abs(det.imag())<<_T("\n");
	//log_my::CoutMessage(ostr.str());
	//}
	return mat.det();
}

t_Complex t_ODES::minorGS(const t_MatCmplx& vectors, const int& dim, const int& nExcludeCol) const {
	if (dim == 0) return 1.0;
	t_SqMatCmplx ret(dim);
	int vectors_nrows = vectors.nRows();
	t_VecCmplx vec_i(vectors_nrows), vec_j(vectors_nrows);
	for (int i = 0; i<nExcludeCol; i++) {
		for (int j = 0; j<dim; j++) {
			vectors.col_to_vec(i, vec_i);
			vectors.col_to_vec(j, vec_j);
			//ret[i][j] = vector::dot(t_VecCmplx(vectors[i]), t_VecCmplx(vectors[j]));
			ret[i][j] = vector::dot(vec_i, vec_j);
		}
	}
	for (int i = nExcludeCol + 1; i <= dim; i++) {
		for (int j = 0; j<dim; j++) {
			vectors.col_to_vec(i, vec_i);
			vectors.col_to_vec(j, vec_j);
			//ret[i-1][j] = vector::dot(t_VecCmplx(vectors[i]), t_VecCmplx(vectors[j]));
			ret[i - 1][j] = vector::dot(vec_i, vec_j);
		}
	}
	return ret.det();
}