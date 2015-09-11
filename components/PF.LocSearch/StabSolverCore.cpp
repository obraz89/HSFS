#include "stdafx.h"

#include "StabSolver.h"
#include "log.h"

using namespace pf;


void t_StabSolver::_setStabMatrix3D(const t_ProfRec& rec){
	t_CompVal imagUnity(0.0, 1.0);

	const double& stabRe = _profStab.scales().ReStab;
	const double& Me = _profStab.scales().Me;

	const mf::t_FldParams& Params = _rFldNS.get_mf_params();

	const double gMaMa = Params.Gamma*Me*Me;
	const double g_1MaMa = (Params.Gamma-1.0)*Me*Me;
	const double Pr = Params.Pr;

	const double inv_t = 1.0/rec.t;
	const double inv_mu = 1.0/rec.mu;

	const double vCoefL = 2.0/3.0*(0.0+2.0); //TODO: second Visc coef instead of 0.0
	const double vCoefS = 2.0/3.0*(0.0-1.0); // 
	const double vCoefM = 1.0+vCoefS;

	const t_CompVal& alpha = _waveChars.a;
	const t_CompVal& beta = _waveChars.b;
	const t_CompVal& freq = _waveChars.w;

	const t_CompVal dEicon_dt = alpha*rec.u + beta*rec.w - freq;

	const double dMu1 = rec.mu1*rec.t1;
	const double dMu2 = rec.mu2*rec.t1*rec.t1 + rec.mu1*rec.t2;
	const double dMu_U= rec.mu2*rec.t1*rec.u1 + rec.mu1*rec.u2;

	const t_CompVal xi = 1.0/(stabRe*inv_mu+imagUnity*vCoefL*gMaMa*dEicon_dt);

	// compose stab matrix for direct problem first

	// first
	_stab_matrix[1][0]=1.0;
	// second
	_stab_matrix[0][1] = imagUnity*stabRe*inv_t*inv_mu*dEicon_dt +
		pow(alpha,2) + pow(beta,2);

	_stab_matrix[1][1] =  - inv_mu*dMu1;

	_stab_matrix[2][1] = stabRe*inv_t*inv_mu*rec.u1 - 
		imagUnity*alpha*(inv_mu*dMu1+vCoefM*inv_t*rec.t1);

	_stab_matrix[3][1] = alpha*(imagUnity*stabRe*inv_mu -
		vCoefM*gMaMa*dEicon_dt);

	_stab_matrix[4][1] = vCoefM*alpha*inv_t*dEicon_dt - inv_mu*dMu_U;

	_stab_matrix[5][1] = -inv_mu*rec.mu1*rec.u1;
	// third
	_stab_matrix[0][2] = -imagUnity*alpha;

	_stab_matrix[2][2] = inv_t*rec.t1;

	_stab_matrix[3][2] = -imagUnity*gMaMa*dEicon_dt;

	_stab_matrix[4][2] = imagUnity*inv_t*dEicon_dt;

	_stab_matrix[6][2] = -imagUnity*beta;
	// fourth
	_stab_matrix[0][3] = -imagUnity*xi*alpha*
		(2.0*inv_mu*dMu1 + vCoefL*inv_t*rec.t1);

	_stab_matrix[1][3] = -imagUnity*xi*alpha;

	_stab_matrix[2][3] = xi*(-alpha*alpha-beta*beta+
		vCoefL*inv_t*inv_mu*dMu1*rec.t1+
		vCoefL*inv_t*rec.t2-
		imagUnity*stabRe*inv_t*inv_mu*dEicon_dt);

	_stab_matrix[3][3] = -imagUnity*xi*vCoefL*gMaMa*
		(
		(inv_mu*dMu1+inv_t*rec.t1)*dEicon_dt+
		alpha*rec.u1+beta*rec.w1
		);

	_stab_matrix[4][3] = imagUnity*xi*
		(
		(inv_mu*rec.mu1+vCoefL*inv_t)*(alpha*rec.u1+beta*rec.w1)+
		vCoefL*inv_t*inv_mu*dMu1*dEicon_dt
		);

	_stab_matrix[5][3] = imagUnity*xi*vCoefL*inv_t*dEicon_dt;

	_stab_matrix[6][3] = -imagUnity*xi*beta*
		(2.0*inv_mu*dMu1 + vCoefL*inv_t*rec.t1);

	_stab_matrix[7][3] = -imagUnity*xi*beta;
	// fifth
	_stab_matrix[5][4]=1.0;
	// sixth
	_stab_matrix[1][5] = -2.0*Pr*g_1MaMa*rec.u1;

	_stab_matrix[2][5] = Pr*
		(
		stabRe*inv_t*inv_mu*rec.t1 - 
		2.0*imagUnity*g_1MaMa*(alpha*rec.u1+beta*rec.w1)
		);

	_stab_matrix[3][5] = -imagUnity*stabRe*Pr*inv_mu*g_1MaMa*dEicon_dt;

	_stab_matrix[4][5] = imagUnity*stabRe*Pr*inv_t*inv_mu*dEicon_dt +
		alpha*alpha + beta*beta - 
		g_1MaMa*Pr*inv_mu*rec.mu1*(rec.u1*rec.u1+rec.w1*rec.w1)-
		inv_mu*dMu2;

	_stab_matrix[5][5] = -2.0*inv_mu*dMu1;

	_stab_matrix[7][5] = -2.0*Pr*g_1MaMa*rec.w1;
	// seventh
	_stab_matrix[7][6]=1.0;
	// last
	_stab_matrix[2][7] = -imagUnity*beta*(inv_mu*dMu1+vCoefM*inv_t*rec.t1)+
		stabRe*inv_mu*inv_t*rec.w1;

	_stab_matrix[3][7] = imagUnity*stabRe*beta*inv_mu-
		vCoefM*beta*gMaMa*dEicon_dt;

	_stab_matrix[4][7] = vCoefM*beta*inv_t*dEicon_dt-
		inv_mu*(rec.mu2*rec.t1*rec.w1+rec.mu1*rec.w2);

	_stab_matrix[5][7] = -inv_mu*rec.mu1*rec.w1;

	_stab_matrix[6][7] = imagUnity*stabRe*inv_t*inv_mu*dEicon_dt +
		alpha*alpha + beta*beta;

	_stab_matrix[7][7] = -inv_mu*dMu1;

	// option for conjugate problem
	// dpsi/dy = H_conj*psi, 
	// H_conj = -1*(H_dir)_transposed_conjugate, H_dir - matrix of direct problem
	// we will solve dpsi_conjugate/dy = A*psi_conjugate
	// then A = -1*H_dir_transposed
	// the result is psi_conjugate : complex conjugate vetor of explicit conjugate vector psi 

	if (_ls_mode.is_flag_on(stab::t_LSMode::CONJUGATE)){

		_stab_matrix.transpose();
		_stab_matrix.mul_by_factor(-1.0);

	}

}

void t_StabSolver::_setStabMatrix3D(const double& a_y){

	const t_ProfRec& rec = _profStab.get_rec(a_y);
	_setStabMatrix3D(rec);
};


// TODO: to speed up - rewrite to "_setAsymptotcs" 
// and write directly to destination formal param

t_MatCmplx t_StabSolver::_getAsymptotics3D
(const t_WCharsLoc& a_waveChars, t_ASYM_MODE mode){
	const int dim = getTaskDim();
	t_MatCmplx initial_vectors(dim,2*dim);
	t_SqMatCmplx b_coef(dim);
	t_VecCmplx lambda(dim,0.0);

	// TODO: function for simplified asymp: u=1.0, u'=0, u''=0, ... ?
	t_ProfRec out_rec = _profStab.get_last_rec();

	if (_ls_mode.is_flag_on(stab::t_LSMode::ASYM_HOMOGEN)){
		// modify ... a little
		// this is how AVF did it
		out_rec.t1=0.0;
		out_rec.t2=0.0;
		out_rec.u1=0.0;
		out_rec.u2=0.0;
		out_rec.w1=0.0;
		out_rec.w2=0.0;
	}
	_setStabMatrix3D(out_rec);
	// to shorten 
	t_SqMatCmplx& _sm = _stab_matrix;

	// det(H - lambda*E)=0 <=> lambda1, lambda2...
	// eigen vectors h1, h2, h3, h4 are used as initials

	b_coef[0][0]=_stab_matrix[0][1];
	b_coef[1][0]=_stab_matrix[3][1];
	b_coef[2][0]=_stab_matrix[4][1];

	b_coef[1][1]=_sm[3][1]*_sm[1][3]+_sm[3][2]*_sm[2][3]+
		_sm[3][5]*_sm[5][3]+_sm[3][7]*_sm[7][3];

	b_coef[2][1]=_sm[4][1]*_sm[1][3]+_sm[4][2]*_sm[2][3]+
		_sm[5][3]*_sm[4][5]+_sm[7][3]*_sm[4][7];

	b_coef[1][2]=_sm[3][5];
	b_coef[2][2]=_sm[4][5];

	b_coef[1][3]=_sm[3][7];
	b_coef[2][3]=_sm[4][7];

	b_coef[3][3] = _stab_matrix[0][1];
	//

	t_CompVal s1 = 0.5*(b_coef[1][1]+b_coef[2][2]);
	t_CompVal s2 = sqrt(
		0.25*std::pow(b_coef[1][1]-b_coef[2][2],2)+
		b_coef[2][1]*b_coef[1][2]
	);
	lambda[0] = -sqrt(b_coef[0][0]);
	lambda[1] = -sqrt(s1+s2);
	lambda[2] = -sqrt(s1-s2);
	lambda[3] = lambda[0];

	b_coef[0][0] = 1.0;
	b_coef[1][0] = 0.0;
	b_coef[2][0] = 0.0;
	b_coef[3][0] = 0.0;
	for (int i=1; i<3; i++){
		t_CompVal L2 = std::pow(lambda[i],2);
		t_CompVal denom = _stab_matrix[0][1] - L2;
		b_coef[0][i] = (
			(L2 - _stab_matrix[4][5])*_stab_matrix[3][1]+
			_stab_matrix[4][1]*_stab_matrix[3][5]
		)/denom;

		b_coef[1][i] = _stab_matrix[4][5] - L2;
		b_coef[2][i] = -_stab_matrix[3][5];
		b_coef[3][i] = (
			_stab_matrix[3][5]*_stab_matrix[4][7]+
			(L2 - _stab_matrix[4][5])*_stab_matrix[3][7]
		)/denom;
	};

	b_coef[0][3]=0.0;
	b_coef[1][3]=0.0;
	b_coef[2][3]=0.0;
	b_coef[3][3]=1.0;

	for (int i=0; i<4; i++){	
		initial_vectors[i][0] = b_coef[0][i];
		initial_vectors[i][1] = lambda[i]*b_coef[0][i];
		initial_vectors[i][2] = (
			_stab_matrix[0][2]*b_coef[0][i]+
			_stab_matrix[3][2]*b_coef[1][i]+
			_stab_matrix[4][2]*b_coef[2][i]+
			_stab_matrix[6][2]*b_coef[3][i]
		)/lambda[i];

		initial_vectors[i][3] = b_coef[1][i];
		initial_vectors[i][4] = b_coef[2][i];
		initial_vectors[i][5] = lambda[i]*b_coef[2][i];
		initial_vectors[i][6] = b_coef[3][i];
		initial_vectors[i][7] = (
			_stab_matrix[3][7]*b_coef[1][i]+
			_stab_matrix[4][7]*b_coef[2][i]+
			_stab_matrix[6][7]*b_coef[3][i]
		)/lambda[i];
	}

	// verification of asymptotics, enable when debug needed
	/*
	t_VecCmplx asym_resid_v(2*dim, 0.0), init_vec(2*dim, 0.0);
	t_VecCmplx v1(2*dim, 0.0), v2(2*dim, 0.0), v3(2*dim, 0.0);
	double resid=0.0;

	for (int i=0; i<dim; i++){

			initial_vectors.col_to_vec(i, init_vec);
			//asym_resid_v = _stab_matrix*init_vec - lambda[i]*init_vec;		
			matrix::base::mat_mul<t_Complex, t_Complex>(_stab_matrix, init_vec, v1);
			matrix::base::mul<t_Complex, t_Complex>(lambda[i], init_vec, v2);
			matrix::base::minus<t_Complex, t_Complex>(v1, v2, asym_resid_v);

			resid = resid + asym_resid_v.norm().real();
	}

	wxLogMessage(_T("Verify asymptotics: resid = %f"), resid);
	//if (resid>ASYM_TOL_DEFAULT)
	//	ssuGENTHROW(_T("StabSolver Error: Verification of Asymptotics failed"));
	*/
	return initial_vectors;
}

t_Complex t_StabSolver::_calcResidual(t_Complex* out_resid_coefs) const{

	const std::vector<t_MatCmplx>& solution = _math_solver.solution;

	const t_MatCmplx& wall_func = solution[_math_solver.getNNodes()-1];

	t_VecCmplx rhs(4);
	t_VecCmplx resid_coefs(4);
	t_SqMatCmplx mat(4);

	t_Complex resid(0.0);

	if (!_ls_mode.is_flag_on(stab::t_LSMode::CONJUGATE)){

		// direct problem
		// get resid by temperature residual
		// construct matrix of 4 cols (u,u',v,w)
		// of solutions at wall
		// solve with rhs (0,1,0,0)
		// and construct resid of temperature

		rhs[1]=1.0;
		for (int i=0; i<4; i++){
			mat[i][0] = wall_func[i][0];
			mat[i][1] = wall_func[i][1];
			mat[i][2] = wall_func[i][2];
			mat[i][3] = wall_func[i][6];
		}
		matrix::base::mat_mul<t_CompVal, t_CompVal>(mat.inverse(),rhs,resid_coefs);

		for (int i=0; i<4; i++){
			resid+=resid_coefs[i]*wall_func[i][4];
		}

	}else{

		// conjugate problem
		// get resid by pressure residual
		// construct matrix of 4 cols (v,u',r,w')
		// of solutions at wall
		// solve with rhs (1,0,0,0)
		// and construct resid of pressure

		rhs[0]=1.0;
		for (int i=0; i<4; i++){
			mat[i][0] = wall_func[i][2];
			mat[i][1] = wall_func[i][1];
			mat[i][2] = wall_func[i][5];
			mat[i][3] = wall_func[i][7];
		}
		matrix::base::mat_mul<t_CompVal, t_CompVal>(mat.inverse(),rhs,resid_coefs);

		for (int i=0; i<4; i++){
			resid+=resid_coefs[i]*wall_func[i][3];
		}

	}

	// when needed write out residual coefs
	if (out_resid_coefs!=NULL)
		for (int i=0; i<4; i++)
			out_resid_coefs[i] = resid_coefs[i];

	return resid;

}

t_Complex t_StabSolver::solve(t_WCharsLoc& stab_point){
	this->_waveChars = stab_point;
	this->_math_solver.clean();
	this->_math_solver.setInitials(_getAsymptotics3D(stab_point));
	_math_solver.solve();
	t_Complex resid = _calcResidual();
	_waveChars.resid = resid;
	stab_point.resid = resid;
	return resid;
}

void t_StabSolver::dumpEigenFuctions(const std::string& fname){
	std::wofstream fstr(&fname[0]);
	fstr<<_T("u_re\tu_im\tu'_re\tu'_im\tv_re\tv_im\tp_re\tp_im\tt_re\tt_im\tr_re\tr_im\tw_re\tw_im\tw'_re\tw'_im\tY\n");

	std::vector<t_MatCmplx> solutions = _math_solver.reconstruct();
	int nvecs = getTaskDim();
	// TODO: ask AVF how to normalize
/*
	for (int i=0; i<nvecs;i++){
		for (int j=0; j<getNNodes(); j++){
			for (int k=0; k<2*nvecs; k++){
				if (std::norm(solutions[j][i][k])>std::norm(max_val)){
					max_val = solutions[j][i][k];
				};

			}
		}
	}
*/

	/*
	for (int i=0; i<nvecs;i++){
		for (int j=0; j<getNNodes(); j++){
			t_Complex val;
			for (int k=0; k<2*nvecs; k++){
				val = solutions[j][i][k];
				fstr<<std_manip::std_format_sci<double>(val.real())
					<<_T("\t")
					<<std_manip::std_format_sci<double>(val.imag())
					<<_T("\t");
			}
			fstr<<_math_solver.varRange[j];
			fstr<<_T("\n");
		}
		// separate solutions to simplify origin export
		fstr<<_T("\n\n\n\n\n\n");
	}*/

	// finally write down eigen solution satisfying bc conditions
	t_VecCmplx wall_coefs(4);
	_calcResidual(&wall_coefs[0]);

	t_VecCmplx cur_sol(2*nvecs);// ok

	for (int j=0; j<getNNodes(); j++){
		const t_VecCmplx* pVecs[4];

		//t_VecCmplx cur_sol = solutions[j]*wall_coefs;;
		matrix::base::mat_mul(solutions[j], wall_coefs, cur_sol);

		for (int k=0; k<2*nvecs; k++){
			fstr<<std_manip::std_format_sci<double>(cur_sol[k].real())
				<<_T("\t")
				<<std_manip::std_format_sci<double>(cur_sol[k].imag())
				<<_T("\t");
		}
		fstr<<_math_solver.varRange[j];
		fstr<<_T("\n");
	}

}

void t_StabSolver::dumpProfileStab(const std::string& fname) const{
	_profStab.dump(fname);
}


