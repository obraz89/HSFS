#include "stdafx.h"

#include "StabSolver.h"
#include "common_data.h"
#include "log.h"

static const int STAB_MATRIX_DIM = 8;

using namespace hsstab;
using namespace hsstab::cmpnts;

using namespace pf;

t_StabSolver::t_StabSolver(const mf::t_Block& a_rFldNS):
_rFldNS(a_rFldNS), _profStab(), 
_math_solver(), _stab_matrix(STAB_MATRIX_DIM), _params(){

	_math_solver._pStab_solver = this;

};

void t_StabSolver::init(const hsstab::TPlugin& g_plug){

	const hsstab::TPluginParamsGroup& g = g_plug.get_settings_grp_const("");

	ls::_init_ortho_ls_base_params(_params, g);

	//_init();	// init solver

}

const t_StabScales& t_StabSolver::get_stab_scales() const{
	return _profStab.scales();
};

// private, to be used in 3D context

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
	const double dMu2 = rec.mu2*rec.t1*rec.u1 + rec.mu1*rec.u2;
	const t_CompVal xi = 1.0/(stabRe*inv_mu+imagUnity*vCoefL+gMaMa*dEicon_dt);

	// first
	_stab_matrix[1][0]=1.0;
	// second
	_stab_matrix[0][1] = imagUnity*stabRe*inv_t*inv_mu*dEicon_dt +
		pow(alpha,2) + pow(beta,2);

	_stab_matrix[1][1] =  - inv_mu*rec.mu1*rec.t1;

	_stab_matrix[2][1] = stabRe*inv_t*inv_mu*rec.u1 - 
		imagUnity*alpha*(inv_mu*rec.mu1*rec.t1+vCoefM*inv_t*rec.t1);

	_stab_matrix[3][1] = alpha*(imagUnity*stabRe*inv_mu -
		vCoefM*gMaMa*dEicon_dt);

	_stab_matrix[4][1] = vCoefM*alpha*inv_t*dEicon_dt - inv_mu*dMu2;

	_stab_matrix[5][1] = -inv_mu*rec.mu1*rec.u1;
	// third
	_stab_matrix[0][2] = -imagUnity*alpha;

	_stab_matrix[2][2] = inv_t*rec.t1;

	_stab_matrix[3][2] = -imagUnity*gMaMa*dEicon_dt;

	_stab_matrix[4][2] = imagUnity*inv_t*dEicon_dt;

	_stab_matrix[6][2] = -imagUnity*beta;
	// fourth
	_stab_matrix[0][3] = -imagUnity*xi*alpha*rec.t1*
		(2.0*inv_mu*rec.mu1 + vCoefL*inv_t);

	_stab_matrix[1][3] = -imagUnity*xi*alpha;

	_stab_matrix[2][3] = xi*(-alpha*alpha-beta*beta+
		vCoefL*inv_t*inv_mu*rec.mu1*rec.t1*rec.t1+
		vCoefL*inv_t*rec.t2-
		imagUnity*stabRe*inv_t*inv_mu*dEicon_dt);

	_stab_matrix[3][3] = -imagUnity*xi*vCoefL*gMaMa*
		(
		(inv_mu*rec.mu1*rec.t1+inv_t*rec.t1)*dEicon_dt+
		alpha*rec.u1+beta*rec.w1
		);

	_stab_matrix[4][3] = imagUnity*xi*
		(
		(inv_mu*rec.mu1+vCoefL*inv_t)*(alpha*rec.u1+beta*rec.w1)+
		vCoefL*inv_t*inv_mu*rec.mu1*rec.t1*dEicon_dt
		);

	_stab_matrix[5][3] = imagUnity*xi*vCoefL*inv_t*dEicon_dt;

	_stab_matrix[6][3] = -imagUnity*xi*beta*
		(2.0*inv_mu*rec.mu1*rec.t1 + vCoefL*inv_t*rec.t1);

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
		inv_mu*(rec.mu2*rec.t1*rec.t1+rec.mu1*rec.t2);
	_stab_matrix[5][5] = -2.0*inv_mu*rec.mu1*rec.t1;

	_stab_matrix[7][5] = -2.0*Pr*g_1MaMa*rec.w1;
	// seventh
	_stab_matrix[7][6]=1.0;
	// last
	_stab_matrix[2][7] = -imagUnity*beta*(inv_mu*rec.mu1*rec.t1+vCoefM*inv_t*rec.t1)+
		stabRe*inv_mu*inv_t*rec.w1;

	_stab_matrix[3][7] = imagUnity*stabRe*beta*inv_mu-
		vCoefM*beta*gMaMa*dEicon_dt;

	_stab_matrix[4][7] = vCoefM*beta*inv_t*dEicon_dt-
		inv_mu*(rec.mu2*rec.t1*rec.w1+rec.mu1*rec.w2);

	_stab_matrix[5][7] = -inv_mu*rec.mu1*rec.w1;

	_stab_matrix[6][7] = imagUnity*stabRe*inv_t*inv_mu*dEicon_dt +
		alpha*alpha + beta*beta;

	_stab_matrix[7][7] = -inv_mu*rec.mu1*rec.t1;
}

void t_StabSolver::_setStabMatrix3D(const double& a_y){

	t_ProfRec& rec = _profStab.get_rec(a_y);
	_setStabMatrix3D(rec);
};

/*const t_SqMatrix& t_StabSolver::_getStabMatrix3D() const{
	return _stab_matrix;
}*/

// old - can be used only in local spectral procedures
/*void t_StabSolver::setMatSplitByW3D_88(const double& a_y, t_SqMatrix& mat_no_w, t_SqMatrix& mat_w) const{
	t_CompVal imagUnity(0.0, 1.0);

	const double& stabRe = _profStab.stabRe;
	const double& Me = _profStab.Me;
	const double gMaMa = MF_Field::Gamma*Me*Me;
	const double g_1MaMa = (MF_Field::Gamma-1.0)*Me*Me;

	const double u = _profStab.getValue(a_y, _profStab.u);
	const double u1 = _profStab.getValue(a_y, _profStab.u1);
	const double u2 = _profStab.getValue(a_y, _profStab.u2);

	const double w = _profStab.getValue(a_y, _profStab.w);
	const double w1 = _profStab.getValue(a_y, _profStab.w1);
	const double w2 = _profStab.getValue(a_y, _profStab.w2);

	const double t = _profStab.getValue(a_y, _profStab.t);
	const double t1 = _profStab.getValue(a_y, _profStab.t1);
	const double t2 = _profStab.getValue(a_y, _profStab.t2);

	const double mu = _profStab.getValue(a_y, _profStab.mu);
	const double mu1 = _profStab.getValue(a_y, _profStab.mu1);
	const double mu2 = _profStab.getValue(a_y, _profStab.mu2);

	const double inv_t = 1.0/t;
	const double inv_mu = 1.0/mu;

	const double vCoefL = 2.0/3.0*(0.0+2.0); //TODO: second Visc coef instead of 0.0
	const double vCoefS = 2.0/3.0*(0.0-1.0); // 
	const double vCoefM = 1.0+vCoefS;

	const t_CompVal& alpha = _waveChars.a;
	const t_CompVal& beta = _waveChars.b;

	const t_CompVal dEicon_no_W = alpha*u + beta*w ;
	const t_CompVal dEicon_W = -1.0;
	const double dMu2 = mu2*t1*u1 + mu1*u2;
	const t_CompVal xi = 1.0/(stabRe*inv_mu+imagUnity*vCoefL+gMaMa*dEicon_no_W);
// set no_W matrix
// of size 5rowsx8cols
	// first <--> "2"
	mat_no_w[0][0] = imagUnity*stabRe*inv_t*inv_mu*dEicon_no_W +
		pow(alpha,2) + pow(beta,2);

	mat_no_w[1][0] =  - inv_mu*mu1*t1;

	mat_no_w[2][0] = stabRe*inv_t*inv_mu*u1 - 
		imagUnity*alpha*(inv_mu*mu1*t1+vCoefM*inv_t*t1);

	mat_no_w[3][0] = alpha*(imagUnity*stabRe*inv_mu -
		vCoefM*gMaMa*dEicon_no_W);

	mat_no_w[4][0] = vCoefM*alpha*inv_t*dEicon_no_W - inv_mu*dMu2;

	mat_no_w[5][0] = -inv_mu*mu1*u1;
	// second <--> "3"
	mat_no_w[0][1] = -imagUnity*alpha;

	mat_no_w[2][1] = inv_t*t1;

	mat_no_w[3][1] = -imagUnity*gMaMa*dEicon_no_W;

	mat_no_w[4][1] = imagUnity*inv_t*dEicon_no_W;

	mat_no_w[6][1] = -imagUnity*beta;
	// third <--> "4"
	mat_no_w[0][2] = -imagUnity*xi*alpha*t1*
		(2.0*inv_mu*mu1 + vCoefL*inv_t);

	mat_no_w[1][2] = -imagUnity*xi*alpha;

	mat_no_w[2][2] = xi*(-alpha*alpha-beta*beta+
		vCoefL*inv_t*inv_mu*mu1*t1*t1+
		vCoefL*inv_t*t2-
		imagUnity*stabRe*inv_t*inv_mu*dEicon_no_W);

	mat_no_w[3][2] = -imagUnity*xi*vCoefL*gMaMa*
		(
		(inv_mu*mu1*t1+inv_t*t1)*dEicon_no_W+
		alpha*u1+beta*w1
		);

	mat_no_w[4][2] = imagUnity*xi*
		(
		(inv_mu*mu1+vCoefL*inv_t)*(alpha*u1+beta*w1)+
		vCoefL*inv_t*inv_mu*mu1*t1*dEicon_no_W
		);

	mat_no_w[5][2] = imagUnity*xi*vCoefL*inv_t*dEicon_no_W;

	mat_no_w[6][2] = -imagUnity*xi*beta*
		(2.0*inv_mu*mu1*t1 + vCoefL*inv_t*t1);

	mat_no_w[7][2] = -imagUnity*xi*beta;
	// fourth <--> "6"
	mat_no_w[1][3] = -2.0*MF_Field::Pr*g_1MaMa*u1;

	mat_no_w[2][3] = MF_Field::Pr*
		(
		stabRe*inv_t*inv_mu*t1 - 
		2.0*imagUnity*g_1MaMa*(alpha*u1+beta*w1)
		);

	mat_no_w[3][3] = -imagUnity*stabRe*MF_Field::Pr*inv_mu*g_1MaMa*dEicon_no_W;

	mat_no_w[4][3] = imagUnity*stabRe*MF_Field::Pr*inv_t*inv_mu*dEicon_no_W +
		alpha*alpha + beta*beta - 
		g_1MaMa*MF_Field::Pr*inv_mu*mu1*(u1*u1+w1*w1)-
		inv_mu*(mu2*t1*t1+mu1*t2);
	mat_no_w[5][3] = -2.0*inv_mu*mu1*t1;

	mat_no_w[7][3] = -2.0*MF_Field::Pr*g_1MaMa*w1;
	// fifth <--> "8"
	mat_no_w[2][4] = -imagUnity*beta*(inv_mu*mu1*t1+vCoefM*inv_t*t1)+
		stabRe*inv_mu*inv_t*w1;

	mat_no_w[3][4] = imagUnity*stabRe*beta*inv_mu-
		vCoefM*beta*gMaMa*dEicon_no_W;

	mat_no_w[4][4] = vCoefM*beta*inv_t*dEicon_no_W-
		inv_mu*(mu2*t1*w1+mu1*w2);

	mat_no_w[5][4] = -inv_mu*mu1*w1;

	mat_no_w[6][4] = imagUnity*stabRe*inv_t*inv_mu*dEicon_no_W +
		alpha*alpha + beta*beta;

	mat_no_w[7][4] = -inv_mu*mu1*t1;
// set W matrix
// of the same sizes
// 5rows x 8 cols
	// first <--> part of "2", 3 non-zeros
	mat_w[0][0] = imagUnity*stabRe*inv_t*inv_mu*dEicon_W;

	mat_w[3][0] = alpha*( -	vCoefM*gMaMa*dEicon_W);

	mat_w[4][0] = vCoefM*alpha*inv_t*dEicon_W;

	// second <--> part of "3", 2 non-zeros

	mat_w[3][1] = -imagUnity*gMaMa*dEicon_W;

	mat_w[4][1] = imagUnity*inv_t*dEicon_W;

	// third <--> "4", 4 non-zeros, a lot)

	mat_w[2][2] = xi*(-imagUnity*stabRe*inv_t*inv_mu*dEicon_W);

	mat_w[3][2] = -imagUnity*xi*vCoefL*gMaMa*
		(
		(inv_mu*mu1*t1+inv_t*t1)*dEicon_W
		);

	mat_w[4][2] = imagUnity*xi*
		(
		vCoefL*inv_t*inv_mu*mu1*t1*dEicon_W
		);

	mat_w[5][2] = imagUnity*xi*vCoefL*inv_t*dEicon_W;

	// fourth <--> part of "6", 2 non-zeros

	mat_w[3][3] = -imagUnity*stabRe*MF_Field::Pr*inv_mu*g_1MaMa*dEicon_W;

	mat_w[4][3] = imagUnity*stabRe*MF_Field::Pr*inv_t*inv_mu*dEicon_W;

	// fifth <--> "8", 3 non-zeros

	mat_w[3][4] = -vCoefM*beta*gMaMa*dEicon_W;

	mat_w[4][4] = vCoefM*beta*inv_t*dEicon_W;

	mat_w[6][4] = imagUnity*stabRe*inv_t*inv_mu*dEicon_W;

	return;
};*/

t_VecCmplx t_StabSolver::_formRHS3D(const double& a_y, const t_VecCmplx& a_vars){
	_setStabMatrix3D(a_y);
	return _stab_matrix*a_vars; 
};

t_MatCmplx t_StabSolver::_getAsymptotics3D
(const t_WCharsLoc& a_waveChars, t_ASYM_MODE mode){
	const int dim = getTaskDim();
	t_MatCmplx initial_vectors(dim,2*dim);
	t_SqMatCmplx b_coef(dim);
	t_VecCmplx lambda(dim,0.0);
	const double& y_e = _profStab.get_thick();
	// TODO: function for simplified asymp: u=1.0, u'=0, u''=0, ... ?
	t_ProfRec out_rec = _profStab.get_rec(y_e);
	if (mode==ASYM_FORCE_SELF_SIM){
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
	_verifyAsymptotics3D(initial_vectors, lambda);
	return initial_vectors;
}

void t_StabSolver::_verifyAsymptotics3D
(const t_MatCmplx& init_vecs, const t_VecCmplx& lambda) const{
	const int dim = getTaskDim();
	t_MatCmplx asym_resid(dim, 2*dim, 0.0);
	t_SqMatCmplx lambda_mat(dim);
	for (int i=0; i<dim; i++){
		lambda_mat[i][i]=lambda[i];
	}
	asym_resid = _stab_matrix*init_vecs - init_vecs*lambda_mat;
	double resid=0.0;
	for (int i=0; i<asym_resid.nCols(); i++){
		for (int j=0; j<asym_resid.nRows(); j++){
			double cur_resid = smat::norm(asym_resid[i][j]);
			if (cur_resid>resid) resid = cur_resid;
		}
	}
	//if (resid>ASYM_TOL_DEFAULT)
	//	ssuGENTHROW(_T("StabSolver Error: Verification of Asymptotics failed"));
}


/************************************************************************/   
// Restore wave chars in global reference frame
// Use time approach
/************************************************************************/

t_WCharsGlob t_StabSolver::popGlobalWCharsTime(const t_ProfileNS& a_rProfNS){

	t_WCharsLoc restore_wave = _waveChars;
	t_WCharsLoc adjust_wave = _waveChars;

	if (_waveChars.get_treat()==stab::t_TaskTreat::SPAT){

		calcGroupVelocity(adjust_wave);

		adjust_wave.to_time();

		adjustLocal(adjust_wave, t_StabSolver::t_MODE::W_MODE);

	}

	t_WCharsGlob glob_wave(adjust_wave, a_rProfNS.get_mtr(), _profStab.scales());

	_waveChars = restore_wave;

	return glob_wave;

};


/************************************************************************/   
// Restore wave chars in global reference frame
// Use spat approach
/************************************************************************/

t_WCharsGlob t_StabSolver::popGlobalWCharsSpat(const t_ProfileNS& a_rProfNS){

	t_WCharsLoc restore_wave = _waveChars;
	t_WCharsLoc adjust_wave = _waveChars;

	if (_waveChars.get_treat()==stab::t_TaskTreat::TIME){

		calcGroupVelocity(adjust_wave);

		adjust_wave.to_spat();

		// TODO: here is another uncertainty:
		// A_MODE or B_MODE or mixed - how to choose
		adjustLocal(adjust_wave, t_StabSolver::t_MODE::A_MODE);

	}

	t_WCharsGlob glob_wave(adjust_wave, a_rProfNS.get_mtr(), _profStab.scales());

	_waveChars = restore_wave;

	return glob_wave;
};

/************************************************************************/   
// Set up the required context for a given point
// 
/************************************************************************/

void t_StabSolver::setContext
(const mf::t_BlkInd a_ind){

	int nnodes_stab = _params.NNodes;

	_math_solver.setContext(nnodes_stab);

	t_ProfileNS profNS(_rFldNS);
	profNS.initialize(a_ind, _params.ThickCoef);
	_profStab.initialize(profNS, nnodes_stab);

	if (_stab_matrix.nCols()!=STAB_MATRIX_DIM){
		_stab_matrix.resize(STAB_MATRIX_DIM);
	}

	for (int j=0; j<nnodes_stab; j++){
		_math_solver.varRange[j] = _profStab.get_y(nnodes_stab-1-j);
	}

	return;
}

/************************************************************************/   
// Set up context using AVF profiles
// 
/************************************************************************/

void t_StabSolver::setContext(const t_ProfileStab* a_prof_stab){

	_profStab = *a_prof_stab;

	int nnodes = _profStab.get_nnodes();

	_math_solver.setContext(nnodes);

	if (_stab_matrix.nCols()!=STAB_MATRIX_DIM){
		_stab_matrix.resize(STAB_MATRIX_DIM);
	}

	for (int j=0; j<nnodes; j++){
		_math_solver.varRange[j] = _profStab.get_y(nnodes-1-j);
	}

	return;
	
}

/************************************************************************/   
//COMPUTATION OF GROUP VELOCITY (VA,VB)=(DW/DA,DW/DB) 
//
/************************************************************************/

void t_StabSolver::calcGroupVelocity
(t_WCharsLoc &a_wave_chars){
	// empiric constant - tolerance
	const double dd = 0.001*a_wave_chars.a.real();//DELTA_SMALL;
	// ensure that we are in eigen 
	adjustLocal(a_wave_chars, W_MODE);

	t_WCharsLoc left_chars = a_wave_chars;
	t_WCharsLoc rght_chars = a_wave_chars;
	// dalpha
	left_chars.a-=dd;
	rght_chars.a+=dd;
    adjustLocal(left_chars, W_MODE);
	adjustLocal(rght_chars, W_MODE);
	if (left_chars.w == rght_chars.w){
		throw t_UnPhysWave();
	};
	a_wave_chars.vga = 0.5*(rght_chars.w - left_chars.w)/dd;

	// reset left, right
	left_chars = a_wave_chars;
	rght_chars = a_wave_chars;

	left_chars.b-=dd;
	rght_chars.b+=dd;
    adjustLocal(left_chars, W_MODE);
	adjustLocal(rght_chars, W_MODE);
	a_wave_chars.vgb = 0.5*(rght_chars.w - left_chars.w)/dd;
}

/************************************************************************/   
// Search maximum instability with keeping Re(w)=const
// Use time approach
/************************************************************************/

t_WCharsLoc t_StabSolver::_getStationaryMaxInstabTime
(const t_WCharsLoc& initial_guess){
	t_WCharsLoc nxt_wave = initial_guess;
	t_WCharsLoc prv_wave;
	t_Complex gv;
	double da, db; 

	int ndebug=0;

	do{
		prv_wave = nxt_wave;
		adjustLocal(prv_wave, W_MODE);   
		calcGroupVelocity(prv_wave);
		// IMPORTANT TODO: empirics 1% old val
		db = 0.01*prv_wave.a.real();
		double coef = prv_wave.vgb.real()/prv_wave.vga.real();
		da = -coef*db;

		nxt_wave.a = prv_wave.a + da;
		nxt_wave.b = prv_wave.b + db;
		nxt_wave.w = prv_wave.w;
		adjustLocal(nxt_wave, W_MODE);

		if (nxt_wave.w.imag()<prv_wave.w.imag()){
			db = -db;
			da = -coef*db;
			nxt_wave.a = prv_wave.a + da;
			nxt_wave.b = prv_wave.b + db;
			nxt_wave.w = prv_wave.w;
			adjustLocal(nxt_wave, W_MODE);	
		}	
		//dbg
		if (ndebug++==100){
			ndebug=0;
			std::wcout<<_T("getStatMaxInstabTime 100 iters done, cur wave:\n")<<nxt_wave<<_T("\n");
		}
	} while(prv_wave.w.imag()<nxt_wave.w.imag());

	double wr_error = abs(1.0-initial_guess.w.real()/prv_wave.w.real());
	std::wcout<<_T("getMaxWave Fixed Wr:err~")<<wr_error<<_T("\n")<<prv_wave;

	return prv_wave;
};

/************************************************************************/   
// Search maximum instability without any restrictions 
// Use time approach
/************************************************************************/

t_WCharsLoc t_StabSolver::_getMaxInstabTime_Stat
(const t_WCharsLoc &init_guess){
	// empiric half percent per iteration
	double dar = 0.005*init_guess.a.real();
	t_WCharsLoc base_wave = _getStationaryMaxInstabTime(init_guess);
	calcGroupVelocity(base_wave);
	t_WCharsLoc next_wave = base_wave;
	next_wave.a+=dar;
	next_wave.w+=base_wave.vga*dar;
	next_wave = _getStationaryMaxInstabTime(next_wave);
	// do we move in proper direction ?
	if (next_wave.w.imag()<base_wave.w.imag()){
		dar = -dar;
	};
	
	do{
		base_wave = next_wave;
		calcGroupVelocity(base_wave);
		next_wave.a+=dar;
		next_wave.w+=base_wave.vga*dar;
		next_wave = _getStationaryMaxInstabTime(next_wave);
	} while (base_wave.w.imag()>next_wave.w.imag());
	return base_wave;
};


/************************************************************************/   
// Search for exact eigen solution given initial approach a_wave_chars
// input: mode defines the value to be adjusted
// Use time approach
/************************************************************************/

bool t_StabSolver::adjustLocal
(t_WCharsLoc &a_wave_chars, t_StabSolver::t_MODE a_mode){
	// parameters
	const t_Complex d_arg = _params.AdjustStep;
	const int max_iters = _params.AdjustMaxIter;
	const double tol = _params.AdjustTol;
	// ensure that there is a base residual
	solve(a_wave_chars);
	t_WCharsLoc backup = a_wave_chars;
	t_Complex *choose_arg;
	if (a_mode==A_MODE){
		choose_arg = &(a_wave_chars.a);
	}else{
		if (a_mode==B_MODE){
			choose_arg = &(a_wave_chars.b);
		}else{
			choose_arg = &(a_wave_chars.w);
		};
	};
	t_Complex& arg = *choose_arg;
	t_Complex arg_base = arg;
	t_Complex res_base = a_wave_chars.resid;
	t_Complex deriv;
	for (int i=0; i<max_iters; i++){
		// check if we are converged
		if (smat::norm(a_wave_chars.resid)<tol){
			//std::wcout<<_T("Adjust Converged")<<res_base<<std::endl;
			return true;
		};
		// newton iteration
		arg+=d_arg;
		solve(a_wave_chars);
		deriv = d_arg/(a_wave_chars.resid - res_base);
		arg = arg_base - res_base*deriv;
		solve(a_wave_chars);
		arg_base = arg;
		res_base = a_wave_chars.resid;
	};
	wxLogMessage(_T("Local Search Error: no convergence\n"));
	a_wave_chars = backup;
	return false;
};

void t_StabSolver::setInitWaves(const std::vector<t_WCharsLoc>& a_inits){
	this->_initWaves = a_inits;
};

std::vector<t_WCharsLoc> t_StabSolver::filterInitWaves
(const std::vector<t_WCharsLoc>& all_initials){
	const double close_tol = 0.1;
	std::vector<t_WCharsLoc> good_initials;
	std::vector<t_WCharsLoc>::const_iterator iter;
	t_WCharsLoc init_wave;
	for (iter = all_initials.begin(); iter<all_initials.end(); iter++){
		const t_WCharsLoc& init_wave = *iter;
		t_WCharsLoc adj_wave = init_wave;
		try{
			adjustLocal(adj_wave, t_MODE::W_MODE);
			if (abs(adj_wave.w.real()-init_wave.w.real())<
				close_tol*abs(adj_wave.w.real()+init_wave.w.real())){
					good_initials.push_back(adj_wave);
			};
		}catch(t_UnPhysWave){
			// oh shit!
		};
	};
	return good_initials;
};


/************************************************************************/   
//Search for nearest eigen solution
// with keeping Re(w)=const
/************************************************************************/
void t_StabSolver::getEigenWFixed
(double wr_fixed, t_WCharsLoc& wave_chars, t_MODE mode){
	t_WCharsLoc backup = wave_chars;
	t_Complex *pVg, *pArg;
	if (mode==A_MODE){
		pArg = &wave_chars.a;
		pVg  = &wave_chars.vga;
	}else{
		if (mode==B_MODE){
			pArg = &wave_chars.b;
			pVg  = &wave_chars.vgb;
		}else{
			throw t_WrongMode();
		}
	}
	double resid=1.0;
	double dwr, darg;
	t_Complex vg;
	// IMPORTANT TODO: do we think Vg change is small??
	// or is it a mistake - Vg is not recalculated!!!
	adjustLocal(wave_chars, W_MODE);
	calcGroupVelocity(wave_chars);
	int n_iter=0;
	do{
		dwr = wr_fixed - wave_chars.w.real();
		vg = *pVg;
		darg = (dwr/vg).real();
		*pArg+= darg;
		wave_chars.w+= dwr;
		adjustLocal(wave_chars, W_MODE);
		resid = abs((dwr)/wr_fixed);
		if (n_iter++>_params.AdjustMaxIter){
			wxLogMessage(_T("In GetEigenWFixed: no convergence"));
			return;
		};
	} while (resid>=_params.AdjustTol);
};

/************************************************************************/   
// search for a wave chars with conditions provided
// with keeping Re(w)=const
// TODO: this is to be developed futher 
/************************************************************************/

// IMPORTANT TODO: task_mode switch, for now only time!
bool t_StabSolver::searchWave
(t_WCharsLoc& wchars, stab::t_LSCond cond, stab::t_TaskTreat task_mode)
{
	if (task_mode==stab::t_TaskTreat::SPAT){
		throw t_WrongMode();
	}

	switch (cond.get_mode())
	{
	case stab::t_LSCond::W_FIXED:
		// TODO: maybe in some cases B_MODE, how to switch ?
		getEigenWFixed(wchars.w.real(), wchars, t_StabSolver::A_MODE);
		return true;

	case (stab::t_LSCond::A_FIXED|stab::t_LSCond::B_FIXED):
		//return adjustLocal(wchars, W_MODE);
		return adjustLocal_Grad(wchars, cond);

	case (stab::t_LSCond::A_FIXED|stab::t_LSCond::W_FIXED):
		//return adjustLocal(wchars, B_MODE);
		return adjustLocal_Grad(wchars, cond);

	case (stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED):
		//return adjustLocal(wchars, A_MODE);
		return adjustLocal_Grad(wchars, cond);

	default:
		throw t_NotImplemented();
	}

};

void t_StabSolver::searchMaxWave(t_WCharsLoc& wchars, stab::t_LSCond cond, stab::t_TaskTreat task_mode){

	if (task_mode==stab::t_TaskTreat::SPAT){
		throw t_WrongMode();
	}
	switch (cond.get_mode())
	{

	case stab::t_LSCond::FREE:
		wchars = _getMaxInstabTime_Grad(wchars);
		return;

	case stab::t_LSCond::W_FIXED:
		wchars = _getStationaryMaxInstabTime(wchars);
		return;

	case stab::t_LSCond::A_FIXED:
		throw t_NotImplemented();

	case stab::t_LSCond::B_FIXED:
		throw t_NotImplemented();

	default:
		throw t_NotImplemented();
	}
}


t_WCharsLoc t_StabSolver::getMaxWave(const mf::t_BlkInd a_ind, 
									 const std::vector<t_WCharsLoc>& a_inits, const int& a_nnodesStab){
    // TODO:to implement somehow "default" parameter
	int nnodes_stab = (a_nnodesStab==0) ? _rFldNS.get_Ny() : a_nnodesStab;
	setContext(a_ind);
	setInitWaves(a_inits);
	t_WCharsLoc res_wave;
	//ensure
	res_wave.w = 0.0;
	t_WCharsLoc cur_wave;
	for (int i=0; i<_initWaves.size(); i++){
		cur_wave = _getMaxInstabTime_Grad(_initWaves[i]);
		if (cur_wave.w.imag()>res_wave.w.imag()){
			res_wave = cur_wave;
		}
	}
	return res_wave;
}
