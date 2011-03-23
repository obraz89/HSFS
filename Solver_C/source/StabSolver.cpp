#include "StabSolver.h"
//#include "SolverCore.h"
t_StabSolver::t_StabSolver(const MF_Field& a_rFldNS, t_StabField& a_rFldStab):
_rFldNS(a_rFldNS), _rFldStab(a_rFldStab), _profStab(0), _math_solver(){
	_math_solver._pStab_solver = this;
};
// obsolete
/*
void t_StabSolver::setParameters(t_ProfileStab& a_profStab){
	_profile = a_profStab;
	// globals -flow
	DNS.AMINF = MF_Field::Mach;
	DNS.REINF = MF_Field::Re;
	DNS.TINF = MF_Field::T_inf;
	DNS.XLL = MF_Field::L_ref;
	// -gas parameters
	OUT.G = 1.4;	// specific heat ratio
	OUT.SIG = 0.72;	// Prandtl number
	OUT.K = 0.0;	  // second viscosity
	OUT.VISC_FLAG = 1; // for viscosity: iupt=0 - Sutherland, 1 - power
	OUT.VISC_POWER = 0.75; // if VISC_FLAG=1, VISC_POWER - viscosity power coef  
	// locals
	NS.XST[0] = fld_ref.fld[i_ind][0][k_ind].x;
	CF_PROFILE_DATA.MF_WE = profile.w_e;
	CF_PROFILE_DATA.MF_ANG = fld_ref.get_cf_wave_dir(i_ind, k_ind);
	ADDITIONAL_NS.Y_DIM = profile.ny;
	MACKY1.NO = 150;	// number of grid points for stability computations
	HADY1.NS = 1;		// coef increasing number of points for stability 
	//computations (adaptation by NAVSTOK on uniform grid) 
	OUT.XI = 0.0;	// notused?
	OUT.PB = 0.0;	// ??
	ABOCT.MAB = 1;	// direct problem
};*/

// private, to be used in 3D context

t_SqMatrix t_StabSolver::getStabMatrix3D(const double& a_y) const{

	t_SqMatrix stab_matrix(8);
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
	const t_CompVal& freq = _waveChars.w;

	const t_CompVal dEicon_dt = alpha*u + beta*w - freq;
	const double dMu2 = mu2*t1*u1 + mu1*u2;
	const t_CompVal xi = 1.0/(stabRe*inv_mu+imagUnity*vCoefL+gMaMa*dEicon_dt);

	// first
	stab_matrix[1][0]=1.0;
	// second
	stab_matrix[0][1] = imagUnity*stabRe*inv_t*inv_mu*dEicon_dt +
		pow(alpha,2) + pow(beta,2);

	stab_matrix[1][1] =  - inv_mu*mu1*t1;

	stab_matrix[2][1] = stabRe*inv_t*inv_mu*u1 - 
		imagUnity*alpha*(inv_mu*mu1*t1+vCoefM*inv_t*t1);

	stab_matrix[3][1] = alpha*(imagUnity*stabRe*inv_mu -
		vCoefM*gMaMa*dEicon_dt);

	stab_matrix[4][1] = vCoefM*alpha*inv_t*dEicon_dt - inv_mu*dMu2;

	stab_matrix[5][1] = -inv_mu*mu1*u1;
	// third
	stab_matrix[0][2] = -imagUnity*alpha;

	stab_matrix[2][2] = inv_t*t1;

	stab_matrix[3][2] = -imagUnity*gMaMa*dEicon_dt;

	stab_matrix[4][2] = imagUnity*inv_t*dEicon_dt;

	stab_matrix[6][2] = -imagUnity*beta;
	// fourth
	stab_matrix[0][3] = -imagUnity*xi*alpha*t1*
		(2.0*inv_mu*mu1 + vCoefL*inv_t);

	stab_matrix[1][3] = -imagUnity*xi*alpha;

	stab_matrix[2][3] = xi*(-alpha*alpha-beta*beta+
		vCoefL*inv_t*inv_mu*mu1*t1*t1+
		vCoefL*inv_t*t2-
		imagUnity*stabRe*inv_t*inv_mu*dEicon_dt);

	stab_matrix[3][3] = -imagUnity*xi*vCoefL*gMaMa*
		(
		(inv_mu*mu1*t1+inv_t*t1)*dEicon_dt+
		alpha*u1+beta*w1
		);

	stab_matrix[4][3] = imagUnity*xi*
		(
		(inv_mu*mu1+vCoefL*inv_t)*(alpha*u1+beta*w1)+
		vCoefL*inv_t*inv_mu*mu1*t1*dEicon_dt
		);

	stab_matrix[5][3] = imagUnity*xi*vCoefL*inv_t*dEicon_dt;

	stab_matrix[6][3] = -imagUnity*xi*beta*
		(2.0*inv_mu*mu1*t1 + vCoefL*inv_t*t1);

	stab_matrix[7][3] = -imagUnity*xi*beta;
	// fifth
	stab_matrix[5][4]=1.0;
	// sixth
	stab_matrix[1][5] = -2.0*MF_Field::Pr*g_1MaMa*u1;

	stab_matrix[2][5] = MF_Field::Pr*
		(
		stabRe*inv_t*inv_mu*t1 - 
		2.0*imagUnity*g_1MaMa*(alpha*u1+beta*w1)
		);

	stab_matrix[3][5] = -imagUnity*stabRe*MF_Field::Pr*inv_mu*g_1MaMa*dEicon_dt;

	stab_matrix[4][5] = imagUnity*stabRe*MF_Field::Pr*inv_t*inv_mu*dEicon_dt +
		alpha*alpha + beta*beta - 
		g_1MaMa*MF_Field::Pr*inv_mu*mu1*(u1*u1+w1*w1)-
		inv_mu*(mu2*t1*t1+mu1*t2);
	stab_matrix[5][5] = -2.0*inv_mu*mu1*t1;

	stab_matrix[7][5] = -2.0*MF_Field::Pr*g_1MaMa*w1;
	// seventh
	stab_matrix[7][6]=1.0;
	// last
	stab_matrix[2][7] = -imagUnity*beta*(inv_mu*mu1*t1+vCoefM*inv_t*t1)+
		stabRe*inv_mu*inv_t*w1;

	stab_matrix[3][7] = imagUnity*stabRe*beta*inv_mu-
		vCoefM*beta*gMaMa*dEicon_dt;

	stab_matrix[4][7] = vCoefM*beta*inv_t*dEicon_dt-
		inv_mu*(mu2*t1*w1+mu1*w2);

	stab_matrix[5][7] = -inv_mu*mu1*w1;

	stab_matrix[6][7] = imagUnity*stabRe*inv_t*inv_mu*dEicon_dt +
		alpha*alpha + beta*beta;

	stab_matrix[7][7] = -inv_mu*mu1*t1;
	return stab_matrix;
};

t_Vec t_StabSolver::formRHS3D(const double& a_y, const t_Vec& a_vars) const{
	// TODO: check input vector dimension and
	// if it is not of 8 elems throw mega-exception)
	t_SqMatrix stab_matrix = getStabMatrix3D(a_y);
	t_Matrix input(1, 8);
	for (int i=0; i<8; i++){
		input[0][i] = a_vars[i];
	};
	// vars used in stability matrix 
	
	// after multiplication we have a rhs vector - matrix 1x8
	t_Matrix output = stab_matrix.mul(input); 
	return output[0];
};

t_Matrix t_StabSolver::getAsymptotics3D(const t_WaveChars& a_waveChars) const{
	t_Matrix initial_vectors(4,8);
	t_SqMatrix b_coef(4);
	t_Vec lambda(4,0.0);
	const double& y_e = _profStab.y.back();
	// TODO: function for simplified asymp: u=1.0, u'=0, u''=0, ... ?
	t_SqMatrix outer_matrix = getStabMatrix3D(y_e);	

	// TODO: what is all this about?

	b_coef[0][0]=outer_matrix[0][1];
	b_coef[1][0]=outer_matrix[3][1];
	b_coef[2][0]=outer_matrix[4][1];
	for (int i=1; i<3; i++){
		b_coef[i][1] =  outer_matrix[i+2][1]*outer_matrix[1][3]+
						outer_matrix[i+2][2]*outer_matrix[2][3]+
						outer_matrix[i+2][5]*outer_matrix[5][3]+
						outer_matrix[i+2][7]*outer_matrix[7][3];
		b_coef[i][2] = outer_matrix[i+2][5];
		b_coef[i][3] = outer_matrix[i+2][7];
	};
	b_coef[3][3] = outer_matrix[0][1];
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
		t_CompVal l2 = std::pow(lambda[i],2);
		t_CompVal denom = outer_matrix[0][1] - l2;
		b_coef[0][i] = (
							(l2 - outer_matrix[4][5])*outer_matrix[3][1]+
							outer_matrix[4][1]*outer_matrix[3][5]
					   )/denom;

		b_coef[1][i] = outer_matrix[4][5] - l2;
		b_coef[2][i] = -outer_matrix[3][5];
		b_coef[3][i] = (
							outer_matrix[3][5]*outer_matrix[4][7]+
							(l2 - outer_matrix[4][5])*outer_matrix[3][7]
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
									outer_matrix[0][2]*b_coef[0][i]+
									outer_matrix[3][2]*b_coef[1][i]+
									outer_matrix[4][2]*b_coef[2][i]+
									outer_matrix[6][2]*b_coef[3][i]
								)/lambda[i];

		initial_vectors[i][3] = b_coef[1][i];
		initial_vectors[i][4] = b_coef[2][i];
		initial_vectors[i][5] = lambda[i]*b_coef[2][i];
		initial_vectors[i][6] = b_coef[3][i];
		initial_vectors[i][7] = (
									outer_matrix[3][7]*b_coef[1][i]+
									outer_matrix[4][7]*b_coef[2][i]+
									outer_matrix[6][7]*b_coef[3][i]
								)/lambda[i];
	}
	return initial_vectors;
}


void t_StabSolver::set3DContext(const int& i_ind, const int& k_ind, const int& a_nnodesStab){
	_math_solver.set3DContext();
	_math_solver.resizeGrid(a_nnodesStab);

	t_ProfileNS profNS(_rFldNS);
	profNS.setProfiles(i_ind, k_ind);

	_profStab.resize(a_nnodesStab);
	_profStab.setProfiles(profNS);

	for (int j=0; j<a_nnodesStab; j++){
		_math_solver.varRange[j] = _profStab.y[a_nnodesStab-1-j];
	}

	return;
}

t_WaveChars t_StabSolver::searchMaxInstability(const t_WaveChars& initial_guess){
	// rewrite fortran |
	//SEARCH_MAX_INSTAB_TIME();
	t_WaveChars max_instab = initial_guess;
	return max_instab;
};
