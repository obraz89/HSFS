#include "StabSolver.h"
//#include "SolverCore.h"
t_StabSolver::t_StabSolver(const MF_Field& a_rFldNS, t_StabField& a_rFldStab, t_ODES& a_odes):
_rFldNS(a_rFldNS), _rFldStab(a_rFldStab), _odes(a_odes), _profStab(0){};
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
t_Vec t_StabSolver::rhsStabMatrix3D(const double& a_y, const t_Vec& a_vars){
	// TODO: check input vector dimension and
	// if it is not of 8 elems throw mega-exception)
	t_Matrix input(1, 8);
	for (int i=0; i<8; i++){
		input[0][i] = a_vars[i];
	};
	// vars used in stability matrix 
	double y = a_y;
	//double u = -_profStab.
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

	t_SqMatrix stab_matrix(8);
// first
	stab_matrix[2][1]=1.0;
// second
	stab_matrix[1][2] = imagUnity*stabRe*inv_t*inv_mu*dEicon_dt +
					   pow(alpha,2) + pow(beta,2);

    stab_matrix[2][2] =  - inv_mu*mu1*t1;

	stab_matrix[3][2] = stabRe*inv_t*inv_mu*u1 - 
						imagUnity*alpha*(inv_mu*mu1*t1+vCoefM*inv_t*t1);

	stab_matrix[4][2] = alpha*(imagUnity*stabRe*inv_mu -
						vCoefM*gMaMa*dEicon_dt);

	stab_matrix[5][2] = vCoefM*alpha*inv_t*dEicon_dt - inv_mu*dMu2;

	stab_matrix[6][2] = -inv_mu*mu1*u1;
// third
	stab_matrix[1][3] = -imagUnity*alpha;
	
	stab_matrix[3][3] = inv_t*t1;

	stab_matrix[4][3] = -imagUnity*gMaMa*dEicon_dt;

	stab_matrix[5][3] = imagUnity*inv_t*dEicon_dt;

	stab_matrix[7][3] = -imagUnity*beta;
// fourth
	stab_matrix[1][4] = -imagUnity*xi*alpha*t1*
						(2.0*inv_mu*mu1 + vCoefL*inv_t);

	stab_matrix[2][4] = -imagUnity*xi*alpha;

	stab_matrix[3][4] = xi*(-alpha*alpha-beta*beta+
							vCoefL*inv_t*inv_mu*mu1*t1*t1+
							vCoefL*inv_t*t2-
							imagUnity*stabRe*inv_t*inv_mu*dEicon_dt);

	stab_matrix[4][4] = -imagUnity*xi*vCoefL*gMaMa*
						 (
							(inv_mu*mu1*t1+inv_t*t1)*dEicon_dt+
							alpha*u1+beta*w1
						 );
	
	stab_matrix[5][4] = imagUnity*xi*
						(
							(inv_mu*mu1+vCoefL*inv_t)*(alpha*u1+beta*w1)+
							vCoefL*inv_t*inv_mu*mu1*t1*dEicon_dt
						);

	stab_matrix[6][4] = imagUnity*xi*vCoefL*inv_t*dEicon_dt;
	
	stab_matrix[7][4] = -imagUnity*xi*beta*
						 (2.0*inv_mu*mu1*t1 + vCoefL*inv_t*t1);
	
	stab_matrix[8][4] = -imagUnity*xi*beta;
// fifth
    stab_matrix[6][5]=1.0;
// sixth
	stab_matrix[2][6] = -2.0*MF_Field::Pr*g_1MaMa*u1;

	stab_matrix[3][6] = MF_Field::Pr*
						(
							stabRe*inv_t*inv_mu*t1 - 
							2.0*imagUnity*g_1MaMa*(alpha*u1+beta*w1)
						);

	stab_matrix[4][6] = -imagUnity*stabRe*MF_Field::Pr*inv_mu*g_1MaMa*dEicon_dt;

	stab_matrix[5][6] = imagUnity*stabRe*MF_Field::Pr*inv_t*inv_mu*dEicon_dt +
						alpha*alpha + beta*beta - 
						g_1MaMa*MF_Field::Pr*inv_mu*mu1*(u1*u1+w1*w1)-
						inv_mu*(mu2*t1*t1+mu1*t2);
	stab_matrix[6][6] = -2.0*inv_mu*mu1*t1;
	
	stab_matrix[8][6] = -2.0*MF_Field::Pr*g_1MaMa*w1;
// seventh
	stab_matrix[8][7]=1.0;
// last
	stab_matrix[3][8] = -imagUnity*beta*(inv_mu*mu1*t1+vCoefM*inv_t*t1)+
						 stabRe*inv_mu*inv_t*w1;

	stab_matrix[4][8] = imagUnity*stabRe*beta*inv_mu-
						vCoefM*beta*gMaMa*dEicon_dt;

	stab_matrix[5][8] = vCoefM*beta*inv_t*dEicon_dt-
						inv_mu*(mu2*t1*w1+mu1*w2);

	stab_matrix[6][8] = -inv_mu*mu1*w1;
	
	stab_matrix[7][8] = imagUnity*stabRe*inv_t*inv_mu*dEicon_dt +
						alpha*alpha + beta*beta;

	stab_matrix[8][8] = -inv_mu*mu1*t1;
	// after multiplication we have a rhs vector - matrix 1x8
	t_Matrix output = stab_matrix.mul(input); 
	return output[0];
};

void t_StabSolver::set3DContext(const int& i_ind, const int& k_ind, const int& a_nnodesNS, const int& a_nnodesStab){
	t_ProfileNS profNS(_rFldNS, a_nnodesNS);
	profNS.setProfiles(i_ind, k_ind);
	t_ProfileStab profStab(a_nnodesStab);
	profStab.setProfiles(profNS);
	_profStab = profStab;
	// TODO: set odes context ??? 2D vs 3D

	return;
}

t_WaveChars t_StabSolver::searchMaxInstability(const t_WaveChars& intial_guess){
	// rewrite fortran |
	//SEARCH_MAX_INSTAB_TIME();
};
void t_StabSolver::searchGlobal(){
	//SEARCH_INITIAL_INSTAB_TIME();
	double w = 0.06;
	double phs = 0.8;
	double angle = 0.0;
	//GLOBAL_TIME_GRAD(w, phs, angle);
};