#include "StabSolver.h"
#include "SolverCore.h"
StabSolver::StabSolver(const MF_Field& _fld_ref, int _i, int _k):
fld_ref(_fld_ref), i_ind(_i), k_ind(_k), profile(_fld_ref, _i, _k){};
void StabSolver::setParameters(){
	// globals -flow
	DNS.AMINF = MF_Field::Mach;
	DNS.REINF = MF_Field::Re;
	DNS.TINF = MF_Field::T_inf;
	DNS.XLL = MF_Field::L_ref;
	// -gas parameters
	OUT.G = 1.4;	// specific heat ratio
	OUT.SIG = 0.72;	// Prandtl number
	IUPT.IUPT = 1;  // for viscosity: iupt=0 - Sutherland, 1 - power
	BASA6.BVB = 0.75; // if iupt=1, BVB - viscosity power coef
	OUT.K = 0.0;	  // second viscosity
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
};

void StabSolver::smoothProfile(){
	profile.smooth();
};

void StabSolver::adaptProfile(){
	profile.adapt();
};

void StabSolver::searchMaxInstability(){
	SEARCH_MAX_INSTAB_TIME();
};
void StabSolver::searchGlobal(){
	SEARCH_INITIAL_INSTAB_TIME();
};