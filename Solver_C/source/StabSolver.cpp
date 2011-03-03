#include "StabSolver.h"
#include "SolverCore.h"
t_StabSolver::t_StabSolver(const MF_Field& a_rFldNS, t_StabField& a_rFldStab, t_ODES& a_odes):
_rFldNS(a_rFldNS), _rFldStab(a_rFldStab), _odes(a_odes){};
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
t_Vec t_StabSolver::rhs3DStabMatrix(const double& a_y, const t_Vec& a_vars){
	// TODO: check input vector dimension and
	// if it is not of 8 elems throw mega-exception)
	t_Matrix input(1, 8);
	for (int i=0; i<8; i++){
		input[0][i] = a_vars[i];
	};
	t_SqMatrix stab_matrix(8);
	// first row, remember that all matrices are column-driven
	stab_matrix[2][1]=1.0;
	// second
	/*
	stab_matrix[1,2] = E*R*TD*M1*WA + A*A + B*B
    stab_matrix[2][2] =  - M1*MY
    stab = R*TD*M1*U1 - E*A*(M1*MY+F*TD*T1)
    H(4,2) = E*A*R*M1 - F*M2*A*WA
    H(5,2) = F*A*TD*WA - M1*UM
    H(6,2) =  - M1*MU1*U1*/

	// after multiplication we have a rhs vector - matrix 1x8
	t_Matrix output = stab_matrix.mul(input); 
	return output[0];
};

void t_StabSolver::set3DContext(t_ProfileStab &a_profStab){

}

void t_StabSolver::searchMaxInstability(){
	SEARCH_MAX_INSTAB_TIME();
};
void t_StabSolver::searchGlobal(){
	//SEARCH_INITIAL_INSTAB_TIME();
	double w = 0.06;
	double phs = 0.8;
	double angle = 0.0;
	GLOBAL_TIME_GRAD(w, phs, angle);
};