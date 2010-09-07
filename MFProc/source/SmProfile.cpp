#include "SmProfile.h"
#include "Smooth.h"
#include "SolverCore.h"
#include <cmath>
SmProfile::SmProfile(const MF_Field &__fld_ref, int __i_ind, int __k_ind)
:fld_ref(__fld_ref), i_ind(__i_ind), k_ind(__k_ind), ny(__fld_ref.ny){
	y_prof = new double[ny];
	u_prof = new double[ny];
	u1_prof = new double[ny];
	u2_prof = new double[ny];
	w_prof = new double[ny];
	w1_prof = new double[ny];
	w2_prof = new double[ny];
	t_prof = new double[ny];
	t1_prof = new double[ny];
	t2_prof = new double[ny];
	p_prof = new double[ny];
	r_prof = new double[ny];
	__fld_ref.get_profiles(__i_ind, __k_ind, y_prof, u_prof, w_prof, p_prof, t_prof, r_prof, u_e, w_e);
};

SmProfile::~SmProfile(){
	delete[] y_prof, u_prof, u1_prof, u2_prof, 
		     w_prof, w1_prof, w2_prof,
			 t_prof, t1_prof, t2_prof,
			 p_prof, r_prof;
}

void SmProfile::smooth(){
	int size = ny;
	SMOOTH_3D_PROFILES(y_prof, u_prof, &size, u1_prof, u2_prof);
	SMOOTH_3D_PROFILES(y_prof, w_prof, &size, w1_prof, w2_prof);
	SMOOTH_3D_PROFILES(y_prof, t_prof, &size, t1_prof, t2_prof);
}
void SmProfile::setSolverParameters(){
	DNS.AMINF = MF_Field::Mach;
	DNS.REINF = MF_Field::Re;
	DNS.TINF = MF_Field::T_inf;
	DNS.XLL = MF_Field::L_ref;

	NS.XST[0] = fld_ref.fld[i_ind][0][k_ind].x;
	CF_PROFILE_DATA.MF_WE = w_e;
	CF_PROFILE_DATA.MF_ANG = fld_ref.get_cf_wave_dir(i_ind, k_ind);
	ADDITIONAL_NS.Y_DIM = ny;

	for (int i=0; i<this->ny;i++){
		NS.YNS[i] = y_prof[i]/sqrt(NS.XST[0]);
		NS.UNS[i] = u_prof[i];
		NS.UNS1[i] = u1_prof[i]*sqrt(NS.XST[0]);
		NS.UNS2[i] = u2_prof[i]*NS.XST[0];
		NS.WNS[i] = w_prof[i];
		NS.WNS1[i] = w1_prof[i]*sqrt(NS.XST[0]);
		NS.WNS2[i] = w2_prof[i]*NS.XST[0];
		NS.TNS[i] = t_prof[i];
		NS.TNS1[i] = t1_prof[i]*sqrt(NS.XST[0]);
		NS.TNS2[i] = t2_prof[i]*NS.XST[0];
		RHONS.RHO[i] = r_prof[i];
	}
}

void SmProfile::searchMaxInstability(){
	SEARCH_MAX_INSTAB_TIME();
}
//void SmProfile::globalTS(){
//	TS_GLOBAL_TIME();
//};