#include "TaskParameters.h"
#include "MF_Field.h"
#include "SolverCore.h"
#include "SmProfile.h"
#include "Smooth.h"
#include <iostream>
void set_solver_parameters();
void set_solver_profiles(double* y, double* u, double* u1, double* u2, 
						            double* w, double* w1, double* w2, 
									double* t, double* t1, double* t2, 
						            double* r, int size);
int main(){
	int nx=81, ny=161, nz=51;
	std::string NSFieldName;
	// std::cin>> NSFieldName;
	// now it is for al=1, initialize NS field struct
	NSFieldName = "input/new/07.85000.dat";		
	MF_Field field(NSFieldName,nx,ny,nz);
	double* y_prof = new double[ny];
	double* u_prof = new double[ny];
	double* u1_prof = new double[ny];
	double* u2_prof = new double[ny];
	double* w_prof = new double[ny];
	double* w1_prof = new double[ny];
	double* w2_prof = new double[ny];
	double* p_prof = new double[ny];
	double* t_prof = new double[ny];
	double* t1_prof = new double[ny];
	double* t2_prof = new double[ny];
	double* r_prof = new double[ny];
	double current_ue, current_we;
	field.trans_to_cyl();
	field.get_profiles(nx-2,nz/2, y_prof, u_prof, w_prof, p_prof, t_prof, r_prof, current_ue, current_we);
	SMOOTH_3D_PROFILES(y_prof, u_prof, &ny, u1_prof, u2_prof);
	SMOOTH_3D_PROFILES(y_prof, w_prof, &ny, w1_prof, w2_prof);
	SMOOTH_3D_PROFILES(y_prof, t_prof, &ny, t1_prof, t2_prof);

	// set "global" vars for stability calculations
	set_solver_parameters();
	getchar();
	delete[] u_prof, u1_prof, u2_prof,
		     w_prof, w1_prof, w2_prof,
			 t_prof, t1_prof, t2_prof,
			 p_prof, r_prof;
	return 0;
}

void set_solver_parameters(){
	DNS.AMINF = MF_Field::Mach;
	DNS.REINF = MF_Field::Re;
	DNS.TINF = MF_Field::T_inf;
	DNS.XLL = MF_Field::L_ref;
};

void set_solver_profiles(double* y, double* u, double* u1, double* u2, 
						            double* w, double* w1, double* w2, 
									double* t, double* t1, double* t2, 
						            double* r, int size){

	for (int i=0; i<size;i++){
		NS.YNS[i] = y[i];
		NS.UNS[i] = u[i];
		NS.UNS1[i] = u1[i];
		NS.UNS2[i] = u2[i];
		NS.WNS[i] = w[i];
		NS.WNS1[i] = w1[i];
		NS.WNS2[i] = w2[i];
		NS.TNS[i] = t[i];
		NS.TNS1[i] = t1[i];
		NS.TNS2[i] = t2[i];
		RHONS.RHO[i] = r[i];
		// also xst!
	}
	
};