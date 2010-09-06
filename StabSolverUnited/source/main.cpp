#include "TaskParameters.h"
#include "MF_Field.h"
#include "SolverCore.h"
#include "SmProfile.h"
#include "Smooth.h"
#include <iostream>
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
	field.trans_to_cyl();
	SmProfile cur_profile(field, nx-2, nz/2);
	cur_profile.smooth();
	// set "global" vars for stability calculations
	cur_profile.setSolverParameters();
	SEARCH_MAX_INSTAB_TIME();
	getchar();
	return 0;
}