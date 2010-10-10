#include "MF_Field.h"
#ifndef __SMPROFILE
#define __SMPROFILE
class SmProfile {
public:
	int i_ind, k_ind;
	const int ny;
	double *y_prof, *u_prof, *u1_prof, *u2_prof,
					*w_prof, *w1_prof, *w2_prof,
					*t_prof, *t1_prof, *t2_prof,
					*p_prof, *r_prof;
	double x_st, u_e, w_e, cf_angle;
	const MF_Field& fld_ref;

	SmProfile( const MF_Field& _fld_ref, int i_ind, int k_ind);
	~SmProfile();
	void smooth();	// smooth profiles and get derivatives using fortran imsl lib
	void setSolverParameters();	// couples DataPoint with StabilitySolver
	void searchMaxInstability(); // wraps SEARCH_MAX_INSTAB_TIME
	//void globalTS();
};
#endif // __SMPROFILE