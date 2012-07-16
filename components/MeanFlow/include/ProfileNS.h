#ifndef __PROF_NS
#define __PROF_NS
#include "MeanFlow.h"
#include "Profile.h"
#include "smooth.h"



class  t_ProfileNS : public t_Profile{
	friend class t_ProfileStab;
	int _i, _k, _bl_bound_ind;
public:
	//int iMF, kMF;
	void initialize(int a_i, int a_k, double a_thick_coef);
	t_Rec get_bl_bound_rec();
	double xDist;
	t_ProfileNS(const t_MeanFlow& rFld);
	~t_ProfileNS();
};
#endif // __PROF_NS