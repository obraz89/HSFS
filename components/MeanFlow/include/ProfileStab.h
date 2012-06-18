#ifndef __PROFILE_STAB
#define __PROFILE_STAB
#include "ProfileNS.h"

class  t_ProfileStab : public t_Profile{	
public:
	double stabRe, Me;
	// dels is dimensional y_scale: 
	// sqrt(nu_e*x_e/u_e)
	// *_e are all dimensional in the formula
	double dels; 
	// dimensional u_e to pass to wave chars characteristics
	double ue_dim;
	t_ProfileStab(const t_MeanFlow& a_rFld);
	void initialize(t_ProfileNS& a_rProfNS, int nnodes=0);
	void initialize(int a_i, int a_k, int nnodes);
	void initialize(int a_i, int a_k);
	inline t_SqMat3 getJac() const{return _jacToLocalRF;};
	~t_ProfileStab();
};

#endif // __SM_PROFILE