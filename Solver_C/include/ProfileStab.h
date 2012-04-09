#include "ProfileNS.h"
#ifndef __PROFILE_STAB
#define __PROFILE_STAB

class t_ProfileStab : public t_Profile{	
public:
	double stabRe, Me;
	// dels is dimensional y_scale: 
	// sqrt(nu_e*x_e/u_e)
	// *_e are all dimensional in the formula
	double dels; 
	t_ProfileStab(const t_MeanFlow& a_rFld);
	void initialize(t_ProfileNS& a_rProfNS, int nnodes=0);
	inline t_SqMat3 getJac() const{return _jacToLocalRF;};
	~t_ProfileStab();
};

#endif // __SM_PROFILE