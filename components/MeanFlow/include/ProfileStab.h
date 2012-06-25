#ifndef __PROFILE_STAB
#define __PROFILE_STAB
#include "ProfileNS.h"

struct t_StabScales{
	double ReStab, Me;
	// dels is dimensional y_scale: 
	// sqrt(nu_e*x_e/u_e)
	// *_e are all dimensional in the formula
	double Dels; 
	// dimensional and nondim ue to pass to wave chars characteristics
	double UeDim;
	double Ue;
	double FreqScale() const{return UeDim/Dels;};
};

class  t_ProfileStab : public t_Profile{	
	t_StabScales _scales;
public:
	inline const t_StabScales& scales() const{return _scales;};
	t_ProfileStab(const t_MeanFlow& a_rFld);
	void initialize(t_ProfileNS& a_rProfNS, int nnodes=0);
	void initialize(int a_i, int a_k, int nnodes);
	void initialize(int a_i, int a_k);
	inline t_SqMat3 getJac() const{return _jacToLocalRF;};
	~t_ProfileStab();
};

#endif // __SM_PROFILE