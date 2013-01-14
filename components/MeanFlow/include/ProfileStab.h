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

	t_ProfileStab();

	void initialize(t_ProfileNS& a_rProfNS, int nnodes=0);

	// for testing with AVF code
	void initialize(const std::wstring fname, const t_StabScales& a_scales);

	~t_ProfileStab();
};

#endif // __SM_PROFILE