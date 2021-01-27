#ifndef __PROFILE_STAB
#define __PROFILE_STAB

#include "stab_shared.h"

#include "ProfileNS.h"

/************************************************************************/
//
// configure profile stab initialization
// NondimScaleType - use bl thick scale from cfd or use self-similar scale
/************************************************************************/
struct IMPEXP_PROFILE t_ProfStabCfg{

	int NNodes;

	enum t_Nondim{
		NONDIM_BY_BL_BOUND_SCALE=0,
		NONDIM_BY_DISP_THICK,
		NONDIM_BY_X_SELFSIM,
		NONDIM_BY_FIXED_VAL
	};
	
	t_Nondim NondimScaleType;

	t_ProfStabCfg():NNodes(0), NondimScaleType(NONDIM_BY_BL_BOUND_SCALE){}

};

/************************************************************************/
//
// configure profile stab initialization
// for DNS disturbance profile
/************************************************************************/
struct IMPEXP_PROFILE t_ProfStabCfgDNSDisturb {

	int NNodes;

	double Ymax;

	double Dels_gndim;

	double Ue, Rhoe, Te;

};

class IMPEXP_PROFILE t_ProfileStab : public t_Profile{	

	t_StabScales _scales;

	void _initialize(t_ProfMFLoc& a_rProfNS, 
		const std::vector<double>& y_distrib, t_ProfStabCfg cfg);

public:

	const t_StabScales& scales() const;

	t_ProfileStab();

	void initialize(t_ProfMFLoc& a_rProfNS, t_ProfStabCfg cfg);

	void initialize(t_ProfMFLoc& a_rProfNS,
		const std::vector<double>& y_distrib, t_ProfStabCfg cfg);

	void initialize_dist_DNS(t_ProfMFLoc& a_rProfNS, t_ProfStabCfgDNSDisturb cfg);

	// for testing with AVF code
	// w=0
	void initialize_2D(const std::string& fname, const t_StabScales& a_scales);
	// w!=0
	void initialize_3D(const std::string& fname, const t_StabScales& a_scales);

	virtual void dump(const std::string& fname) const;

	~t_ProfileStab();
};

#endif // __SM_PROFILE