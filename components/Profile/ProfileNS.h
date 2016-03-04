#ifndef __PROF_NS
#define __PROF_NS

#include "dll_impexp-profile.h"
#include "Profile.h"


/************************************************************************/
//
// Profiles extracted from t_MeanFlow structure
// this profile is in external streamline reference frame
// vertical coord is multipied by sqrt(Re)
// t_ProfileNS <=> t_ProfMFLoc (lazy to rewrite)
/************************************************************************/

class IMPEXP_PROFILE t_ProfileNS : public t_ProfMF{

	// true if this is a disturbance profile
	bool _is_disturb_profile;

	void _initialize_extract(const mf::t_GeomPoint& xyz, const mf::t_ProfDataCfg& data_cfg);
	void _initialize_interpolate(const mf::t_GeomPoint& xyz, const mf::t_ProfDataCfg& data_cfg);

public:

	t_ProfileNS(const mf::t_DomainBase& rDomain);

	bool is_disturbance_profile(){return _is_disturb_profile;};

	// tmp, to load DNS profiles
	void init_from_uvwpt_vec(const std::vector<double>& y_vec, 
		const std::vector<std::vector<double>>& v_uvwpt,
		double a_bl_thick_scale);

};

/************************************************************************/
//
// Profiles extracted from t_MeanFlow structure
// this profile is in Global Reference Frame
// i.e. profile extracted from MF Domain 
// with no additional transformations
/************************************************************************/

class IMPEXP_PROFILE t_ProfMFGlob : public t_ProfMF{

	void _initialize_extract(const mf::t_GeomPoint& xyz, const mf::t_ProfDataCfg& data_cfg);
	void _initialize_interpolate(const mf::t_GeomPoint& xyz, const mf::t_ProfDataCfg& data_cfg);

public:

	t_ProfMFGlob(const mf::t_DomainBase& rDomain);

};

/************************************************************************/
//
// Profiles extracted from t_MeanFlow structure
// this profile is in Local Reference Frame
// i.e. profile extracted from MF Domain 
// with vectors transformed to Local RF
/************************************************************************/

class IMPEXP_PROFILE t_ProfMFLoc : public t_ProfMF{

	void _initialize_extract(const mf::t_GeomPoint& xyz, const mf::t_ProfDataCfg& data_cfg);
	void _initialize_interpolate(const mf::t_GeomPoint& xyz, const mf::t_ProfDataCfg& data_cfg);

public:

	t_ProfMFLoc(const mf::t_DomainBase& rDomain);

	void dump(const std::string& fname) const;

};
#endif // __PROF_NS