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

	void _initialize_extract(const mf::t_GeomPoint& xyz, const mf::t_ProfDataCfg& data_cfg);
	void _initialize_interpolate(const mf::t_GeomPoint& xyz, const mf::t_ProfDataCfg& data_cfg);

public:

	t_ProfileNS(const mf::t_DomainBase& rDomain);

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

};
#endif // __PROF_NS