#ifndef __PROF_NS
#define __PROF_NS

#include "dll_impexp-profile.h"
#include "Profile.h"

#include "mf_shared.h"
#include "MFDomainBase.h"



// TODO: wrap all profiles in blp namespace {blp = Boundary Layer Profile}
namespace blp{
	enum t_NSInit{EXTRACT=0, INTERPOLATE};
}


/************************************************************************/
//
// Profiles extracted from t_MeanFlow structure
/************************************************************************/

class IMPEXP_PROFILE t_ProfileNS : public t_Profile{

	mf::t_GeomPoint _xyz;

	// TODO:for now this is not used, may be used later?
	int __bl_bound_ind;

	// TODO: remove?
	const mf::t_DomainBase& _rDomain;

	double _xScale;

	void _initialize_extract(const mf::t_GeomPoint& xyz, const mf::t_ProfDataCfg& init_cfg );
	void _initialize_interpolate(const mf::t_GeomPoint& xyz, const mf::t_ProfDataCfg& init_cfg );

public:

	t_ProfileNS(const mf::t_DomainBase& rDomain);

	void initialize(const mf::t_GeomPoint& xyz, const mf::t_ProfDataCfg& init_cfg, blp::t_NSInit init_type);

	const mf::t_DomainBase& getMFDomain() const;

	double get_x_Scale() const;

	int get_bound_ind() const;

	t_Rec get_bound_rec();

};
#endif // __PROF_NS