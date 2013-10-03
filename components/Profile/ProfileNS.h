#ifndef __PROF_NS
#define __PROF_NS

#include "dll_impexp-profile.h"
#include "Profile.h"

#include "mf_shared.h"
#include "MFDomainBase.h"

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

public:

	t_ProfileNS(const mf::t_DomainBase& rDomain);

	void initialize(const mf::t_GeomPoint xyz, double a_thick_coef, int nnodes=0);

	const mf::t_DomainBase& getMFDomain() const;

	double get_x_Scale() const;

	int get_bound_ind() const;

	t_Rec get_bound_rec();

};
#endif // __PROF_NS