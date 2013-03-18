#ifndef __PROF_NS
#define __PROF_NS

#include "dll_impexp-profile.h"
#include "Profile.h"

#include "mf_shared.h"
#include "MFBlockBase.h"

/************************************************************************/
//
// Profiles extracted from t_MeanFlow structure
/************************************************************************/
class IMPEXP_PROFILE t_ProfileNS : public t_Profile{

	mf::t_BlkInd _mf_ind;

	int _bl_bound_ind;

	// to speed up: avoid multiple calls _rBlk.mtr(_mf_ind)
	mf::t_Mtr _mtr;

	// TODO: remove?
	const mf::t_Block& _rBlk;

	// TODO: better rename to x_scale
	double _xDist;

public:

	t_ProfileNS(const mf::t_Block& rBlk);

	void initialize(const mf::t_BlkInd ind, double a_thick_coef);

	const mf::t_Block& getBlk() const;

	mf::t_BlkInd getMFInd() const;

	double get_xDist() const;

	int get_bl_bound_ind() const;

	t_Rec get_bl_bound_rec();

	const mf::t_Mtr& get_mtr() const;

	t_SqMat3Dbl getJac() const;
};
#endif // __PROF_NS