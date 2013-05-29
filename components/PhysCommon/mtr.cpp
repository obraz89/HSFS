#include "stdafx.h"
#include "mf_shared.h"
#include "MFBlockBase.h"

using namespace mf;

// for future
t_Mtr t_Block::get_mtr(const t_BlkInd ind) const{

	// HAck, works only for a flat plate
	t_Mtr mtr;

	mtr.jac.setToUnity();

	return mtr;
};