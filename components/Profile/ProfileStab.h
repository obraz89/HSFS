#ifndef __PROFILE_STAB
#define __PROFILE_STAB

#include "stab_shared.h"

#include "ProfileNS.h"



class IMPEXP_PROFILE t_ProfileStab : public t_Profile{	

	t_StabScales _scales;

public:

	const t_StabScales& scales() const;

	t_ProfileStab();

	void initialize(t_ProfileNS& a_rProfNS, int nnodes=0);

	// for testing with AVF code
	void initialize(const std::wstring fname, const t_StabScales& a_scales);

	~t_ProfileStab();
};

#endif // __SM_PROFILE