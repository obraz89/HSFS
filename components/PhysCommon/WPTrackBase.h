#ifndef __WPTRACK_BASE
#define __WPTRACK_BASE

#include "PluginBase.h"
#include "MFBlockBase.h"
#include "LocSearchBase.h"

#include "WaveChars.h"

#include "dll_impexp-phys_common.h"

namespace stab{

	class IMPEXP_PHYSCOMMON t_WPTrackBase: public hsstab::TPlugPhysPart{
	public:

		t_WPTrackBase();

		virtual void retrace(mf::t_BlkInd start_from, t_WCharsLoc init_wave, 
			stab::t_LSBase& loc_solver)=0;

		virtual ~t_WPTrackBase();
	};

};

#endif	// __WPTRACK_BASE