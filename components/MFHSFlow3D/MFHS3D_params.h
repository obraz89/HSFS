#ifndef __MFHS3D_PARAMS
#define __MFHS3D_PARAMS

#include "PluginBase.h"
#include "mf_shared.h"

#include "wx/fileconf.h"

namespace mfhs{

	namespace hsf3d{
		void _init_fld_base_params(mf::t_FldParams& params, const hsstab::TPluginParamsGroup& g);
		void _hsflow_default_settings(hsstab::TPluginParamsGroup& g);
	}

}		// ~namespace hsflow

#endif  // __MFHS3D_PARAMS