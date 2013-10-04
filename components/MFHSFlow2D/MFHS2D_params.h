#ifndef __HSMF2D_PARAMS
#define __HSMF2D_PARAMS

#include "PluginBase.h"
#include "mf_shared.h"

#include "wx/fileconf.h"

namespace mfhs{

	class t_HSFlowParams2D: public mf::t_FldParams{
		void _init_params_map();
	public:
		t_HSFlowParams2D():t_FldParams(){};
		mf::t_AxeSym MFSym;
		double ZSpan;
	};

	// To reduce code
	namespace hsf2d{
		void _init_fld_base_params(mf::t_FldParams& params, const hsstab::TPluginParamsGroup& g);
		void _hsflow_default_settings(hsstab::TPluginParamsGroup& g);
	};

}		// ~namespace mfhs

#endif  // __HSMF_PARAMS