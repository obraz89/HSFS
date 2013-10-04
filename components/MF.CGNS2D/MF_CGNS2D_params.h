#ifndef __MFCGNS2D_PARAMS
#define __MFCGNS2D_PARAMS

#include "PluginBase.h"
#include "mf_shared.h"

#include "wx/fileconf.h"

namespace mf{

	class t_CGNS2DParams: public t_FldParams{

		int _Nz;

		void _init_params_map();
	public:
		t_CGNS2DParams():t_FldParams(){};
		enum{AxeSym=0, Plane} MFSym;
		double ZSpan;
	};

	// To reduce code
	namespace cgns2d{

		//tmp, read later from cgns domain 
		void _init_fld_base_params(t_FldParams& params, hsstab::TPluginParamsGroup& g);
		void _plug_default_settings(hsstab::TPluginParamsGroup& g);
	};

}		// ~namespace hsflow

#endif  // __MFCGNS2D_PARAMS