#ifndef __MFCGNS2D_PARAMS
#define __MFCGNS2D_PARAMS

#include "PluginBase.h"
#include "mf_shared.h"

#include "wx/fileconf.h"

namespace mf{

	class t_CGNS2DParams: public t_FldParams{

		void _init_params_map();
	public:
		t_CGNS2DParams():t_FldParams(){};
		int Nz;
		t_AxeSym MFSym;
		double ZSpan;
		double ThetaSpan;
	};

	// To reduce code
	namespace cg{
		namespace hsf2d{
			void _init_fld_base_params(t_FldParams& params, const hsstab::TPluginParamsGroup& g);
			void _plug_default_settings(hsstab::TPluginParamsGroup& g);
		}
	};

}		// ~namespace hsflow

#endif  // __MFCGNS2D_PARAMS