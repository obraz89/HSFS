#ifndef __MFCGNS2D_PARAMS
#define __MFCGNS2D_PARAMS

#include "PluginBase.h"
#include "mf_shared.h"

#include "wx/fileconf.h"

namespace mf{

	struct t_CGNS2DParams: public t_FldParams{

		typedef std::map<wxString,int> t_MapWxStrInt; 

		static t_MapWxStrInt AXESYM_MODES_STR;

		static void init_supported_options();

		t_CGNS2DParams():t_FldParams(){};
		int Nz;
		t_AxeSym MFSym;
		double ZSpan;
		double ThetaSpan;

		static void plug_default_settings(hsstab::TPluginParamsGroup& g);
		static void init_fld_base_params(t_CGNS2DParams& params, const hsstab::TPluginParamsGroup& g);
	};

	// To reduce code
	namespace cg{
		namespace hsf2d{
		}
	};

}		// ~namespace hsflow

#endif  // __MFCGNS2D_PARAMS