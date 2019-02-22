#ifndef __MFHS3D_PARAMS
#define __MFHS3D_PARAMS

#include "PluginBase.h"
#include "mf_shared.h"

#include "wx/fileconf.h"

typedef std::map<wxString, int> t_MapWxStrInt;

namespace mf{

	struct t_CGNS3DParams: public t_FldParams{

		t_MapWxStrInt VD_TYPES_STR;
		t_MapWxStrInt VD_PLACES_STR;

		cg::t_VDParams vd_params;

		t_CGNS3DParams();
		static void plug_default_settings(hsstab::TPluginParamsGroup& g);
		static void init_fld_base_params(t_CGNS3DParams& params, const const hsstab::TPluginParamsGroup& g);
	};


}		// ~namespace hsflow

#endif  // __MFHS3D_PARAMS