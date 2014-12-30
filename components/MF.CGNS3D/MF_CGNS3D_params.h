#ifndef __MFHS3D_PARAMS
#define __MFHS3D_PARAMS

#include "PluginBase.h"
#include "mf_shared.h"

#include "wx/fileconf.h"

namespace mf{

	class t_CGNS3DParams: public t_FldParams{

		void _init_params_map();
	public:
		t_CGNS3DParams():t_FldParams(){};
	};

	namespace cg{
		namespace hsf3d{
			void _init_fld_base_params(t_FldParams& params, const hsstab::TPluginParamsGroup& g);
			void _plug_default_settings(hsstab::TPluginParamsGroup& g);
		}
	}

}		// ~namespace hsflow

#endif  // __MFHS3D_PARAMS