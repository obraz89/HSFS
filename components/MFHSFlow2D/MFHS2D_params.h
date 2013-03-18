#ifndef __HSMF2D_PARAMS
#define __HSMF2D_PARAMS

#include "PluginBase.h"
#include "mf_shared.h"

#include "wx/fileconf.h"

namespace mf{

	class t_HSFlowParams2D: public t_FldParams{
		void _init_params_map();
	public:
		t_HSFlowParams2D():t_FldParams(){};
		class t_AxeSym: public t_Enum{
		public:
			static const int AxeSym/*=0*/, Plane;
			t_AxeSym(){_init_map_vals();set_value(AxeSym);};
			void operator=(const int& val){t_Enum::operator =(val);};
			bool operator==(const int& val) const{return t_Enum::operator ==(val);};
		protected:
			void _init_map_vals();
		};
		t_AxeSym MFSym;
		double ZSpan;
	};

	// To reduce code
	namespace hsf2d{
		void _init_fld_base_params(t_FldParams& params, const hsstab::TPluginParamsGroup& g);
		void _hsflow_default_settings(hsstab::TPluginParamsGroup& g);
	};

}		// ~namespace hsflow

#endif  // __HSMF_PARAMS