#ifndef __PF_GS_PARAMS
#define __PF_GS_PARAMS

#include "PluginBase.h"
#include "mf_shared.h"

namespace pf{

	struct t_EigenGSParams{

		static hsstab::TPluginParamsGroup default_settings();
		static void init_base_params(t_EigenGSParams& params, const hsstab::TPluginParamsGroup& g);

		int NVars, NNodes;
		double ThickCoef;
		double W_Threshold;
		int NSProfInit;
	};

	namespace gs{
		void _init_eigen_gs_base_params(t_EigenGSParams& params, const hsstab::TPluginParamsGroup& g);
		void _eigen_gs_default_settings(hsstab::TPluginParamsGroup& g);
	};

};

#endif // __PF_GS_PARAMS