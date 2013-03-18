#ifndef __PF_LS_PARAMS
#define __PF_LS_PARAMS

#include "PluginBase.h"

namespace pf{

	struct t_StabSolverParams{

		static hsstab::TPluginParamsGroup default_settings();
		static void init_base_params(t_StabSolverParams& params, const hsstab::TPluginParamsGroup& g);

		int NVars, NNodes;
		double ThickCoef;
		double AdjustTol, AdjustStep;
		int AdjustMaxIter;
	};

	namespace ls{
		void _init_ortho_ls_base_params(t_StabSolverParams& params, const hsstab::TPluginParamsGroup& g);
		void _ortho_ls_default_settings(hsstab::TPluginParamsGroup& g);
	}

}

#endif	// __PF_LS_PARAMS