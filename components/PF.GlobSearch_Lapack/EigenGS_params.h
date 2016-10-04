#ifndef __PF_GS_PARAMS
#define __PF_GS_PARAMS

#include "PluginBase.h"
#include "mf_shared.h"

#include "ProfileStab.h"

typedef std::map<wxString,int> t_MapWxStrInt; 

namespace pf{

	struct t_EigenGSParams{

		int NVars, NNodes;
		double ThickCoef;
		double ThickHalfNodesCoef;
		double Arg_Threshold;
		double SecondViscRatio;

		static t_MapWxStrInt PROFNS_INIT_TYPES_STR;
		blp::t_NSInit NSProfInit;

		static t_MapWxStrInt PROFSTAB_NONDIM_TYPES_STR;
		t_ProfStabCfg::t_Nondim NondimScaleType;

		static void init_supported_options();

		void init(const hsstab::TPluginParamsGroup& g);

		static void default_settings(hsstab::TPluginParamsGroup& g);
	};

};

#endif // __PF_GS_PARAMS