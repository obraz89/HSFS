#ifndef __PF_LS_PARAMS
#define __PF_LS_PARAMS

#include "PluginBase.h"

#include "ProfileStab.h"

typedef std::map<wxString,int> t_MapWxStrInt; 
typedef std::map<wxString, bool> t_MapWxStrBool;

namespace pf{

	struct t_StabSolverParams{

		int NVars, NNodes;
		double ThickCoef;
		double AdjustTol, AdjustStep;
		int AdjustMaxIter;

		t_MapWxStrBool CURV_TERMS_FLAGS_STR;
		bool CurvTermsOn;

		t_MapWxStrInt PROFNS_INIT_TYPES_STR;
		blp::t_NSInit NSProfInit;

		t_MapWxStrInt PROFSTAB_NONDIM_TYPES_STR;
		t_ProfStabCfg::t_Nondim NondimScaleType;

		double MaxNonDimIncrement;

		bool bCheckAlphaPositive;
		bool bCheckPhaseSpeedSSonic;

		enum t_WallBC {
			WALL_HOMOGEN = 0,
			WALL_SLIP
		};
		t_WallBC WallBC;
		t_MapWxStrInt WALL_BCS_OPTS_STR;

		double WallBC_EtaU, WallBC_EtaW;

		t_StabSolverParams();
		void init_supported_options();

		void init(const hsstab::TPluginParamsGroup& g);

		static void default_settings(hsstab::TPluginParamsGroup& g);

	};

}

#endif	// __PF_LS_PARAMS