#include "stdafx.h"

#include "WavePackLine_params.h"
#include "common_data.h"
#include "log.h"

using namespace pf;

// ---------------------------------------

static const double TIME_STEP_DEFAULT = 0.01;

typedef std::map<wxString,int> t_MapWxStrInt; 

static t_MapWxStrInt  RETRACE_MODES_STR;
wxString RETRACE_MODE_DEFAULT_STR;

static t_MapWxStrInt MARCH_OPTS_STR;
wxString MARCH_OPT_DEFAULT_STR;

struct t_InitSupportedOptions{
	t_InitSupportedOptions(){
		RETRACE_MODES_STR.insert(std::make_pair(_T("W_FIXED"), stab::t_WPRetraceMode::W_FIXED));
		RETRACE_MODES_STR.insert(std::make_pair(_T("WB_FIXED"), stab::t_WPRetraceMode::WB_FIXED));
		RETRACE_MODES_STR.insert(std::make_pair(_T("ENVELOPE"), stab::t_WPRetraceMode::ENVELOPE));

		RETRACE_MODE_DEFAULT_STR = _T("W_FIXED");

		MARCH_OPTS_STR.insert(std::make_pair(_T("GROUP_VELO"), t_WPLineParams::GROUP_VELO));
		MARCH_OPTS_STR.insert(std::make_pair(_T("STREAMLINE"), t_WPLineParams::STREAMLINE));
		
		MARCH_OPT_DEFAULT_STR = _T("GROUP_VELO");
		
	}
};

// initialize supported modes before entering any wptrack method
t_InitSupportedOptions init_supported_modes;

void pf::_wpline_default_settings(hsstab::TPluginParamsGroup& g){

	g.add("TimeStep", TIME_STEP_DEFAULT , _T("dr=V*dt, set dt"));

	g.add("RetraceMode", RETRACE_MODE_DEFAULT_STR, _T("Retrace Mode"));

	g.add("MarchAlong", MARCH_OPT_DEFAULT_STR, _T("Retrace direction"));

}

void pf::_init_wpline_base_params(t_WPLineParams& params, const hsstab::TPluginParamsGroup& g){

	params.TimeStep = g.get_real_param("TimeStep");

	wxString rmode_str = g.get_string_param("RetraceMode");

	t_MapWxStrInt::iterator it = RETRACE_MODES_STR.find(rmode_str);

	if (it==RETRACE_MODES_STR.end()) 
		ssuGENTHROW(_T("Unknown value provided for option RetraceMode!"));

	int rmode = RETRACE_MODES_STR.find(rmode_str)->second;

	switch (rmode)
	{

	case stab::t_WPRetraceMode::W_FIXED:
		params.RetraceMode = stab::t_WPRetraceMode::W_FIXED;
		break;

	case stab::t_WPRetraceMode::WB_FIXED:
		params.RetraceMode = stab::t_WPRetraceMode::WB_FIXED;
		break;

	default:
		ssuGENTHROW(_T("Retracing Mode not supported!"));
	}

	wxString rdir_str = g.get_string_param("MarchAlong");

	it = MARCH_OPTS_STR.find(rdir_str);

	if (it==MARCH_OPTS_STR.end()) 
		ssuGENTHROW(_T("Unknown value provided for option MarchAlong!"));

	int rdir = MARCH_OPTS_STR.find(rdir_str)->second;

	switch (rdir)
	{
	case t_WPLineParams::GROUP_VELO:
		params.RetraceDir = t_WPLineParams::GROUP_VELO;
		break;

	case t_WPLineParams::STREAMLINE:
		params.RetraceDir = t_WPLineParams::STREAMLINE;
		break;
	default:
		ssuGENTHROW(_T("Retracing Direction option not supported!"));
	}
}

// ~t_WPLineParams
// ---------------------------------------
