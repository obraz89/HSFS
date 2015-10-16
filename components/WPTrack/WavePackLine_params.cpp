#include "stdafx.h"

#include "WavePackLine_params.h"
#include "common_data.h"
#include "log.h"

using namespace pf;

#include "wx/tokenzr.h"

// ---------------------------------------

static const double TIME_STEP_DEFAULT = 0.01;

typedef std::map<wxString,int> t_MapWxStrInt; 

t_MapWxStrInt  t_WPLineParams::RETRACE_MODES_STR;
#define RETRACE_MODE_DEFAULT_STR _T("W_FIXED")

t_MapWxStrInt t_WPLineParams::MARCH_OPTS_STR;
#define MARCH_OPT_DEFAULT_STR _T("GROUP_VELO")


void t_WPLineParams::init_supported_options(){
		RETRACE_MODES_STR.clear();
		RETRACE_MODES_STR.insert(std::make_pair(wxString(_T("W_FIXED")), stab::t_WPRetraceMode::W_FIXED));
		RETRACE_MODES_STR.insert(std::make_pair(wxString(_T("WB_FIXED")), stab::t_WPRetraceMode::WB_FIXED));
		RETRACE_MODES_STR.insert(std::make_pair(wxString(_T("WBRAD_FIXED")), stab::t_WPRetraceMode::WBRAD_FIXED));
		RETRACE_MODES_STR.insert(std::make_pair(wxString(_T("ENVELOPE")), stab::t_WPRetraceMode::ENVELOPE));

		MARCH_OPTS_STR.clear();
		MARCH_OPTS_STR.insert(std::make_pair(wxString(_T("GROUP_VELO")), t_WPLineParams::GROUP_VELO));
		MARCH_OPTS_STR.insert(std::make_pair(wxString(_T("STREAMLINE")), t_WPLineParams::STREAMLINE));
		MARCH_OPTS_STR.insert(std::make_pair(wxString(_T("FIXED_DIRECTION")), t_WPLineParams::FIXED_DIRECTION));
		
		
};

void t_WPLineParams::wpline_default_settings(hsstab::TPluginParamsGroup& g){

	init_supported_options();

	g.add("TimeStep", TIME_STEP_DEFAULT , _T("dr=V*dt, set dt"));

	g.add("RetraceMode", RETRACE_MODE_DEFAULT_STR, _T("Retrace Mode"));

	g.add("MarchAlong", MARCH_OPT_DEFAULT_STR, _T("Retrace direction"));

	g.add("RetraceVec", _T("1.000; 0.000; 0.000"), _T("Retrace vector [when FIXED_DIRECTION option chosen]"));

}

void t_WPLineParams::read_parse_retrace_vec(const hsstab::TPluginParamsGroup& g){

	wxString rvec_str = g.get_string_param("RetraceVec");

	wxArrayString wxNames = wxStringTokenize(rvec_str , _T(';'));

	if (wxNames.Count()!=3) wxLogError(_T("Retrace vector parse error: type in x;y;z"));

	bool ok = true;

	for (int i=0; i<3; i++) {

		wxString& rStr = wxNames[i];

		// trim from both left and right
		rStr.Trim(true);rStr.Trim(false);

		ok =  ok && rStr.ToDouble(&(RetraceVec[i]));

	}

	if (!ok) wxLogError(_T("Retrace vector parse error: failed to parse"));

}

void t_WPLineParams::init_wpline_base_params(t_WPLineParams& params, const hsstab::TPluginParamsGroup& g){

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

	case stab::t_WPRetraceMode::WBRAD_FIXED:
		params.RetraceMode = stab::t_WPRetraceMode::WBRAD_FIXED;
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

	case t_WPLineParams::FIXED_DIRECTION:
		params.RetraceDir = t_WPLineParams::FIXED_DIRECTION;
		params.read_parse_retrace_vec(g);
		break;
	default:
		ssuGENTHROW(_T("Retracing Direction option not supported!"));
	}
}

// ~t_WPLineParams
// ---------------------------------------



