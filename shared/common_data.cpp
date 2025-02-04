#include "stdafx.h"
#include "common_data.h"

using namespace hsstab;

// task
wxString hsstab::CASE_SETTINGS_DIR = _T("settings");
wxString hsstab::OUTPUT_DIR = _T("output");
wxString hsstab::LOG_FILE = _T("log.txt");

// component names
const wxString cmpnts::MF_HSFLOW3D_NAME = _T("MFHSFLOW3D");
const wxString cmpnts::MF_HSFLOW2D_NAME = _T("MFHSFLOW2D");
const wxString cmpnts::PF_LOCSRCH_NAME = _T("PF.LocSearch");
const wxString cmpnts::PF_GLOBSRCH_NAME = _T("PF.GlobSearch");

const wxString cmpnts::WPTRACK_DEFAULT_NAME = _T("WPTrack");

// config domains
const wxString cmpnts::MF_CONF_DOMAIN = _T("MeanFlow");
const wxString cmpnts::STABSOLVER_CONF_DOMAIN = _T("LocalSearch");
const wxString cmpnts::EIGEN_CONF_DOMAIN= _T("GlobalSearch");

// solving state
TState G_State;