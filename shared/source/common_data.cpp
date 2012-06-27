#include "common_data.h"

using namespace common;

// component names
const wxString cmpnts::MF_HSFLOW3D_NAME = _T("MF_HSFLOW3D");
const wxString cmpnts::MF_HSFLOW2D_NAME = _T("MF_HSFLOW2D");
const wxString cmpnts::STABSOLVER3D_NAME= _T("StabSolver");
const wxString cmpnts::EIGEN3D_NAME = _T("EigenGS");

// config domains
const wxString cmpnts::MF_CONF_DOMAIN = _T("MeanFlow");
const wxString cmpnts::STABSOLVER_CONF_DOMAIN = _T("StabSolver");
const wxString cmpnts::EIGEN_CONF_DOMAIN= _T("EigenGS");

wxString io::LOG_FILE = _T("log.txt");
wxString io::LOG_DIR = _T("");