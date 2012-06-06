#include <cmath>
#include <sstream>
#include <fstream>

#include "MFParams.h"
#include "ParseHelper.h"
#include <boost/algorithm/string.hpp>

#include "wx/fileconf.h"
#include "wx/string.h"
#include "common_data.h"

const int t_MFParams::t_ViscType::ViscPower=0;
const int t_MFParams::t_ViscType::ViscSuther=1;

const int t_MFParamsHS2D::t_AxeSym::AxeSym = 0;
const int t_MFParamsHS2D::t_AxeSym::Plane  = 1;

t_MFParams::t_MFParams():t_ComponentParamsGroup(MF_CONF_DOMAIN){};

t_MFParams::t_MFParams(wxString configfile):t_ComponentParamsGroup(MF_CONF_DOMAIN){
	_init_params_map();
	wxFileConfig conf = _get_config_handle(configfile);
	_load_via_params(conf);
}

void t_MFParams::_init_params_map(){
	t_CompParamStr* pMfBin = 
		new t_CompParamStr(mf_bin_path, _T("MF_bin_path"), _T("Absolute path to a Mean Flow Binary"));
	pMfBin->set_default(_T("")); // hmmm 
	_add_param(pMfBin);

	t_CompParamInt* pNx = 
		new t_CompParamInt(Nx, _T("Nx"), _T("Nx"));
	pNx->set_default(1);
	_add_param(pNx);

	t_CompParamInt* pNy = 
		new t_CompParamInt(Ny, _T("Ny"), _T("Ny"));
	pNy->set_default(1);
	_add_param(pNy);

	t_CompParamInt* pNz = 
		new t_CompParamInt(Nz, _T("Nz"), _T("Nz"));
	pNz->set_default(1);
	_add_param(pNz);

	t_CompParamDbl* pMach = 
		new t_CompParamDbl(Mach, _T("Mach"), _T("Free stream Mach number"));
	pMach->set_default(1.0e+00);
	_add_param(pMach);
	
	t_CompParamDbl* pRe = 
		new t_CompParamDbl(Re, _T("Re"), _T("nondim Reynolds number calculated by a ref length L_ref"));
	pRe->set_default(1.0e+06);
	_add_param(pRe);

	t_CompParamDbl* pAlpha = 
		new t_CompParamDbl(Alpha, _T("Alpha"), _T("Angle of Attack"));
	pAlpha->set_default(0.0e+00);
	_add_param(pAlpha);

	t_CompParamDbl* pPr = 
		new t_CompParamDbl(Pr, _T("Pr"), _T("Prandtl number"));
	pPr->set_default(0.72);
	_add_param(pPr);

	t_CompParamDbl* pGamma = 
		new t_CompParamDbl(Gamma, _T("Gamma"), _T("Specific heat ratio"));
	pGamma->set_default(1.4);
	_add_param(pGamma);

	t_CompParamInt* pViscLaw = 
		new t_CompParamInt(ViscType, _T("ViscType"), _T("Viscosity Law"));
	pViscLaw->set_default(0l);
	_add_param(pViscLaw);
	
	t_CompParamDbl* pLRef = 
		new t_CompParamDbl(L_ref, _T("L_ref"), _T("Dimensional Reference Length"));
	pLRef->set_default(1.0e+00);
	_add_param(pLRef);

	t_CompParamDbl* pTInf = 
		new t_CompParamDbl(T_inf, _T("T_inf"), _T("Freestream static temperature"));
	pTInf->set_default(100.0e+00);
	_add_param(pTInf);

	t_CompParamDbl* pTWall = 
		new t_CompParamDbl(T_wall, _T("T_wall"), _T("Dimensional wall temperature"));
	pTWall->set_default(100.0e+00);
	_add_param(pTWall);

	t_CompParamDbl* pTMju = 
		new t_CompParamDbl(T_mju, _T("T_mju"), _T("Sutherland law reference temperature"));
	pTMju->set_default(110.4e+00);
	_add_param(pTMju);

	t_CompParamDbl* pMjuPow = 
		new t_CompParamDbl(Mju_pow, _T("Mju_pow"), _T("power viscosity law coef"));
	pMjuPow->set_default(0.75);
	_add_param(pMjuPow);
};

void t_MFParams::_load_direct(wxFileConfig& conf){
	conf.SetRecordDefaults(); //write defaults to config file
	conf.SetPath(ConfigDomain);
	this->mf_bin_path = conf.Read(_T("MF binary path"), _T(""));
	this->Nx = conf.Read(_T("Nx"), 1);
	this->Ny = conf.Read(_T("Ny"), 1);
	this->Nz = conf.Read(_T("Nz"), 1);
	this->Mach = conf.Read(_T("Mach"), 1.0);
	this->Re = conf.Read(_T("Re"), 1.0e+06);
	this->Alpha = conf.Read(_T("Alpha"), 0.0);
	this->Pr = conf.Read(_T("Pr"), 0.72);
	this->Gamma = conf.Read(_T("Gamma"), 1.4);
	this->T_inf = conf.Read(_T("T_inf"), 90.318);
	this->Mju_pow = conf.Read(_T("Mju_pow"), 0.75);
	this->L_ref = conf.Read(_T("L_ref"), 0.381);
	if (conf.Read(_T("ViscType"), 0l)){
		int vvv = t_ViscType::ViscPower;
		ViscType.operator=(vvv);
	}else{
		ViscType = t_ViscType::ViscSuther; 
	}
	conf.SetPath(_T("/"));
}

void t_MFParams::load_direct(wxString configfile){
	_load_direct(_get_config_handle(configfile));
};

void t_MFParams::_load_via_params(wxFileConfig& handle){
	t_ComponentParamsGroup::_load_via_params(handle);
};

void t_MFParams::load_via_params(wxString configfile){
	_load_via_params(_get_config_handle(configfile));
};

void t_MFParams::save(wxString configfile){
	// for future
};


/*
std::string t_MeanFlow::t_Params::_get_conf_dir(std::string conf_path){
	int found=conf_path.find_last_of("/\\");
	return conf_path.substr(0, found);
}
void t_MeanFlow::t_Params::_init(const std::string a_conf_file_path){
	// start parsing config file
	std::ifstream f_str(&a_conf_file_path[0]);
	std::stringstream s_buf;
	int line_number=0;
	std::string cur_line;

	while (std::getline(f_str, cur_line)){
		s_buf<<cur_line;
		switch (line_number){
			// define local scope
			case 0:{
				// ttl string - omit ?
				std::string ttl_str = parse::get_val<std::string>(cur_line);
				std::vector<std::string> split_cont;
				boost::split(split_cont, ttl_str, boost::is_any_of("."),boost::token_compress_on);
				split_cont.back()="dat";
				_mf_bin_path = _get_conf_dir(a_conf_file_path).append("/").append(boost::join(split_cont, "."));
				// TODO:
				// if (!SearchPath(...)) throw BadMFFile;
				std::cout<<_mf_bin_path<<std::endl;
				break;
			}
			case 1:
			// grid file - omit ?
				break;
			case 2:
			// nx
				Nx = parse::get_val<int>(cur_line);
				break;
			case 3:
			// ny
				Ny = parse::get_val<int>(cur_line);
				break;
			case 4:
			// nz
				Nz = parse::get_val<int>(cur_line);
				break;
			case 5:
			// nread
				break;
			case 6:
			// wall parameter
				break;
			case 7:
			// number of time steps
				break;
			case 8:
			// tau - temp step
				break;
			case 9:
			// internal solver tolerance
				break;
			case 10:
			// Mach
				Mach = parse::get_val<double>(cur_line);
				break;
			case 11:
			// Re
				Re = parse::get_val<double>(cur_line);
				break;
			case 12:
			// angle of attack
				break;
			case 13:
			// non-dim wall temperature (if non-adiabatic wall-condition used)
				break;
			case 14:
			// max non-linear iters
				break;
			case 15:
			// linear solver regime parameter ???
				break;
			case 16:
			// max iters without Jac recalcs
				break;
			case 17:
			// Newton iteration start parameter
				break;
			case 18:
			// minmod rounding 
				break;
			case 19:
			// grid-type 
				break;
			case 20:
			// time-order
				break;
			case 21:
			// 2-time-order parameter
				break;
			case 22:
			// write out data every -- iterations
				break;
			case 23:
			// number of GMRES iters
				break;
			// for a while, later decide what to do with .ttl
			default:
				break;
		}
		line_number++;
	}
	// check necessary constants initialization
	// TODO: must be in config NS file
	// for now it is necessary to change by hand
	Pr = 0.72;
	Gamma = 1.4;
	T_inf = 90.318;
	T_mju = 110.4/T_inf;
	Mju_pow = 0.75;
	L_ref = 0.381;
	Alpha = 2.0;
	ViscType = ViscPower; // power
};
// initializing from 3D field
t_MeanFlow::t_Params::t_Params(const std::string conf_path){
	_init(conf_path);
}
// initializing from 2D field
t_MeanFlow::t_Params::t_Params(const std::string conf_path, int kk){
	_init(conf_path);
	Nz = kk;
};

t_MeanFlow::t_Params::~t_Params(){};
*/

t_MFParamsHS2D::t_MFParamsHS2D(wxString configfile):t_MFParams(){
	_init_params_map();
	load_via_params(configfile);
};

void t_MFParamsHS2D::_init_params_map(){
	t_MFParams::_init_params_map();
	
	t_CompParamInt* pAxeSym = 
		new t_CompParamInt(MFSym, _T("Gamma"), _T("Specific heat ratio"));
	pAxeSym->set_default(0l);
	_add_param(pAxeSym);
};

void t_MFParamsHS2D::load_via_params(wxString configfile){
	t_MFParams::load_via_params(configfile);
};

void t_MFParamsHS2D::load_direct(wxString configfile){
	t_MFParams::load_direct(configfile);
};

void t_MFParamsHS3D::load_via_params(wxString configfile){
	// no additional params
	t_MFParams::load_via_params(configfile);
};

void t_MFParamsHS3D::load_direct(wxString configfile){
	// no additional params
	t_MFParams::load_direct(configfile);
};