#include "stdafx.h"

#include "MF_CGNS3D.h"
#include "MF_CGNS3D_params.h"

using namespace mf;
using namespace hsstab;

//---------------------------------------------------------------------3D params

const t_CGNS3DParams& t_MFCGNS3D::get_params() const{
	return _base_params;
};

const t_FldParams& t_MFCGNS3D::get_mf_params() const{return _base_params;}

//----------------------------------------------------------------shared init

void mf::cg::hsf3d::_plug_default_settings(TPluginParamsGroup& g){

	// TODO: read all these params from cgns db

	g.add("nu", 5, _("Number of funcs written in cgns fld file"));

	g.add("FldBinPath", _(""), _("Path to binary cgns file containing field")); 

	g.add("FuncNames", _("VelocityX, VelocityY, VelocityZ, Pressure, Temperature"));

	g.add("GrdBinPath", _(""), _("Path to binary cgns file containing grid")); 

	g.add("Mach", 1.0, _("Free stream Mach number")); // pMach

	g.add("Re", 1.0e+06, _T("nondim Reynolds number calculated by a ref length L_ref")); // pRe

	g.add("Alpha", 0.0e+00, _T("Angle of Attack")); // pAlpha

	g.add("Pr", 0.72e+00, _T("Prandtl number")); // pPr

	g.add("Gamma", 1.4e+00, _T("Specific heat ratio")); // pGamma

	g.add("ViscType", 0, _T("Viscosity Law")); // pViscLaw

	g.add("LRef", 1.0e+00, _T("Dimensional reference Length")); // pLRef

	g.add("TInf", 1.0e+02, _T("Dimensional freestream static temperature"));  // pTInf

	g.add("TWall", 1.0e+02, _T("Dimensional wall temperature")); // pTWall

	g.add("TMju", 1.104e+02, _T("Sutherland law reference temperature")); //pT_Mju

	g.add("MjuPow", 0.75e+00, _T("power viscosity law coef")); // pMjuPow

	g.add("MolWeight", 2.7e-02, _T("Dimensional molecular weight of the gas[kg/mol*K]")); 

	g.add("RGas", 8.31e+00, _T("Dimensional universal gas constant [J/mol*K]"));

	g.add("BCWallFamilyNames", _T("Ymin, wall"), _T("BC Family names for viscous wall"));

	g.add("BLCalcType", _T("BY_VELO_DERIV"), _T("Method of computing Boundary Layer thickness"));

	g.add("BLThickTol", 0.1, _T("Parameter-tolerance for a specified BLCalcTypeMethod"));

	g.add("BLThickCoefDefault", 3.0, _T("Thick Coef to use in GetProfiles task"));

	g.add("BBox_Xmin", -1.0, _T("Bounding box Xmin"));
	g.add("BBox_Xmax", 1.0, _T("Bounding box Xmax"));
	g.add("BBox_Ymin", -1.0, _T("Bounding box Ymin"));
	g.add("BBox_Ymax", 1.0, _T("Bounding box Ymax"));
	g.add("BBox_Zmin", -1.0, _T("Bounding box Zmin"));
	g.add("BBox_Zmax", 1.0, _T("Bounding box Zmax"));
}

void mf::cg::hsf3d::_init_fld_base_params(t_FldParams& params, const TPluginParamsGroup& g){


	params.Alpha = g.get_real_param("Alpha");

	params.Gamma = g.get_real_param("Gamma");

	params.L_ref = g.get_real_param("LRef");

	params.Mach = g.get_real_param("Mach");

	params.Mju_pow = g.get_real_param("MjuPow");

	params.Mol_weight = g.get_real_param("MolWeight");

	params.Pr = g.get_real_param("Pr");

	params.R_Gas = g.get_real_param("RGas");

	params.Re = g.get_real_param("Re");

	params.T_inf = g.get_real_param("TInf");

	params.T_mju = g.get_real_param("TMju");

	params.T_wall = g.get_real_param("TWall");

	params.ViscType = g.get_int_param("ViscType");

}

/*
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
		int vvv = mf::t_ViscType::ViscPower;
		ViscType.operator=(vvv);
	}else{
		ViscType = mf::t_ViscType::ViscSuther; 
	}
	conf.SetPath(_T("/"));
}
*/

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