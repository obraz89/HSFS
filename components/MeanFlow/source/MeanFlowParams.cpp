#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>

#include "MeanFlow.h"
#include "ParseHelper.h"
#include <boost/algorithm/string.hpp>
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
