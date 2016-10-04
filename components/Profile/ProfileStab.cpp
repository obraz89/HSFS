#include <iostream>
#include "fstream"

#include <wx/log.h>

#include "ProfileStab.h"

#include <cmath>

#include "io_helpers.h"

using namespace mf;

t_ProfileStab::t_ProfileStab():t_Profile(0){};
t_ProfileStab::~t_ProfileStab(){};

const t_StabScales& t_ProfileStab::scales() const{return _scales;};

/************************************************************************/
// interpolate to uniform grid
/************************************************************************/
void t_ProfileStab::initialize(t_ProfileNS& a_rProfNS, t_ProfStabCfg cfg){

	if (cfg.NNodes>0){
		_resize(cfg.NNodes);
	}else{
		_resize(a_rProfNS.get_nnodes());	
	};

	std::vector<double> y_distrib(_nnodes);

	double dy = a_rProfNS.get_thick()/(_nnodes-1);
	for (int i=0; i<_nnodes; i++) y_distrib[i] = i*dy;

	_initialize(a_rProfNS, y_distrib, cfg);

}

/************************************************************************/
// interpolate to a grid with given distribution of nodes
// y_distrib is scaled as all mf data (just as all values in t_ProfileNS)
/************************************************************************/
void t_ProfileStab::initialize(t_ProfileNS& a_rProfNS, 
			const std::vector<double>& y_distrib ,t_ProfStabCfg cfg){

	if (cfg.NNodes>0){
		_resize(cfg.NNodes);
	}else{
		_resize(a_rProfNS.get_nnodes());	
	};

	_initialize(a_rProfNS, y_distrib, cfg);

}

/************************************************************************/
void t_ProfileStab::_initialize(t_ProfileNS& a_rProfNS, 
			const std::vector<double>& a_y_distrib, t_ProfStabCfg cfg){

	t_Profile::t_Rec ns_outer_rec;

	if (!a_rProfNS.is_disturbance_profile()){

		ns_outer_rec = a_rProfNS.get_bound_rec();

	}else{

		// tmp, to initialize disturbance "stability profiles" on a flat plate
		wxLogMessage(_T("Warning: Initialization for disturbance profiles only!"));
		ns_outer_rec.u = 1.0;
		ns_outer_rec.v = 0.0;
		ns_outer_rec.w = 0.0;

		ns_outer_rec.r = 1.0;
		ns_outer_rec.t = 1.0;
		ns_outer_rec.mu = 1.0;


	}

	double mu_e = ns_outer_rec.mu;

	double u_e = ns_outer_rec.u;

	double rho_e = ns_outer_rec.r;

	double t_e = ns_outer_rec.t;

	const mf::t_DomainBase& rMF = a_rProfNS.getMFDomain();

	const t_FldParams& Params = rMF.get_mf_params();

	double bl_thick_scale = a_rProfNS.get_bl_thick_scale();

	// old selfsim scale, TODO: keep as option to nondim ?
	double x_scale = a_rProfNS.get_x_scale();
	double y_scale_selfsim = sqrt(mu_e*x_scale/(u_e*rho_e));
	double y_scale_bl = bl_thick_scale*sqrt(Params.Re);

	double y_scale;
	// TODO: make an option when adequate x_scale calculations added
	// variant 1
	// use scaling with bl_thick_scale provided by CFD 

	switch (cfg.NondimScaleType)
	{
	case t_ProfStabCfg::NONDIM_BY_CFD_SCALE:
		y_scale = y_scale_bl;
		_scales.ReStab = rho_e*u_e*bl_thick_scale/mu_e*Params.Re;
		_scales.Dels = Params.L_ref*bl_thick_scale;
		break;
	case t_ProfStabCfg::NONDIM_BY_X_SELFSIM:
		y_scale = y_scale_selfsim;
		_scales.ReStab = sqrt(Params.Re*u_e*rho_e*x_scale/mu_e);
		_scales.Dels = Params.L_ref*y_scale/sqrt(Params.Re);
		break;
	default:
		wxLogError(_T("Unsupported option for profile stab non dim!"));
		break;
	}

	_scales.Me = Params.Mach*u_e/sqrt(t_e);

	_scales.UeDim = rMF.calc_c_dim(t_e)*_scales.Me;

	_scales.Ue = u_e;

	// order important - first interpolate then nondim
	for (int i=0; i<_nnodes; i++){

		double cur_y = a_y_distrib[i];

		// interpolate
		set_rec(a_rProfNS.get_rec(cur_y), i);

		//nondim
		_y[i] = _y[i]/y_scale;

		_u[i] =_u[i]/u_e;
		_u1[i]=_u1[i]*y_scale/u_e;
		_u2[i]=_u2[i]*pow(y_scale,2)/u_e;

		_t[i]=_t[i]/t_e;
		_t1[i]=_t1[i]*y_scale/t_e;
		_t2[i]=_t2[i]*pow(y_scale,2)/t_e;

		// IMPORTANT TODO: why isn't w nondim by u_e in orig solver?
		_w[i]=_w[i]/u_e;
		_w1[i]=_w1[i]*y_scale/u_e;
	    _w2[i]=_w2[i]*pow(y_scale,2)/u_e;

		// for viscosity we store dmu/dt and d2mu/dt2
		_mu[i]=_mu[i]/mu_e;

		// _v is not used (assuming 0 in parallel flow assumption)
		// rescale for verification purpose
		_v[i] = _v[i]/u_e;

		if (Params.ViscType==mf::t_ViscType::ViscPower){

			const double& visc_power = Params.Mju_pow;

			_mu1[i]=visc_power*pow(_t[i],visc_power-1.0);
			_mu2[i]=visc_power*(visc_power-1.0)*pow(_t[i],visc_power-2.0);

		}
		else{

		    const double lt=_t[i];

			const double t_suth=Params.T_mju/(Params.T_inf*t_e);
			const double d = 1.5/lt-1.0/(lt+t_suth);

			_mu1[i] = _mu[i]*d;
			_mu2[i] = _mu1[i]*d-_mu[i]/lt*(d+t_suth/pow(lt+t_suth,2));

		}
	}
};

void t_ProfileStab::initialize_2D(const std::string& wfname, const t_StabScales& a_scales)
{

	std::wifstream ifstr(&wfname[0]);
	std::wstringstream istr;

	int n_line=0;
	int nnodes=0;
	const int max_lsize = 1000;
	wxChar line[max_lsize];
	wxChar ch;
	
	_scales = a_scales;

	// read-process profiles size
	ifstr.get(line, max_lsize, '\n');
	ifstr.get(ch);
	istr.clear();
	istr<<line;
	istr>>nnodes;
	_resize(nnodes);

	// read profiles
	while(ifstr.get(line, max_lsize, '\n')){
		if (ifstr.get(ch) && ch!='\n'){
			wxString msg = _T("failed to initialize stability profile from file: line exceeded");
			wxLogMessage(msg);
			ssuGENTHROW(msg);
		}
		istr.clear();
		istr<<line;
		
		io_hlp::write_to_val<double>(istr, _y[n_line]);

		io_hlp::write_to_val<double>(istr, _u[n_line]);
		io_hlp::write_to_val<double>(istr, _u1[n_line]);
		io_hlp::write_to_val<double>(istr, _u2[n_line]);

		io_hlp::write_to_val<double>(istr, _t[n_line]);
		io_hlp::write_to_val<double>(istr, _t1[n_line]);
		io_hlp::write_to_val<double>(istr, _t2[n_line]);

		io_hlp::write_to_val<double>(istr, _mu[n_line]);
		io_hlp::write_to_val<double>(istr, _mu1[n_line]);
		io_hlp::write_to_val<double>(istr, _mu2[n_line]);

		_w[n_line] = _w1[n_line] = _w2[n_line] = 0.0;

		n_line++;
	}
}

void t_ProfileStab::initialize_3D(const std::string& wfname, 
							   const t_StabScales& a_scales)
{

	std::wifstream ifstr(&wfname[0]);
	std::wstringstream istr;

	int n_line=0;
	int nnodes=0;
	const int max_lsize = 1000;
	wxChar line[max_lsize];
	wxChar ch;

	_scales = a_scales;

	// read-process profiles size
	ifstr.get(line, max_lsize, '\n');
	ifstr.get(ch);
	istr.clear();
	istr<<line;
	istr>>nnodes;
	_resize(nnodes);

	// read profiles
	while(ifstr.get(line, max_lsize, '\n')){
		if (ifstr.get(ch) && ch!='\n'){
			wxString msg = _T("failed to initialize stability profile from file: line exceeded");
			wxLogMessage(msg);
			ssuGENTHROW(msg);
		}
		istr.clear();
		istr<<line;

		io_hlp::write_to_val<double>(istr, _y[n_line]);

		io_hlp::write_to_val<double>(istr, _u[n_line]);
		io_hlp::write_to_val<double>(istr, _u1[n_line]);
		io_hlp::write_to_val<double>(istr, _u2[n_line]);

		io_hlp::write_to_val<double>(istr, _t[n_line]);
		io_hlp::write_to_val<double>(istr, _t1[n_line]);
		io_hlp::write_to_val<double>(istr, _t2[n_line]);

		io_hlp::write_to_val<double>(istr, _w[n_line]);
		io_hlp::write_to_val<double>(istr, _w1[n_line]);
		io_hlp::write_to_val<double>(istr, _w2[n_line]);

		io_hlp::write_to_val<double>(istr, _mu[n_line]);
		io_hlp::write_to_val<double>(istr, _mu1[n_line]);
		io_hlp::write_to_val<double>(istr, _mu2[n_line]);

		n_line++;
	}
}

