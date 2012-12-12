#include "MeanFlow.h"
#include "ProfileNS.h"
#include "ProfileStab.h"

#include <cmath>
#include "io_helpers.h"

#include "log.h"


t_ProfileStab::t_ProfileStab(const t_MeanFlow& a_rFld):t_Profile(a_rFld, 0){};
t_ProfileStab::~t_ProfileStab(){};

void t_ProfileStab::initialize(t_ProfileNS& a_rProfNS, int nnodes/* =0*/){
	if (nnodes>0){
		_resize(nnodes);
	}else{
		_resize(a_rProfNS._nnodes);	// TODO: fix : should be stab_params nnodes default
	};
	// interpolate to uniform grid and
	// non-dimensionalize y and all derivs
	//to A = sqrt(nu_e*x/u_e)
	// all values in eq. for A are dimensional
	// NS profiles are nondim by sqrt(Re)
	const t_Profile::t_Rec& bl_outer_rec = a_rProfNS.get_bl_bound_rec();
	double mu_e = bl_outer_rec.mu;
	double x = a_rProfNS.xDist;
	double u_e = bl_outer_rec.u;
	double rho_e = bl_outer_rec.r;
	double t_e = bl_outer_rec.t;
	double y_scale = sqrt(mu_e*x/(u_e*rho_e));
	const t_MFParams& Params = _rFld.base_params();
	// keep Jacobian to local RF
	// in order to finally pass it to wavepack line
	_jacToLocalRF = a_rProfNS._jacToLocalRF;

	_scales.ReStab = sqrt(Params.Re*u_e*rho_e*x/mu_e);
	_scales.Me = Params.Mach*u_e/sqrt(t_e);
	_scales.Dels = Params.L_ref*y_scale/sqrt(Params.Re);
	_scales.UeDim = 
		_rFld.calc_c_dim(a_rProfNS._i, a_rProfNS._bl_bound_ind, a_rProfNS._k)*_scales.Me;
	_scales.Ue = u_e;
	double dy = (a_rProfNS._y[a_rProfNS.size()-1])/((double)this->size());
	for (int i=0; i<size(); i++){
		// order important - first interpolate then nondim
		double cur_y = (double)i*dy;
		int sizeNS = a_rProfNS.size();
		set_rec(a_rProfNS.get_rec(cur_y), i);
		/*
		_y[i] = _interpolate(cur_y, a_rProfNS._y, a_rProfNS._y, sizeNS);

		_u[i] = _interpolate(cur_y, a_rProfNS._y, a_rProfNS._u, sizeNS);
		_u1[i] = _interpolate(cur_y, a_rProfNS._y, a_rProfNS._u1, sizeNS);
		_u2[i] = _interpolate(cur_y, a_rProfNS._y, a_rProfNS._u2, sizeNS);

		_t[i] = _interpolate(cur_y, a_rProfNS._y, a_rProfNS._t, sizeNS);
		_t1[i] = _interpolate(cur_y, a_rProfNS._y, a_rProfNS._t1, sizeNS);
		_t2[i] = _interpolate(cur_y, a_rProfNS._y, a_rProfNS._t2, sizeNS);

		_w[i] = _interpolate(cur_y, a_rProfNS._y, a_rProfNS._w, sizeNS);
		_w1[i] = _interpolate(cur_y, a_rProfNS._y, a_rProfNS._w1, sizeNS);
		_w2[i] = _interpolate(cur_y, a_rProfNS._y, a_rProfNS._w2, sizeNS);

		_mu[i] = _interpolate(cur_y, a_rProfNS._y, a_rProfNS._mu, sizeNS);
		_mu1[i] = _interpolate(cur_y, a_rProfNS._y, a_rProfNS._mu1, sizeNS);
		_mu2[i] = _interpolate(cur_y, a_rProfNS._y, a_rProfNS._mu2, sizeNS);

		_p[i] = _interpolate(cur_y, a_rProfNS._y, a_rProfNS._p, sizeNS);
		_r[i] = _interpolate(cur_y, a_rProfNS._y, a_rProfNS._r, sizeNS);
		*/
		//nondim
		_y[i] = _y[i]/y_scale;

		_u[i] =_u[i]/u_e;
		_u1[i]=_u1[i]*y_scale/u_e;
		_u2[i]=_u2[i]*pow(y_scale,2)/u_e;

		_t[i]=_t[i]/t_e;
		_t1[i]=_t1[i]*y_scale/t_e;
		_t2[i]=_t2[i]*pow(y_scale,2)/t_e;

// TODO: why isn't w nondim by u_e in orig solver?
		_w[i]=_w[i]/u_e;
		_w1[i]=_w1[i]*y_scale/u_e;
	    _w2[i]=_w2[i]*pow(y_scale,2)/u_e;
// for viscosity we store dmu/dt and d2mu/dt2
		_mu[i]=_mu[i]/mu_e;
		if (Params.ViscType==t_MFParams::t_ViscType::ViscPower){
			const double& visc_power = Params.Mju_pow;
			_mu1[i]=visc_power*pow(_t[i],visc_power-1.0);
			_mu2[i]=visc_power*(visc_power-1.0)*pow(_t[i],visc_power-2.0);
		}
		else{
		    const double& lt=_t[i];
			const double& t_suth=Params.T_mju/Params.T_inf;
			const double d = 1.5/lt-1.0/(lt+t_suth);
			_mu1[i] = _mu[i]*d;
			_mu2[i] = _mu1[i]*d-_mu[i]/lt*(d+t_suth/pow(lt+t_suth,2));
		}
	}
};

void t_ProfileStab::initialize(int a_i, int a_k, double a_thick_coef,int nnodes){
	t_ProfileNS ns_prof(_rFld);
	ns_prof.initialize(a_i, a_k, a_thick_coef);
	this->initialize(ns_prof, nnodes);
};

void t_ProfileStab::initialize(int a_i, int a_k, double a_thick_coef){
	initialize(a_i, a_k, a_thick_coef, 0);
}

void t_ProfileStab::initialize(const std::wstring wfname){
	std::string fname = wx_to_stdstr(wxString(&wfname[0]));
	std::ifstream ifstr(&fname[0]);
	std::stringstream istr;

	int n_line=0;
	int nnodes=0;
	const int max_lsize = 1000;
	char line[max_lsize];
	char ch;

	// read R
	ifstr.get(line, max_lsize, '\n');
	ifstr.get(ch);
	istr.clear();
	istr<<line;
	istr>>_scales.ReStab;

	// read Me
	ifstr.get(line, max_lsize, '\n');
	ifstr.get(ch);
	istr.clear();
	istr<<line;
	istr>>_scales.Me;

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
			wxString msg = _("failed to initialize stability profile from file: line exceeded");
			Log<<msg;
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

