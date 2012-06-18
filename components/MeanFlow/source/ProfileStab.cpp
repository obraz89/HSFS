#include "MeanFlow.h"
#include "ProfileNS.h"
#include "ProfileStab.h"
#include <cmath>


t_ProfileStab::t_ProfileStab(const t_MeanFlow& a_rFld):t_Profile(a_rFld, 0){};
t_ProfileStab::~t_ProfileStab(){};

void t_ProfileStab::initialize(t_ProfileNS& a_rProfNS, int nnodes/* =0*/){
	if (nnodes>0){
		_resize(nnodes);
	}else{
		_resize(a_rProfNS._nnodes);
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

	this->stabRe = sqrt(Params.Re*u_e*rho_e*x/mu_e);
	this->Me = Params.Mach*u_e/sqrt(t_e);
	this->dels = Params.L_ref*y_scale/sqrt(Params.Re);
	this->ue_dim = _rFld.calc_c_dim(a_rProfNS._i, a_rProfNS._bl_bound_ind, a_rProfNS._k)*Me;
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
// instead of derivs for viscosity special coefs mu1 mu2 are stored:
// dmy/dy = mu1*(dt/dy);
// d2my/dy2 = mu2*(dt/dy)^2.0 + mu1*(d2t/dy2)
		_mu[i]=_mu[i]/mu_e;
		if (Params.ViscType==t_MFParams::t_ViscType::ViscPower){
			const double& visc_power = Params.Mju_pow;
			_mu1[i]=visc_power*pow(_t[i],visc_power-1.0);
			_mu2[i]=visc_power*(visc_power-1.0)*pow(_t[i],visc_power-2.0);
		}
		else{
// TODO: check Stagnation ?
		    const double& lt=_t[i];
			const double& t_coef=Params.T_mju*(1.0+t_e);
			_mu1[i]=1.5*_mu[i]/lt-_mu[i]/(lt+t_coef);
			_mu2[i]=1.5*(_mu1[i]-_mu[i]/lt)/lt-
				   (_mu1[i]-_mu[i]/(lt+t_coef))/(lt+t_coef);
		}
	}
};

void t_ProfileStab::initialize(int a_i, int a_k, int nnodes){
	t_ProfileNS ns_prof(_rFld);
	ns_prof.initialize(a_i, a_k);
	this->initialize(ns_prof, nnodes);
};

void t_ProfileStab::initialize(int a_i, int a_k){
	initialize(a_i, a_k, 0);
}

