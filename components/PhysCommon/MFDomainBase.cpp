#include "stdafx.h"
#include "MFDomainBase.h"

mf::t_DomainBase::t_DomainBase(){};

mf::t_DomainBase::~t_DomainBase(){};

void mf::t_DomainBase::set_bl_thick_calc_type(t_BLThickCalcType v){_bl_thick_ctype = v;}

double mf::t_DomainBase::calc_enthalpy(const mf::t_Rec& ptr) const{

	const t_FldParams& mf_prms = get_mf_params();

	double cp_t, v_2;

	cp_t = ptr.t/((mf_prms.Gamma-1.0)*mf_prms.Mach*mf_prms.Mach);

	v_2 = 0.5*(pow(ptr.u,2.0)+pow(ptr.v,2.0)+pow(ptr.w,2.0));

	return (cp_t + v_2);

}

double mf::t_DomainBase::calc_enthalpy(const t_GeomPoint& xyz) const
{

	const t_Rec& ptr = get_rec(xyz);

	return calc_enthalpy(ptr);

};

double mf::t_DomainBase::calc_enthalpy_freestream() const{

	const t_FldParams& mf_prms = get_mf_params();

	// TODO: is this always correct?
	// is always |u|=1
	return (1.0/((mf_prms.Gamma-1.0)*mf_prms.Mach*mf_prms.Mach) + 0.5);

};

double mf::t_DomainBase::calc_viscosity(const double t) const{

	const t_FldParams& mf_prms = get_mf_params();

	if (mf_prms.ViscType==t_ViscType::ViscPower){

		return pow(t, mf_prms.Mju_pow);

	}
	else{

		double t_suth = mf_prms.T_mju/mf_prms.T_inf;

		return pow(t, 1.5)*(1.0+t_suth)/(t+t_suth);

	}

}

// calc nondim viscosity
double mf::t_DomainBase::calc_viscosity(const t_GeomPoint& xyz) const{

	const t_Rec& rRec = get_rec(xyz);

	return calc_viscosity(rRec.t);
};


double mf::t_DomainBase::calc_cin_visc_inf_dim() const{

	const t_FldParams& mf_prms = get_mf_params();

	double u_inf_dim = calc_u_inf();

	return u_inf_dim*mf_prms.L_ref/mf_prms.Re;

};

double mf::t_DomainBase::calc_cin_visc_dim(const t_GeomPoint& xyz) const{

	const t_FldParams& mf_prms = get_mf_params();

	double nju_inf = calc_cin_visc_inf_dim();

	const t_Rec& rec = get_rec(xyz);
	return nju_inf*calc_viscosity(rec.t)/rec.r;

};

double mf::t_DomainBase::calc_c_dim(const double t) const{

	const t_FldParams& mf_prms = get_mf_params();

	double t_dim = t*mf_prms.T_inf;

	return sqrt(mf_prms.Gamma*mf_prms.R_Gas*t_dim/mf_prms.Mol_weight);

}

double mf::t_DomainBase::calc_c_dim(const t_GeomPoint& xyz) const{

	const t_FldParams& mf_prms = get_mf_params();

	const t_Rec& rRec = get_rec(xyz);
	double t = rRec.t;
	double t_dim = t*mf_prms.T_inf;

	return sqrt(mf_prms.Gamma*mf_prms.R_Gas*t_dim/mf_prms.Mol_weight);

}

double mf::t_DomainBase::calc_u_inf() const{

	const t_FldParams& mf_prms = get_mf_params();

	return mf_prms.Mach*
		sqrt(mf_prms.Gamma*mf_prms.R_Gas*mf_prms.T_inf/mf_prms.Mol_weight);
};

double mf::t_DomainBase::calc_mach(const t_GeomPoint& xyz) const{

	const t_Rec& rRec = get_rec(xyz);

	const t_FldParams& mf_prms = get_mf_params();

	double vAbs = sqrt(pow(rRec.u,2.0)+pow(rRec.v,2.0)+pow(rRec.w,2.0));

	return mf_prms.Mach*vAbs/sqrt(rRec.t);

}


//===============================================================