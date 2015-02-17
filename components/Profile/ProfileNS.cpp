#include "ProfileNS.h"
#include "math_operands.h"
#include "smooth.h"

using namespace mf;

t_ProfileNS::t_ProfileNS(const t_DomainBase& a_rDomain):t_ProfMF(a_rDomain){};

void t_ProfileNS::_initialize_interpolate(const t_GeomPoint& xyz, const mf::t_ProfDataCfg& init_cfg){

	//IMPORTANT TODO: remove all "estimations" when 2-nd order interpolation done!!!
	const double& a_thick_coef = init_cfg.ThickCoef;
	const int& a_nnodes = init_cfg.NNodes;

	int num_bl_nodes_ns = _rDomain.estim_num_bl_nodes(xyz);

	if (num_bl_nodes_ns<=0 ){
		wxLogError(_T("ProfileNS: failed to estimate bl nodes number"));
		if (a_nnodes<=0) ssuGENTHROW(_T("ProfileNS: failed to calculate nnodes"));
	}


	wxLogMessage(_T("ProfileNS: Check _initialize_interpolate"));
	a_nnodes>0 ? _resize(a_nnodes) : _resize(num_bl_nodes_ns);

	t_GeomPoint r_xyz_base;
	t_Vec3Dbl surf_norm;
	_rDomain.calc_surf_point(xyz, r_xyz_base, surf_norm);

	_xyz = r_xyz_base;
	double bl_thick = _rDomain.calc_bl_thick(_xyz);

	double cur_eta = bl_thick;

	double prof_thick = a_thick_coef*bl_thick;

	// transform vector fields and coordinates
	// to a new reference frame
	t_Vec3Dbl u_xyz, u_ked, dr_xyz, r_ked;
	t_GeomPoint r_xyz;

	t_SqMat3Dbl jac = _rDomain.calc_jac_to_loc_rf(r_xyz_base);
	t_SqMat3Dbl inv_jac = jac.inverse();

	double deta = prof_thick/double(_nnodes-1);

	for (int j=0; j<_nnodes; j++){

		dr_xyz = j*deta*surf_norm;
		r_xyz = r_xyz_base + dr_xyz;
		const mf::t_Rec mf_rec = _rDomain.interpolate_to_point(r_xyz);
		
		r_ked = inv_jac*dr_xyz;

		u_xyz = mf_rec.get_uvw();
		u_ked = inv_jac*u_xyz;
		
		// convention : to make y = O(1)
		// it mul by factor sqrt(ReINF)
		// as it is for local parallel task we store only y assuming x=z=0
		// TODO: (good debug check is to assure that)
		const mf::t_FldParams& Params = _rDomain.get_mf_params();
		_y[j] = r_ked[1]*sqrt(Params.Re);
		// again we postulate v==0
		_u[j] = u_ked[0];
		_w[j] = u_ked[2];

		_p[j] = mf_rec.p;
		_t[j] = mf_rec.t;
		// rho is non-dim as follows
		double gMaMa = Params.Gamma*pow(Params.Mach, 2);
		_r[j] = mf_rec.p/mf_rec.t*gMaMa;
		_mu[j]=_rDomain.calc_viscosity(_t[j]);
	}

	// SMOOTHING
	// FORTRAN CALLING CONVENTION

	int nnodes = get_nnodes();
	SMOOTH_3D_PROFILES(&_y[0], &_u[0], &nnodes, &_u1[0], &_u2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_w[0], &nnodes, &_w1[0], &_w2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_t[0], &nnodes, &_t1[0], &_t2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_mu[0], &nnodes, &_mu1[0], &_mu2[0]);
}


void t_ProfileNS::_initialize_extract(const t_GeomPoint& xyz, const mf::t_ProfDataCfg& init_cfg){

	std::vector<mf::t_Rec> raw_profile;

	_rDomain.extract_profile_data(xyz, init_cfg, raw_profile);

	_resize(raw_profile.size());

	t_GeomPoint r_xyz_base, r_xyz;
	t_Vec3Dbl dr_xyz, r_ked;
	t_Vec3Dbl u_xyz, u_ked;

	r_xyz_base.set(raw_profile[0]);
	_xyz = r_xyz_base;

	t_SqMat3Dbl jac = _rDomain.calc_jac_to_loc_rf(r_xyz_base);
	t_SqMat3Dbl inv_jac = jac.inverse();

	for (int j=0; j<_nnodes; j++){

		const mf::t_Rec& mf_rec = raw_profile[j];

		r_xyz.set(mf_rec);
		matrix::base::minus<double, double>(r_xyz, r_xyz_base, dr_xyz);
		matrix::base::mat_mul<double, double>(inv_jac, dr_xyz, r_ked);

		u_xyz.set(mf_rec.u, mf_rec.v, mf_rec.w);
		matrix::base::mat_mul<double, double>(inv_jac, u_xyz, u_ked);

		// convention : to make y = O(1)
		// it mul by factor sqrt(ReINF)
		// as it is for local parallel task we store only y assuming x=z=0
		// TODO: (good debug check is to assure that)
		const mf::t_FldParams& Params = _rDomain.get_mf_params();
		_y[j] = r_ked[1]*sqrt(Params.Re);
		// again we postulate v==0
		_u[j] = u_ked[0];
		_w[j] = u_ked[2];

		_p[j] = mf_rec.p;
		_t[j] = mf_rec.t;
		// rho is non-dim as follows
		double gMaMa = Params.Gamma*pow(Params.Mach, 2);
		_r[j] = mf_rec.p/mf_rec.t*gMaMa;
		_mu[j]=_rDomain.calc_viscosity(_t[j]);
	}

	// SMOOTHING
	// FORTRAN CALLING CONVENTION

	int nnodes = get_nnodes();
	SMOOTH_3D_PROFILES(&_y[0], &_u[0], &nnodes, &_u1[0], &_u2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_w[0], &nnodes, &_w1[0], &_w2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_t[0], &nnodes, &_t1[0], &_t2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_mu[0], &nnodes, &_mu1[0], &_mu2[0]);
}

//--------------------------------------------------------------------------------~t_ProfileNS

//--------------------------------------------------------------------------------t_ProfileNSGlob

t_ProfMFGlob::t_ProfMFGlob(const t_DomainBase& a_rDomain):t_ProfMF(a_rDomain){};

void t_ProfMFGlob::_initialize_interpolate(const t_GeomPoint& xyz, const mf::t_ProfDataCfg& init_cfg){

	//IMPORTANT TODO: remove all "estimations" when 2-nd order interpolation done!!!
	const double& a_thick_coef = init_cfg.ThickCoef;
	const int& a_nnodes = init_cfg.NNodes;

	int num_bl_nodes_ns = _rDomain.estim_num_bl_nodes(xyz);

	if (num_bl_nodes_ns<=0 ){
		wxLogError(_T("ProfileNS: failed to estimate bl nodes number"));
		if (a_nnodes<=0) ssuGENTHROW(_T("ProfileNS: failed to calculate nnodes"));
	}


	wxLogMessage(_T("ProfileNS: Check _initialize_interpolate"));
	a_nnodes>0 ? _resize(a_nnodes) : _resize(num_bl_nodes_ns);

	t_GeomPoint r_xyz_base;
	t_Vec3Dbl surf_norm;
	_rDomain.calc_surf_point(xyz, r_xyz_base, surf_norm);

	_xyz = r_xyz_base;
	double bl_thick = _rDomain.calc_bl_thick(_xyz);

	double cur_eta = bl_thick;

	double prof_thick = a_thick_coef*bl_thick;

	t_Vec3Dbl u_xyz, dr_xyz;
	t_GeomPoint r_xyz;

	t_SqMat3Dbl jac = _rDomain.calc_jac_to_loc_rf(r_xyz_base);
	t_SqMat3Dbl inv_jac = jac.inverse();

	double deta = prof_thick/double(_nnodes-1);

	for (int j=0; j<_nnodes; j++){

		dr_xyz = j*deta*surf_norm;
		r_xyz = r_xyz_base + dr_xyz;
		const mf::t_Rec mf_rec = _rDomain.interpolate_to_point(r_xyz);

		u_xyz = mf_rec.get_uvw();

		const mf::t_FldParams& Params = _rDomain.get_mf_params();

		_y[j] = dr_xyz.norm();

		_u[j] = u_xyz[0];
		_v[j] = u_xyz[1];
		_w[j] = u_xyz[2];

		_p[j] = mf_rec.p;
		_t[j] = mf_rec.t;
		// rho is non-dim as follows
		double gMaMa = Params.Gamma*pow(Params.Mach, 2);
		_r[j] = mf_rec.p/mf_rec.t*gMaMa;
		_mu[j]=_rDomain.calc_viscosity(_t[j]);
	}

	// SMOOTHING
	// FORTRAN CALLING CONVENTION

	int nnodes = get_nnodes();
	SMOOTH_3D_PROFILES(&_y[0], &_u[0], &nnodes, &_u1[0], &_u2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_w[0], &nnodes, &_w1[0], &_w2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_t[0], &nnodes, &_t1[0], &_t2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_mu[0], &nnodes, &_mu1[0], &_mu2[0]);
}


void t_ProfMFGlob::_initialize_extract(const t_GeomPoint& xyz, const mf::t_ProfDataCfg& init_cfg){

	std::vector<mf::t_Rec> raw_profile;

	_rDomain.extract_profile_data(xyz, init_cfg, raw_profile);

	_resize(raw_profile.size());

	t_GeomPoint r_xyz_base, r_xyz;
	t_Vec3Dbl dr_xyz;

	r_xyz_base.set(raw_profile[0]);
	_xyz = r_xyz_base;

	for (int j=0; j<_nnodes; j++){

		const mf::t_Rec& mf_rec = raw_profile[j];

		r_xyz.set(mf_rec);
		matrix::base::minus<double, double>(r_xyz, r_xyz_base, dr_xyz);

		const mf::t_FldParams& Params = _rDomain.get_mf_params();

		_y[j] = dr_xyz.norm();

		_u[j] = mf_rec.u;
		_v[j] = mf_rec.v;
		_w[j] = mf_rec.w;

		_p[j] = mf_rec.p;
		_t[j] = mf_rec.t;
		// rho is non-dim as follows
		double gMaMa = Params.Gamma*pow(Params.Mach, 2);
		_r[j] = mf_rec.p/mf_rec.t*gMaMa;
		_mu[j]=_rDomain.calc_viscosity(_t[j]);
	}

	// SMOOTHING
	// FORTRAN CALLING CONVENTION

	int nnodes = get_nnodes();
	SMOOTH_3D_PROFILES(&_y[0], &_u[0], &nnodes, &_u1[0], &_u2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_w[0], &nnodes, &_w1[0], &_w2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_t[0], &nnodes, &_t1[0], &_t2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_mu[0], &nnodes, &_mu1[0], &_mu2[0]);
}