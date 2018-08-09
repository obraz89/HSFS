#include "ProfileNS.h"

// tmp
#include <fstream>

#include <sstream>

#include "wx/log.h"

using namespace mf;

t_ProfileNS::t_ProfileNS(const t_DomainBase& a_rDomain):t_ProfMF(a_rDomain){};

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

	const mf::t_FldParams& Params = _rDomain.get_mf_params();

	for (int j = 0; j<_nnodes; j++) {

		const mf::t_Rec& mf_rec = raw_profile[j];

		if (j == 0) {
			_y[j] = 0.0;
		}
		else {
			r_xyz.set(raw_profile[j]);
			r_xyz_base.set(raw_profile[j - 1]);
			matrix::base::minus<double, double>(r_xyz, r_xyz_base, dr_xyz);
			// convention : to make y = O(1)
			// it mul by factor sqrt(ReINF)
			// as it is for local parallel task we store only y assuming x=z=0
			// TODO: check r_ked[0] and r_ked[2] are small
			_y[j] = _y[j - 1] + dr_xyz.norm()*sqrt(Params.Re);
		}

		//r_xyz.set(mf_rec);
		//matrix::base::minus<double, double>(r_xyz, r_xyz_base, dr_xyz);
		//matrix::base::mat_mul<double, double>(inv_jac, dr_xyz, r_ked);

		u_xyz.set(mf_rec.u, mf_rec.v, mf_rec.w);
		matrix::base::mat_mul<double, double>(inv_jac, u_xyz, u_ked);

		// TODO: check that v is small
		_u[j] = u_ked[0];
		_v[j] = u_ked[1];
		_w[j] = u_ked[2];

		_p[j] = mf_rec.p;
		_t[j] = mf_rec.t;
		// rho is non-dim as follows
		double gMaMa = Params.Gamma*pow(Params.Mach, 2);
		_r[j] = mf_rec.p / mf_rec.t*gMaMa;
		_mu[j] = _rDomain.calc_viscosity(_t[j]);
	}

	_calc_derivs();

}

void _read_crop_raw_profile_from_file(const std::string fname, const mf::t_ProfDataCfg& init_cfg, 
	std::vector<mf::t_Rec>& raw_profile, t_Vec3Dbl& Ue_dir) {

	std::ifstream ifstr(fname);
	std::stringstream istr;

	int nnodes_ns = 0;
	double bl_thick_scale;
	const int max_lsize = 1000;
	char line[max_lsize];
	char ch;

	// read-process profiles size
	ifstr.get(line, max_lsize, '\n');
	ifstr.get(ch);
	istr.clear();
	istr << line;
	istr >> nnodes_ns;
	istr >> bl_thick_scale;
	wxLogMessage(
		_T("Reading profile from file: %d records to read, bl_thick_scale=%lf"),
		nnodes_ns, bl_thick_scale);

	// second line - direction of mean flow velocity at bl edge
	// read-process profiles size
	ifstr.get(line, max_lsize, '\n');
	ifstr.get(ch);
	istr.clear();
	istr << line;
	istr >> Ue_dir[0];
	istr >> Ue_dir[1];
	istr >> Ue_dir[2];

	wxLogMessage(_T("Ue_dir:%s"), Ue_dir.to_wxstr());

	// y
	std::vector<t_Vec3Dbl> xyz_vec(nnodes_ns);

	// u, v, w, p, T
	std::vector<double> zero_vec(5);
	std::vector<std::vector<double>> vec_data(nnodes_ns, zero_vec);

	// read profiles
	for (int i = 0; i < nnodes_ns; i++) {

		ifstr.get(line, max_lsize, '\n');
		if (ifstr.get(ch) && ch != '\n') {
			wxString msg = _T("failed to initialize stability profile from file: line exceeded");
			wxLogMessage(msg);
			ssuGENTHROW(msg);
		}
		istr.clear();
		istr << line;

		io_hlp::write_to_val<double>(istr, xyz_vec[i][0]);
		io_hlp::write_to_val<double>(istr, xyz_vec[i][1]);
		io_hlp::write_to_val<double>(istr, xyz_vec[i][2]);

		for (int j = 0; j<5; j++) {

			io_hlp::write_to_val<double>(istr, vec_data[i][j]);

		}
	}

	// crop data read from file according to init_cfg:

	double total_thick = bl_thick_scale * init_cfg.ThickCoef;

	int total_nodes = 0;
	t_Vec3Dbl dr;

	bool thick_ok = false;

	for (int m = 0; m<nnodes_ns; m++) {

		matrix::base::minus<double, double>(xyz_vec[m], xyz_vec[0], dr);
		total_nodes++;
		if (dr.norm() > total_thick) {
			thick_ok = true;
			break;
		}

	}

	raw_profile.resize(total_nodes);

	for (int p = 0; p < total_nodes; p++) {

		raw_profile[p].set_xyz(xyz_vec[p]);

		raw_profile[p].u = vec_data[p][0];
		raw_profile[p].v = vec_data[p][1];
		raw_profile[p].w = vec_data[p][2];
		raw_profile[p].p = vec_data[p][3];
		raw_profile[p].t = vec_data[p][4];

	} ;

	if (!thick_ok) {
		wxLogError(_T("Error in ProfileNS::_read_crop_raw_profile_from_file Requested thickness \
						  is greater than provided from file, trying to proceed"));
	}
	else {
		// do some interpolation for the last point
		// we want thickness to be exactly total_thick
		// for now only xyz interpolated
		// TODO: linear interpolation for flow variables
		t_GeomPoint xyz1, xyz2, surf_xyz;

		t_Vec3Dbl rvec;

		surf_xyz = raw_profile[0].get_xyz();

		xyz1 = raw_profile[total_nodes-2].get_xyz();

		matrix::base::minus<double, double>(xyz1, surf_xyz, rvec);

		double r1 = rvec.norm();

		xyz2 = raw_profile[total_nodes - 1].get_xyz();

		matrix::base::minus<double, double>(xyz2, surf_xyz, rvec);

		double r2 = rvec.norm();

		double coef = (total_thick - r1) / (r2 - r1);

		if (coef<0.0 || coef>1.0)
			wxLogMessage(_T("Error in ProfileNS::_read_crop_raw_profile_from_file: interpolation coef is %lf, should be between 0 and 1"), coef);

		t_GeomPoint xyz_int = xyz1 + coef*(xyz2 - xyz1);

		raw_profile[total_nodes - 1].set_xyz(xyz_int);
	}

}

t_SqMat3Dbl _calc_jac_to_loc_rf(const std::vector<mf::t_Rec>& raw_profile, const t_Vec3Dbl& Ue_dir) {

	t_SqMat3Dbl jac;

	t_Rec bound_rec, surface_rec;

	surface_rec = raw_profile[0];
	bound_rec = raw_profile.back();

	// construct transformation matrix
	// construct normal to a surface - e2':
	// for now we need gridline j=const to be normal to surface
	// TODO: make interpolation of field to a normal for arbitrary grid
	t_Vec3Dbl e1, e2, e3, u_e;
	/*
	e2[0] = bound_rec.x - surface_rec.x;
	e2[1] = bound_rec.y - surface_rec.y;
	e2[2] = bound_rec.z - surface_rec.z;
	*/
	e2 = bound_rec.get_xyz() - surface_rec.get_xyz();
	//double norm = two_norm(e2);
	e2.normalize();
	// construct e1': along inviscid streamline
	// be sure e1'*e2'=0:
	// e1' = norm(Ue - (e2'*Ue));

	u_e = Ue_dir;
	e1 = u_e - vector::dot(e2, u_e)*e2;
	e1.normalize();
	// e3' = [e1' x e2']
	e3 = vector::cross(e1, e2);
	// ordinary orthogonal 
	// transformation matrix S
	// e' = eS
	/*
	_jacToLocalRF = e1[0], e2[0], e3[0],
	e1[1], e2[1], e3[1],
	e1[2], e2[2], e3[2];
	*/
	jac.set_col(0, e1);
	jac.set_col(1, e2);
	jac.set_col(2, e3);

	wxLogMessage(_T("bound rec xyz:%s"), bound_rec.get_xyz().to_wxstr());
	wxLogMessage(_T("bound rec UVW:%s"), bound_rec.get_uvw().to_wxstr());
	wxLogMessage(_T("Jac to loc rf:%s"), jac.to_wxstr());

	return jac;

}

void t_ProfileNS::initialize_extract(const std::string fname, const mf::t_ProfDataCfg& init_cfg) {

	std::vector<mf::t_Rec> raw_profile;
	t_Vec3Dbl Ue_dir;

	_read_crop_raw_profile_from_file(fname, init_cfg, raw_profile, Ue_dir);

	// #########initialize
	_resize(raw_profile.size());

	t_GeomPoint r_xyz_base, r_xyz;
	t_Vec3Dbl dr_xyz, r_ked;
	t_Vec3Dbl u_xyz, u_ked;

	r_xyz_base.set(raw_profile[0]);
	_xyz = r_xyz_base;

	t_SqMat3Dbl jac = _calc_jac_to_loc_rf(raw_profile, Ue_dir);
	t_SqMat3Dbl inv_jac = jac.inverse();

	const mf::t_FldParams& Params = _rDomain.get_mf_params();

	for (int j = 0; j<_nnodes; j++) {

		const mf::t_Rec& mf_rec = raw_profile[j];

		if (j == 0) {
			_y[j] = 0.0;
		}
		else {
			r_xyz.set(raw_profile[j]);
			r_xyz_base.set(raw_profile[j - 1]);
			matrix::base::minus<double, double>(r_xyz, r_xyz_base, dr_xyz);
			// convention : to make y = O(1)
			// it mul by factor sqrt(ReINF)
			// as it is for local parallel task we store only y assuming x=z=0
			// TODO: check r_ked[0] and r_ked[2] are small
			_y[j] = _y[j - 1] + dr_xyz.norm()*sqrt(Params.Re);
		}

		//r_xyz.set(mf_rec);
		//matrix::base::minus<double, double>(r_xyz, r_xyz_base, dr_xyz);
		//matrix::base::mat_mul<double, double>(inv_jac, dr_xyz, r_ked);

		u_xyz.set(mf_rec.u, mf_rec.v, mf_rec.w);
		matrix::base::mat_mul<double, double>(inv_jac, u_xyz, u_ked);

		_u[j] = u_ked[0];
		_v[j] = u_ked[1];
		_w[j] = u_ked[2];

		_p[j] = mf_rec.p;
		_t[j] = mf_rec.t;
		// rho is non-dim as follows
		double gMaMa = Params.Gamma*pow(Params.Mach, 2);
		_r[j] = mf_rec.p / mf_rec.t*gMaMa;
		_mu[j] = _rDomain.calc_viscosity(_t[j]);
	}

	_calc_derivs();

}

void t_ProfileNS::_initialize_interpolate(const t_GeomPoint& xyz, const mf::t_ProfDataCfg& init_cfg) {

	//IMPORTANT TODO: remove all "estimations" when 2-nd order interpolation done!!!
	const double& a_thick_coef = init_cfg.ThickCoef;
	const int& a_nnodes = init_cfg.NNodes;

	int num_bl_nodes_ns = _rDomain.estim_num_bl_nodes(xyz);

	if (num_bl_nodes_ns <= 0) {
		wxLogError(_T("ProfileNS: failed to estimate bl nodes number"));
		if (a_nnodes <= 0) ssuGENTHROW(_T("ProfileNS: failed to calculate nnodes"));
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

	double deta = prof_thick / double(_nnodes - 1);

	for (int j = 0; j<_nnodes; j++) {

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
		_y[j] = r_ked[1] * sqrt(Params.Re);
		// again we postulate v==0
		_u[j] = u_ked[0];
		_w[j] = u_ked[2];

		_p[j] = mf_rec.p;
		_t[j] = mf_rec.t;
		// rho is non-dim as follows
		double gMaMa = Params.Gamma*pow(Params.Mach, 2);
		_r[j] = mf_rec.p / mf_rec.t*gMaMa;
		_mu[j] = _rDomain.calc_viscosity(_t[j]);
	}

	int nnodes = get_nnodes();

	_calc_derivs();

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

	_calc_derivs();
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

		if (j==0){
			_y[j] = 0.0;
		}else{
			r_xyz.set(raw_profile[j]);
			r_xyz_base.set(raw_profile[j-1]);
			matrix::base::minus<double, double>(r_xyz, r_xyz_base, dr_xyz);
			_y[j] = _y[j-1] + dr_xyz.norm();
		}

		const mf::t_FldParams& Params = _rDomain.get_mf_params();

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

	_calc_derivs();

}


//--------------------------------------------------------------------------------t_ProfileMFLoc

t_ProfMFLoc::t_ProfMFLoc(const t_DomainBase& a_rDomain):t_ProfMF(a_rDomain){};

void t_ProfMFLoc::_initialize_interpolate(const t_GeomPoint& xyz, const mf::t_ProfDataCfg& init_cfg){

	//IMPORTANT TODO: remove all "estimations" when 2-nd order interpolation done!!!
	const double& a_thick_coef = init_cfg.ThickCoef;
	const int& a_nnodes = init_cfg.NNodes;

	int num_bl_nodes_ns = _rDomain.estim_num_bl_nodes(xyz);

	if (num_bl_nodes_ns<=0 ){
		wxLogError(_T("ProfileNS: failed to estimate bl nodes number"));
		if (a_nnodes<=0) ssuGENTHROW(_T("ProfileNS: failed to calculate nnodes"));
	}


	wxLogMessage(_T("ProfMFLoc: Check _initialize_interpolate"));
	a_nnodes>0 ? _resize(a_nnodes) : _resize(num_bl_nodes_ns);

	t_GeomPoint r_xyz_base;
	t_Vec3Dbl surf_norm;
	_rDomain.calc_surf_point(xyz, r_xyz_base, surf_norm);

	_xyz = r_xyz_base;
	double bl_thick = _rDomain.calc_bl_thick(_xyz);

	double cur_eta = bl_thick;

	double prof_thick = a_thick_coef*bl_thick;

	t_Vec3Dbl u_xyz, dr_xyz, u_ked;
	t_GeomPoint r_xyz;

	t_SqMat3Dbl jac = _rDomain.calc_jac_to_loc_rf(r_xyz_base);
	t_SqMat3Dbl inv_jac = jac.inverse();

	double deta = prof_thick/double(_nnodes-1);

	for (int j=0; j<_nnodes; j++){

		dr_xyz = j*deta*surf_norm;
		r_xyz = r_xyz_base + dr_xyz;
		const mf::t_Rec mf_rec = _rDomain.interpolate_to_point(r_xyz);

		u_xyz = mf_rec.get_uvw();

		matrix::base::mat_mul<double, double>(inv_jac, u_xyz, u_ked);

		const mf::t_FldParams& Params = _rDomain.get_mf_params();

		_y[j] = dr_xyz.norm();

		_u[j] = u_ked[0];
		_v[j] = u_ked[1];
		_w[j] = u_ked[2];

		_p[j] = mf_rec.p;
		_t[j] = mf_rec.t;
		// rho is non-dim as follows
		double gMaMa = Params.Gamma*pow(Params.Mach, 2);
		_r[j] = mf_rec.p/mf_rec.t*gMaMa;
		_mu[j]=_rDomain.calc_viscosity(_t[j]);
	}

	_calc_derivs();
}


void t_ProfMFLoc::_initialize_extract(const t_GeomPoint& xyz, const mf::t_ProfDataCfg& init_cfg){

	std::vector<mf::t_Rec> raw_profile;

	_rDomain.extract_profile_data(xyz, init_cfg, raw_profile);

	_resize(raw_profile.size());

	t_GeomPoint r_xyz_base, r_xyz;
	t_Vec3Dbl u_xyz, dr_xyz, u_ked;

	r_xyz_base.set(raw_profile[0]);
	_xyz = r_xyz_base;

	t_SqMat3Dbl jac = _rDomain.calc_jac_to_loc_rf(r_xyz_base);
	t_SqMat3Dbl inv_jac = jac.inverse();

	for (int j=0; j<_nnodes; j++){

		const mf::t_Rec& mf_rec = raw_profile[j];

		if (j==0){
			_y[j] = 0.0;
		}else{
			r_xyz.set(raw_profile[j]);
			r_xyz_base.set(raw_profile[j-1]);
			matrix::base::minus<double, double>(r_xyz, r_xyz_base, dr_xyz);
			_y[j] = _y[j-1] + dr_xyz.norm();
		}

		const mf::t_FldParams& Params = _rDomain.get_mf_params();

		u_xyz = mf_rec.get_uvw();
		matrix::base::mat_mul<double, double>(inv_jac, u_xyz, u_ked);

		_u[j] = u_ked[0];
		_v[j] = u_ked[1];
		_w[j] = u_ked[2];

		_p[j] = mf_rec.p;
		_t[j] = mf_rec.t;
		// rho is non-dim as follows
		double gMaMa = Params.Gamma*pow(Params.Mach, 2);
		_r[j] = mf_rec.p/mf_rec.t*gMaMa;
		_mu[j]=_rDomain.calc_viscosity(_t[j]);
	}

	_calc_derivs();

}

void t_ProfMFLoc::dump(const std::string& fname) const{
	std::wofstream fstr(&fname[0], std::ios::out);
	t_Rec rec;

	const t_Rec rec_out = get_bound_rec();

	const double rue_1 = 1.0/(rec_out.r*rec_out.u);

	fstr<<_T("y\tu\tu'\tu''\tt\tt'\tt''\tr\tmu\tmu'\tmu''\tw\tw'\tw''\tv\tMach\tdelta**\n");

	double mach, dd;

	mf::t_Rec mf_rec;

	std::vector<double> dd_v(get_nnodes());
	std::vector<double> mthick_v(get_nnodes());

	// for momentum thickness calculations
	for (int i=0; i<get_nnodes(); i++){

		rec = get_rec(i);

		dd_v[i] = rec.r *rec.u*rue_1*(1.0-rec.u/rec_out.u);
	}

	smat::integrate_over_range(_y, dd_v, mthick_v);

	for (int i=0; i<get_nnodes(); i++){

		rec = get_rec(i);

		mf_rec = rec.make_mf_rec();

		mach = _rDomain.calc_mach(mf_rec);

		fstr<<rec.y<<_T("\t")<<
			rec.u<<_T("\t")<<rec.u1<<_T("\t")<<rec.u2<<_T("\t")<<
			rec.t<<_T("\t")<<rec.t1<<_T("\t")<<rec.t2<<_T("\t")<<rec.r<<_T("\t")<<
			rec.mu<<_T("\t")<<rec.mu1<<_T("\t")<<rec.mu2<<_T("\t")<<
			rec.w<<_T("\t")<<rec.w1<<_T("\t")<<rec.w2<<_T("\t")<<rec.v<<_T("\t")<<
			mach<<_T("\t")<<mthick_v[i]<<("\n");


	}
};