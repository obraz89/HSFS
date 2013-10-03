#include "ProfileNS.h"
#include "math_operands.h"
#include "smooth.h"

using namespace mf;

t_ProfileNS::t_ProfileNS(const t_DomainBase& a_rDomain):t_Profile(0), _rDomain(a_rDomain){};

void t_ProfileNS::initialize(const t_GeomPoint xyz, double a_thick_coef, int a_nnodes){

	//if (a_nnodes<2){
	//	ssuGENTHROW(_T("ProfNS initialization error: less than 2 points!"));
	//}

	// so I need interpolator right now or tmp solution
	a_nnodes==0 ? _resize(_rDomain.estim_num_bl_nodes(xyz)) : _resize(a_nnodes);

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

	// VERY IMPORTANT TODO: how to choose _xScale
	// [ to keep "wavechars const" behavior in non-dim like for selfsim case]
	_xScale = _rDomain.calc_x_scale(_xyz);

	// SMOOTHING
	// FORTRAN CALLING CONVENTION

	int nnodes = get_nnodes();
	SMOOTH_3D_PROFILES(&_y[0], &_u[0], &nnodes, &_u1[0], &_u2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_w[0], &nnodes, &_w1[0], &_w2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_t[0], &nnodes, &_t1[0], &_t2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_mu[0], &nnodes, &_mu1[0], &_mu2[0]);
}


const mf::t_DomainBase& t_ProfileNS::getMFDomain() const{return _rDomain;};

double t_ProfileNS::get_x_Scale() const{return _xScale;};

int t_ProfileNS::get_bound_ind() const{return get_nnodes()-1;};

t_ProfileNS::t_Rec t_ProfileNS::get_bound_rec(){
	return get_rec(get_bound_ind());
}