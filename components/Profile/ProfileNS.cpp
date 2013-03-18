#include "ProfileNS.h"
#include "math_operands.h"
#include "smooth.h"

using namespace mf;

t_ProfileNS::t_ProfileNS(const t_Block& a_rBlk):t_Profile(0), _rBlk(a_rBlk){};

void t_ProfileNS::initialize(const t_BlkInd a_ind, double a_thick_coef){

	_mf_ind = a_ind;
	_mf_ind.j = 0;

	_bl_bound_ind = _rBlk.get_bound_index(a_ind);

	const int a_i = a_ind.i;
	const int a_k = a_ind.k;

	double bl_thick = 
		_rBlk.calc_distance(t_BlkInd(a_i, _bl_bound_ind, a_k), 
						    t_BlkInd(a_i, 0            , a_k));

	double cur_y = bl_thick;
	int cur_y_ind = _bl_bound_ind;
	double prof_thick = a_thick_coef*bl_thick;

	const int Nx = _rBlk.get_Nx();
	const int Ny = _rBlk.get_Ny();
	const int Nz = _rBlk.get_Nz();

	while((cur_y<prof_thick)&&(cur_y_ind<Ny)){
		cur_y_ind++;
		cur_y = _rBlk.calc_distance(t_BlkInd(a_i, cur_y_ind, a_k), 
									t_BlkInd(a_i, 0        , a_k));
	};

	// TODO: check for ind>Ny overflow, throw something

	this->_resize(cur_y_ind);

	const mf::t_Rec& bound_rec = _rBlk.get_rec(t_BlkInd(a_i, cur_y_ind, a_k));
	const mf::t_Rec& surface_rec = _rBlk.get_rec(t_BlkInd(a_i, 0, a_k));

// construct transformation matrix
	// construct normal to a surface - e2':
	// for now we need gridline j=const to be normal to surface
	// TODO: make interpolation of field to a normal for arbitrary grid
	t_Vec3Dbl e1,e2,e3, u_e;
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

	u_e = bound_rec.get_uvw();
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
	t_SqMat3Dbl& mtr_jac = _mtr.jac;
	mtr_jac[0] = e1;
	mtr_jac[1] = e2;
	mtr_jac[2] = e3;
	//
// transform vector fields and coordinates
// to a new reference frame
	t_Vec3Dbl u_xyz, u_ked;
	t_GeomPoint r_xyz_base, r_xyz, dr_xyz, r_ked;
	t_SqMat3Dbl inv_jac;

	inv_jac = _mtr.inv_jac();

	r_xyz_base = surface_rec.get_xyz();

	for (int j=0; j<get_nnodes(); j++){
		const mf::t_Rec& mf_rec = _rBlk.get_rec(t_BlkInd(a_i, j, a_k));
		r_xyz = mf_rec.get_xyz();

		// TODO: What The Fuck??? I want 
		// dr_xyz = r_xyz - r_xyz_base
		(t_Vec3Dbl&)dr_xyz = r_xyz - r_xyz_base;
		(t_Vec3Dbl&)r_ked = inv_jac*dr_xyz;

		u_xyz = mf_rec.get_uvw();
		u_ked = inv_jac*u_xyz;
		
		// convention : to make y = O(1)
		// it mul by factor sqrt(ReINF)
		// as it is for local parallel task we store only y assuming x=z=0
		// TODO: (good debug check is to assure that)
		const mf::t_FldParams& Params = _rBlk.get_mf_params();
		_y[j] = r_ked[1]*sqrt(Params.Re);
		// again we postulate v==0
		_u[j] = u_ked[0];
		_w[j] = u_ked[2];

		_p[j] = mf_rec.p;
		_t[j] = mf_rec.t;
		// rho is non-dim as follows
		double gMaMa = Params.Gamma*pow(Params.Mach, 2);
		_r[j] = mf_rec.p/mf_rec.t*gMaMa;
		_mu[j]=_rBlk.calc_viscosity(t_BlkInd(a_i, j, a_k));
	}

	// distance along grid line !!!
	mf::t_BlkInd from(0, 0, a_k), to(a_i, 0, a_k);
	_xDist = _rBlk.calc_gridline_distance(t_Block::ALONG_LINE::I, from, to);

	// SMOOTHING
	// FORTRAN CALLING CONVENTION

	int nnodes = get_nnodes();
	SMOOTH_3D_PROFILES(&_y[0], &_u[0], &nnodes, &_u1[0], &_u2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_w[0], &nnodes, &_w1[0], &_w2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_t[0], &nnodes, &_t1[0], &_t2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_mu[0], &nnodes, &_mu1[0], &_mu2[0]);
	// status check
}


const mf::t_Block& t_ProfileNS::getBlk() const{return _rBlk;};

mf::t_BlkInd t_ProfileNS::getMFInd() const{return _mf_ind;};

double t_ProfileNS::get_xDist() const{return _xDist;};

int t_ProfileNS::get_bl_bound_ind() const{return _bl_bound_ind;};

const mf::t_Mtr& t_ProfileNS::get_mtr() const{return _mtr;};

t_SqMat3Dbl t_ProfileNS::getJac() const{return _mtr.jac;};

t_ProfileNS::t_Rec t_ProfileNS::get_bl_bound_rec(){
	return get_rec(_bl_bound_ind);
}