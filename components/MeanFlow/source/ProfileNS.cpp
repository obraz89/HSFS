#include "ProfileNS.h"

t_ProfileNS::t_ProfileNS(const t_MeanFlow& a_rFld):t_Profile(a_rFld, 0){};
t_ProfileNS::~t_ProfileNS(){};

void t_ProfileNS::initialize(int a_i , int a_k){
	_i = a_i;
	_k = a_k;
	_bl_bound_ind = _rFld.get_bound_index(a_i, a_k);
	// empiric  - 3 thickn of BL to be used in stab comps
	// !!!
	double bl_thick = _rFld.calc_distance(t_Index(a_i, _bl_bound_ind, a_k), 
										  t_Index(a_i, 0, a_k));
	double prof_thick = 3.0*bl_thick;
	double cur_y = bl_thick;
	int cur_y_ind = _bl_bound_ind;
	const t_MFParams& Params = _rFld.base_params();
	while((cur_y<prof_thick)&&(cur_y_ind<Params.Ny)){
		cur_y_ind++;
		cur_y = _rFld.calc_distance(t_Index(a_i, cur_y_ind, a_k), 
									t_Index(a_i, 0, a_k));
	};
	this->_resize(cur_y_ind);
	const t_MeanFlow::t_Rec& bound_rec = _rFld.get_rec(a_i, cur_y_ind, a_k);
	const t_MeanFlow::t_Rec& surface_rec = _rFld.get_rec(a_i, 0, a_k);

// construct transformation matrix
	// construct normal to a surface - e2':
	// for now we need gridline j=const to be normal to surface
	// TODO: make interpolation of field to a normal for arbitrary grid
	dense_vector<double, _param_fixed_vec3> e1,e2,e3, u_e;
	e2[0] = bound_rec.x - surface_rec.x;
	e2[1] = bound_rec.y - surface_rec.y;
	e2[2] = bound_rec.z - surface_rec.z;
	double norm = two_norm(e2);
	for (int i=0; i<mtl::vector::size(e2); i++){
		e2[i] = e2[i]/norm;
	};
	// construct e1': along inviscid streamline
	// be sure e1'*e2'=0:
	// e1' = norm(Ue - (e2'*Ue));

	// I hate this syntax
	u_e = bound_rec.u, bound_rec.v, bound_rec.w;
	e1 = u_e - vector::dot(e2, u_e)*e2;
	norm = two_norm(e1);
	for (int i=0; i<mtl::vector::size(e1); i++){
		e1[i] = e1[i]/norm;
	};
	// e3' = [e1' x e2']
	e3 = vector::cross(e1, e2);
	// ordinary orthogonal 
	// transformation matrix S
	// e' = eS
	_jacToLocalRF = e1[0], e2[0], e3[0],
					e1[1], e2[1], e3[1],
					e1[2], e2[2], e3[2];
	//
// transform vector fields and coordinates
// to a new reference frame
	dense_vector<double, _param_fixed_vec3> u_xyz, u_ked, r_xyz, dr_xyz, r_ked, r_xyz_base;
	t_SqMat3 inv_jac;
	inv_jac = mtl::matrix::inv(_jacToLocalRF);
	r_xyz_base = surface_rec.x, surface_rec.y, surface_rec.z;
	for (int j=0; j<this->size(); j++){
		const t_MeanFlow::t_Rec& mf_rec = _rFld.get_rec(a_i, j, a_k);
		r_xyz = mf_rec.x, mf_rec.y, mf_rec.z;
		dr_xyz = r_xyz - r_xyz_base;
		r_ked = inv_jac*dr_xyz;
		u_xyz = mf_rec.u, mf_rec.v, mf_rec.w;
		u_ked = inv_jac*u_xyz;
		
		// convention : to make y = O(1)
		// it mul by factor sqrt(ReINF)
		// as it is for local parallel task we store only y assuming x=z=0
		// TODO: (good debug check is to assure that)
		_y[j] = r_ked[1]*sqrt(Params.Re);
		// again we postulate v==0
		_u[j] = u_ked[0];
		_w[j] = u_ked[2];

		_p[j] = mf_rec.p;
		_t[j] = mf_rec.t;
		// rho is non-dim as follows
		double gMaMa = Params.Gamma*pow(Params.Mach, 2);
		_r[j] = mf_rec.p/mf_rec.t*gMaMa;
		_mu[j]=_rFld.calc_viscosity(a_i, j, a_k);
	}

	// distance along grid line !!!
	t_MeanFlow::t_GridIndex from(0, 0, a_k), to(a_i, 0, a_k);
	xDist = _rFld.calc_gridline_distance(t_MeanFlow::ALONG_LINE::I, from, to);

	// SMOOTHING
	// FORTRAN CALLING CONVENTION

	int nnodes = this->size();
	SMOOTH_3D_PROFILES(&_y[0], &_u[0], &nnodes, &_u1[0], &_u2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_w[0], &nnodes, &_w1[0], &_w2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_t[0], &nnodes, &_t1[0], &_t2[0]);
	SMOOTH_3D_PROFILES(&_y[0], &_mu[0], &nnodes, &_mu1[0], &_mu2[0]);
	// status check
}

t_ProfileNS::t_Rec t_ProfileNS::get_bl_bound_rec(){
	return get_rec(_bl_bound_ind);
}