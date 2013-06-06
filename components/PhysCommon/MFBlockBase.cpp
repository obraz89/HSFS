#include "stdafx.h"
#include "MFBlockBase.h"

static const double BL_BOUND_VELO_TOL = 0.01;

using namespace mf;

void t_Block::_allocate(){
	if (!_allocated){

		_fld = new t_Rec**[Nx];
		for(int i=0;i<Nx;i++) 
		{
			_fld[i] = new t_Rec*[Ny];
			for(int j=0;j<Ny;j++)
			{
				_fld[i][j] = new t_Rec[Nz];
			};
		}
		_allocated = true;
	}
};

t_Block::t_Block():_allocated(false){};

t_Block::~t_Block()
{
	for (int i=0; i<Nx; i++)
	{
		for (int j=0; j<Ny; j++) delete[] _fld[i][j];
		delete[] _fld[i];
	}
};

int t_Block::get_Nx() const{return Nx;};
int t_Block::get_Ny() const{return Ny;};
int t_Block::get_Nz() const{return Nz;};

const t_Rec& t_Block::get_rec(const t_BlkInd ind) const{
	return _fld[ind.i][ind.j][ind.k];
};

void t_Block::print_entry(const t_BlkInd ind) const
{
	int i = ind.i;
	int j = ind.j;
	int k = ind.k;

	const t_Rec& p = _fld[i][j][k];
	std::cout<<i<<";"<<j<<";"<<k
		<<"\nX:"<<p.x<<"  Y:"<<p.y<<"  Z:"<<p.z
		<<"\nU:"<<p.u<<"  V:"<<p.v<<"  W:"<<p.w
		<<"\nP:"<<p.p<<"  T:"<<p.t<<" Ro:"<<p.r<<"\n";
};

t_SqMat3Dbl mf::t_Block::calc_jac_to_loc_rf(const t_BlkInd ind) const{

	t_SqMat3Dbl jac;

	const t_Rec& bound_rec = get_rec(t_BlkInd(ind.i, get_bound_index(ind), ind.k));
	const t_Rec& surface_rec = get_rec(t_BlkInd(ind.i, 0, ind.k));

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
	jac[0] = e1;
	jac[1] = e2;
	jac[2] = e3;
	
	return jac;

}

double t_Block::calc_enthalpy(const t_BlkInd ind) const
{

	double cp_t, v_2;

	const t_FldParams& mf_prms = get_mf_params();

	const t_Rec& ptr = get_rec(ind);

	cp_t = ptr.t/((mf_prms.Gamma-1.0)*mf_prms.Mach*mf_prms.Mach);

	v_2 = 0.5*(pow(ptr.u,2.0)+pow(ptr.v,2.0)+pow(ptr.w,2.0));

	return (cp_t + v_2);

};

double t_Block::calc_enthalpy_freestream() const{

	const t_FldParams& mf_prms = get_mf_params();

	// TODO: is this always correct?
	// is always |u|=1
	return (1.0/((mf_prms.Gamma-1.0)*mf_prms.Mach*mf_prms.Mach) + 0.5);

};


// for future 
double t_Block::calc_delta(const t_BlkInd ind) const{return 0.0;};

// calc nondim viscosity
double t_Block::calc_viscosity(const t_BlkInd ind) const{

	const t_Rec& rRec = get_rec(ind);

	const t_FldParams& mf_prms = get_mf_params();

	if (mf_prms.ViscType==t_ViscType::ViscPower){

		return pow(rRec.t, mf_prms.Mju_pow);

	}
	else{

		double t_suth = mf_prms.T_mju/mf_prms.T_inf;

		return pow(rRec.t, 1.5)*(1.0+t_suth)/(rRec.t+t_suth);

	}
};

// dimensional infinity cinematic viscosity
double t_Block::calc_cin_visc_inf_dim() const{

	const t_FldParams& mf_prms = get_mf_params();

	double u_inf_dim = calc_u_inf();

	return u_inf_dim*mf_prms.L_ref/mf_prms.Re;

};

double t_Block::calc_cin_visc_dim(const t_BlkInd ind) const{

	const t_FldParams& mf_prms = get_mf_params();

	double nju_inf = calc_cin_visc_inf_dim();

	return nju_inf*calc_viscosity(ind)/get_rec(ind).r;

};

double t_Block::calc_c_dim(const t_BlkInd ind) const{

	const t_FldParams& mf_prms = get_mf_params();

	double t_dim = get_rec(ind).t*mf_prms.T_inf;

	return sqrt(mf_prms.Gamma*mf_prms.R_Gas*t_dim/mf_prms.Mol_weight);

};

double t_Block::calc_u_inf() const{

	const t_FldParams& mf_prms = get_mf_params();

	return mf_prms.Mach*
		sqrt(mf_prms.Gamma*mf_prms.R_Gas*mf_prms.T_inf/mf_prms.Mol_weight);
};

double t_Block::calc_mach(const t_BlkInd ind) const{

	const t_Rec& rRec = get_rec(ind);

	const t_FldParams& mf_prms = get_mf_params();

	double vAbs = sqrt(pow(rRec.u,2.0)+pow(rRec.v,2.0)+pow(rRec.w,2.0));

	return mf_prms.Mach*vAbs/sqrt(rRec.t);

}

int t_Block::get_bound_index(const t_BlkInd ind) const
{
	return get_bound_ind_velo(ind);
};

int t_Block::get_bound_ind_enth(const t_BlkInd ind) const{
		/* for spec. xz_plane_ind. = {i,0,k} computes {i,j,k},
	   j is boundary layer border index;
	   enthalpy criterion is used	*/
int j=0;

// TODO: this is not correct in multiblock
// IMPORTANT TODO: when restructure is done, replace with calc_enth_freestream
double h_inf = calc_enthalpy(t_BlkInd(0,Ny/2,0));

bool brd_rchd = false;
while(!brd_rchd)
{
	brd_rchd = true;
	j++;
	for(int dj=0;dj<20;dj++)
	{
		// TODO : fix this !!! Param?
		if (fabs(calc_enthalpy(t_BlkInd(ind,0,dj,0)) - h_inf)/h_inf>0.001) 
		{
			brd_rchd = false;
			break;
		};
	}
};
return (j);
}

int t_Block::get_bound_ind_velo(const t_BlkInd a_ind) const{
	t_BlkInd cur_ind = a_ind;
	t_BlkInd nxt_ind(cur_ind, 0,1,0);
	double dy, du_dy_wall, du_dy_cur;
	t_Vec3Dbl u, du;

	u = get_rec(nxt_ind).get_uvw();
	dy = calc_distance(cur_ind, nxt_ind);
	du_dy_wall = u.norm()/dy;

	for (int j=1; j++; j<Ny-1){
		cur_ind.j=j;
		nxt_ind.j=j+1;
		du = (get_rec(nxt_ind) - get_rec(cur_ind)).get_uvw();
		dy = calc_distance(nxt_ind, cur_ind);
		du_dy_cur = du.norm()/dy;
		if (du_dy_cur<BL_BOUND_VELO_TOL*du_dy_wall){
			return j;
		}
	}
	return -1;
}

double t_Block::calc_gridline_distance
(ALONG_LINE along_line, t_BlkInd from, t_BlkInd to) const{
	int* pChgInd=NULL;
	int n = 0;
	t_BlkInd cur_ind = from, tmp, prev_ind;
	switch(along_line){
		case ALONG_LINE::I:
			pChgInd = &(cur_ind.i);
			if (to.i<from.i){
				tmp = to;
				to = from;
				from = tmp;
			}
			n = to.i - from.i;
			break;
		case ALONG_LINE::J:
			pChgInd = &(cur_ind.j);
			if (to.j<from.j){
				tmp = to;
				to = from;
				from = tmp;
			}
			n = to.j - from.j;
			break;
		case ALONG_LINE::K:
			if (to.k<from.k){
				tmp = to;
				to = from;
				from = tmp;
			}
			pChgInd = &(cur_ind.k);
			n = to.k - from.k;
			break;
	}
	double distance = 0.0;
	for (int m=0; m<n; m++){
		prev_ind = cur_ind;
		(*pChgInd)++;
		distance+=calc_distance(prev_ind, cur_ind);
	}
	return distance;
};

double t_Block::calc_distance(const t_BlkInd a, const t_BlkInd b) const{
	const t_Rec& r1 = get_rec(a);
	const t_Rec& r2 = get_rec(b);
	double dx = r1.x - r2.x;
	double dy = r1.y - r2.y;
	double dz = r1.z - r2.z;
	return sqrt(dx*dx+dy*dy+dz*dz);
};

void t_Block::_calc_dir_vec(t_VecDbl& vec, t_BlkInd ind, ALONG_LINE along_line) const{
	if (!_check_ind(ind)){
		ssuGENTHROW(_("MF Index Out of Range"));
	};
	t_Rec dd;

	switch (along_line)
	{
	case I:
		if (ind.i==0){
			t_BlkInd adj_ind(ind, 1, 0, 0);
			dd = get_rec(adj_ind) - get_rec(ind);
		}else{
			if (ind.i==Nx-1){
				t_BlkInd adj_ind(ind, -1, 0, 0);
				dd = get_rec(ind) - get_rec(adj_ind);
			}else{
				t_BlkInd adj_ind_p(ind, 1, 0, 0);
				t_BlkInd adj_ind_m(ind, -1, 0, 0);
				dd = get_rec(adj_ind_p) - get_rec(adj_ind_m);
			};
		};
		break;
	case J:
		if (ind.j==0){
			t_BlkInd adj_ind(ind, 0, 1, 0);
			dd = get_rec(adj_ind) - get_rec(ind);
		}else{
			if (ind.j==Ny-1){
				t_BlkInd adj_ind(ind, 0, -1, 0);
				dd = get_rec(ind) - get_rec(adj_ind);
			}else{
				t_BlkInd adj_ind_p(ind, 0, 1, 0);
				t_BlkInd adj_ind_m(ind, 0, -1, 0);
				dd = get_rec(adj_ind_p) - get_rec(adj_ind_m);
			};
		};
		break;
	case K:
		if (ind.k==0){
			t_BlkInd adj_ind(ind, 0, 0, 1);
			dd = get_rec(adj_ind) - get_rec(ind);
		}else{
			if (ind.i==Nz-1){
				t_BlkInd adj_ind(ind, 0, 0, -1);
				dd = get_rec(ind) - get_rec(adj_ind);
			}else{
				t_BlkInd adj_ind_p(ind, 0, 0, 1);
				t_BlkInd adj_ind_m(ind, 0, 0, -1);
				dd = get_rec(adj_ind_p) - get_rec(adj_ind_m);
			}
		}
		break;
	}
	vec = dd.get_xyz();
	vec.normalize();
};

void t_Block::_calc_gridline_dirs
(t_VecDbl &i_dir, t_VecDbl& j_dir, t_VecDbl& k_dir, t_BlkInd ind) const{
	_calc_dir_vec(i_dir, ind, I);
	_calc_dir_vec(j_dir, ind, J);
	_calc_dir_vec(k_dir, ind, K);
};

// this is really interesting
bool t_Block::_is_inside
(const t_Vec3Dbl& point, t_BlkInd diag1, t_BlkInd diag2) const{
	// base index
	// to get "canonical" diagonal
	// [i,j,k] <~> [i+1, j+1, k+1]
	t_BlkInd base = _get_base_ind(diag1, diag2);
	// in most cases di=dj=dk=1
	int di = abs(diag1.i - diag2.i);
	int dj = abs(diag1.j - diag2.j);
	int dk = abs(diag1.k - diag2.k);
	t_VecDbl rbase = get_rec(base).get_xyz();
	t_VecDbl r = point - rbase;
	t_VecDbl bvecs[3], norms;
	t_VecDbl& bvec_i = bvecs[0]; 
	t_VecDbl& bvec_j = bvecs[1];
	t_VecDbl& bvec_k = bvecs[2];
	bvec_i = get_rec(t_BlkInd(base, di, 0,  0)).get_xyz() - rbase;
	bvec_j = get_rec(t_BlkInd(base, 0,  dj, 0)).get_xyz() - rbase;
	bvec_k = get_rec(t_BlkInd(base, 0,  0, dk)).get_xyz() - rbase;
	for (int i=0; i<3; i++){
		norms[i] = bvecs[i].norm();
	};
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			bvecs[i][j]/=norms[i];
		};
	};
	t_SqMatrix<double> gramm, inv_gramm;
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			gramm[i][j] = vector::dot(bvecs[i], bvecs[j]);
		};
	};
	inv_gramm = gramm.inverse();
	t_Vec3Dbl rhs(vector::dot(r, bvec_i),
				  vector::dot(r, bvec_j),
				  vector::dot(r, bvec_k));
				
	t_VecDbl coefs;
	coefs = inv_gramm*rhs;
	bool inside = true;
	for (int i=0; i<3; i++){
		if ((coefs[i]>norms[i])||(coefs[i]<0.0)) inside=false;
	};
	return inside;
};

bool t_Block::_check_ind(const t_BlkInd& ind) const{

	return ((ind.i<Nx)&&(ind.j<Ny)&&(ind.k<Nz)&&
			(ind.i>=0)&&(ind.j>=0)&&(ind.k>=0));
};

t_BlkInd t_Block::_get_nearest_node
(const t_GeomPoint& point, t_BlkInd diag1, t_BlkInd diag2) const{
	t_BlkInd ret;
	// IMPORTANT TODO: fix this
	double min_dst = calc_distance(t_BlkInd(0,0,0),t_BlkInd(Nx-1, 0,0));
	t_BlkInd base = _get_base_ind(diag1, diag2);
	// in most cases di=dj=dk=1
	int di = abs(diag1.i - diag2.i);
	int dj = abs(diag1.j - diag2.j);
	int dk = abs(diag1.k - diag2.k);
	for (int i=0;i<2; i++){
		for (int j=0; j<2; j++){
			for (int k=0; k<2; k++){
				t_BlkInd vertex(base, di*i, dj*j, dk*k);
				t_VecDbl dr;
				dr = get_rec(vertex).get_xyz() - point;
				if (dr.norm()<min_dst){
					ret = vertex;
				};
			};
		};
	};
	return ret;
};

t_BlkInd t_Block::_get_base_ind(t_BlkInd diag1, t_BlkInd diag2) const{
	t_BlkInd base;
	base.i = (diag1.i<diag2.i) ? diag1.i : diag2.i;
	base.j = (diag1.j<diag2.j) ? diag1.j : diag2.j;
	base.k = (diag1.k<diag2.k) ? diag1.k : diag2.k;	
	return base;
};

t_BlkInd t_Block::_get_nearest_index_loc
(t_BlkInd start_from, const t_GeomPoint& point) const{
	if (!_check_ind(start_from)){
		ssuGENTHROW(_("MF:Bad start index in _get_nearest_BlkInd_loc"));
	};
	t_BlkInd cur_ind = start_from;
	double cur_dst=1.0, prev_dst;
	t_VecDbl i_dir, j_dir, k_dir;
	int di, dj, dk;
	do{
		prev_dst = cur_dst;
		t_VecDbl dir;
		dir = point - get_rec(cur_ind).get_xyz();
		dir.normalize();
		_calc_gridline_dirs(i_dir, j_dir, k_dir, cur_ind);
		di = (vector::dot(i_dir, dir)>0.0) ? 1 : -1;
		dj = (vector::dot(j_dir, dir)>0.0) ? 1 : -1;
		dk = (vector::dot(k_dir, dir)>0.0) ? 1 : -1;
		t_BlkInd ind(cur_ind, di, dj, dk);
		if (_check_ind(ind)){
			t_VecDbl dd = get_rec(ind).get_xyz() - point;
			double dst = dd.norm();
			if (dst<prev_dst){
				cur_dst = dst;
				cur_ind = ind;
			};
		}else{
			// "manual" adjustment (check all 27 adj nodes)
			for (int i=-1; i<2; i++)
				for (int j=-1; j<2; j++)
					for (int k=-1; k<2; k++){
						t_BlkInd ind(cur_ind, i, j, k);
						if (_check_ind(ind)){
							double dst = (get_rec(ind).get_xyz() - point).norm();
							if (dst<prev_dst){
								cur_dst = dst;
								cur_ind = ind;
							} 
						}
					}
			
		};
	} while (cur_dst<prev_dst);
	return cur_ind;
};

t_BlkInd t_Block::get_nearest_index_loc(t_BlkInd start_from, t_Rec rec) const{
	return _get_nearest_index_loc(start_from, rec.get_xyz());
};

t_BlkInd t_Block::get_nearest_index_loc(t_BlkInd start_from, t_GeomPoint geom_point) const{
	return _get_nearest_index_loc(start_from, geom_point);
};
//old (bad) version
t_BlkInd t_Block::get_nearest_index_raw(t_GeomPoint geom_point) const{

	// old crappy mess 
	/*
	int pos_lft = 0;
	int pos_rgt = Nx-1;
	while(pos_rgt-pos_lft>1) {
		ind_nrst.i = (pos_lft + pos_rgt)/2;
		const t_Rec& p = get_rec(ind_nrst); //_rFldMF.fld[ind_nrst.i][ind_nrst.j][ind_nrst.k];
		x_cmp = p.x;
		if (x>=x_cmp) 
			pos_lft = ind_nrst.i;
		else 
			pos_rgt = ind_nrst.i;
	};
	if (abs(x - get_rec(t_BlkInd(pos_lft,0,0)).x)
		<abs(x - get_rec(t_BlkInd(pos_rgt,0,0)).x))
		ind_nrst.i = pos_lft;
	else
		ind_nrst.i = pos_rgt;

	double z_cmp = 0.;
	pos_lft = 0;
	pos_rgt = Nz-1;
	while(pos_rgt-pos_lft>1) {
		ind_nrst.k = (pos_lft + pos_rgt)/2;
		const t_Rec& p = get_rec(ind_nrst);
		z_cmp = p.z;
		if (z>=z_cmp) 
			pos_lft = ind_nrst.k;
		else 
			pos_rgt = ind_nrst.k;
	};
	if (abs(z - get_rec(t_BlkInd(ind_nrst.i,0,pos_lft)).z)
		<abs(z - get_rec(t_BlkInd(ind_nrst.i,0,pos_rgt)).z))
		ind_nrst.k = pos_lft;
	else
		ind_nrst.k = pos_rgt;

	double y_cmp = 0.;
	pos_lft = 0;
	pos_rgt = Ny-1;
	while(pos_rgt-pos_lft>1) {
		ind_nrst.j = (pos_lft + pos_rgt)/2;
		const t_Rec& p = get_rec(ind_nrst);
		y_cmp = p.y;
		if (y>=y_cmp) 
			pos_lft = ind_nrst.j;
		else 
			pos_rgt = ind_nrst.j;
	};
	ind_nrst.j = pos_lft;
	*/

	t_BlkInd ind_nrst;

	double cur_dst, min_dst = get_mf_params().L_ref;


	t_BlkInd cur_ind;
	t_GeomPoint surf_node_xyz;


	// TODO: is this too slow?

	for (int i=0; i<get_Nx(); i++){
		for (int k=0; k<get_Nz(); k++)
		{

			cur_ind.i = i;
			cur_ind.j = 0;
			cur_ind.k = k;

			surf_node_xyz = get_rec(cur_ind).get_xyz();
			cur_dst = (geom_point - surf_node_xyz).norm();

			if (cur_dst<=min_dst)
			{
				ind_nrst = cur_ind;
				min_dst = cur_dst;
			}
		}
	}

	return ind_nrst;
}

t_BlkInd t_Block::get_nearest_index_raw(t_Rec rec) const{
	return get_nearest_index_raw(rec.get_xyz());
};

t_Rec t_Block::interpolate_to_point(t_GeomPoint point) const{
	// TODO: good interpolation !!!!!!!!!!!!!!!!!!!!!!!!!
	// for now the worst variant:
	t_BlkInd nrst_raw = get_nearest_index_raw(point);
	t_BlkInd nrst_ind = nrst_raw;//get_nearest_index_loc(nrst_raw, point);
	t_Rec ret = get_rec(nrst_ind);
	ret.set_xyz(point);
	return ret;
}
//##############################################
// old croosflow low level mess
/*
void t_Block::get_merid_distrib(const int i_ind) const
{
	std::ostringstream _to_file_name;
	std::string filename;
	_to_file_name<<"output/merid_i="<<i_ind<<"_al="<<int(Alpha*58.0)<<".dat";
	filename = _to_file_name.str();
	std::ofstream _to_file(&filename[0]);
	_to_file<<"merid_angle \t w_e \t u_e \t m_e \t re1_e \n";
	for (int k=0;k<nz;k++){
	int bound_ind = get_bound_index(i_ind,k);
	std::cout<<k<<" : "<<bound_ind<<"\n";
	const Rec& bound_rec = _fld[i_ind][bound_ind][k];
	double gMaMa = Params.Gamma*Params.Mach*Params.Mach;
	double re1_e = gMaMa*bound_rec.p*bound_rec.u/pow(bound_rec.t,1.75)*Params.Re/Params.L_ref;
	double m_e = bound_rec.u*Params.Mach/sqrt(bound_rec.t);
	_to_file<<double(k)/double(nz)*3.1415<<" \t "
			<<bound_rec.w<<" \t "
			<<bound_rec.u<<" \t "
			<<m_e<<" \t "<<re1_e<<"\n";
	}
	_to_file.close();
}

void t_Block::get_cf_profile(std::vector<ProfileRec>& prof, const int i_ind, const int k_ind) const
{
	int bound_ind = get_bound_index(i_ind,k_ind);
	const Rec& bound_rec = _fld[i_ind][bound_ind][k_ind];
	double u_e = bound_rec.u;
	double w_e = bound_rec.w;
	double ang = atan(w_e/u_e);
	prof.clear();
	for (int j=0; j<ny; j++)
	{
		const Rec& cur_rec = _fld[i_ind][j][k_ind];
		double u_new = cur_rec.u*cos(ang) + cur_rec.w*sin(ang);
		double w_new = -cur_rec.u*sin(ang) + cur_rec.w*cos(ang);
		prof.push_back(ProfileRec(cur_rec.y, w_new));
	};
};

//--------------------PRIVATE PART
void t_Block::get_cf_prof_rotated
(const int i_ind, const int k_ind, const double ang, std::vector<double>& profile) const
{
	int bound_ind = get_bound_index(i_ind,k_ind);
	profile.clear();

	for (int j=0; j<bound_ind; j++)
	{
		const Rec& cur_rec = fld[i_ind][j][k_ind];
		profile.push_back(-cur_rec.u*sin(ang) + cur_rec.w*cos(ang));
	};
	return;
};

int t_Block::get_cf_zero_index(const std::vector<double>& profile) const
{
	for(int j=1; j<profile.size()-1; j++)
		if (profile[j-1]*profile[j+1]<0.0)
			return j;
	return 0;
};

int t_Block::get_cf_infl_index(const std::vector<double>& profile) const
{
	std::vector<double> w2yy_prof(profile.size(), 0); 
	for (int j=1; j<w2yy_prof.size()-2; j++)
		w2yy_prof[j] =	profile[j+1] - 2*profile[j] + profile[j-1];
	w2yy_prof[0] = w2yy_prof[1];
	w2yy_prof[w2yy_prof.size()-1] = w2yy_prof[w2yy_prof.size()-2];
	for (int j=1; j<w2yy_prof.size()-1; j++)
		if (w2yy_prof[j-1]*w2yy_prof[j+1]<0.0) 
			return j;
	return profile.size()-1;
};

double t_Block::get_cf_wave_dir(const int i_ind, const int k_ind) const
{
	const int bound_ind = get_bound_index(i_ind, k_ind);
	const Rec& ptr = fld[i_ind][bound_ind][k_ind];
	std::cout<<"ext streamline tan(psi) = "<<ptr.w/ptr.u<<"\n";
	double ang_def = atan(ptr.w/ptr.u);
	double ang_min = 0.0;
	double ang_max = 2.0*ang_def;
	std::ofstream f2file("output/crossflow_rotated.dat");
	double ang_res;
	int resid=1000;
	std::vector<double> profile(0);
	int niter=100;
	for (int i=0;i<niter;i++)
	{
		double ang = ang_min + double(i)/double(niter)*(ang_max - ang_min);
		get_cf_prof_rotated(i_ind, k_ind,ang,profile);
		int zero = get_cf_zero_index(profile);
		int infl = get_cf_infl_index(profile);
		int cur_resid = abs(zero-infl);
		if (cur_resid<resid )
		{
			resid = cur_resid;
			ang_res = ang;
		}
	};
	//get_cf_prof_rotated(i_ind, k_ind,ang_res,profile);
	for (int j=0; j<profile.size(); j++) f2file<<j<<"	"<<profile[j]<<"\n"; 
	std::cout<<"result angle,resid:"<<ang_res<<" ; "<<resid<<"\n";
	f2file.close();
	// wave phase speed dir in the inviscid streamline coord frame:
	//return (ang_res-ang_def);
	// wave phase speed dir in grid coord frame:
	return ang_res;
};


void t_Block::get_profiles(const int i_ind, const int k_ind) const{
	std::ostringstream to_file_name;
	std::string filename;
	to_file_name<<"output/profiles_i="<<i_ind<<"_k="
		<<k_ind<<"_al="<<int(Alpha*58.0)<<"_new.dat";
	filename = to_file_name.str();
	std::ofstream s_to_file(&filename[0]);
	s_to_file<<" y[ft] \t u[ft/s] \t u/a[] \t w[ft/s] \n";

	double u_inf_dim = t_Block::Mach*sqrt(t_Block::Gamma*8.31*t_Block::T_inf/0.029);
	std::vector<ProfileRec> w_prof;
	get_cf_profile(w_prof, i_ind, k_ind);
	for (int j=0; j<ny; j++){
		const Rec& cur_rec = fld[i_ind][j][k_ind];
	s_to_file<< cur_rec.y*t_Block::L_ref/0.3048 << " \t "	// [ft]
			 << cur_rec.u*u_inf_dim/0.3048<< " \t "			// [ft/s]
			 << t_Block::Mach*cur_rec.u/sqrt(cur_rec.t)<< " \t "  //[]
			 << w_prof[j].val*u_inf_dim/0.3048<<"\n";				// [ft/s]
	};
	s_to_file.close();
//##############################################################

};
*/