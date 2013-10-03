#include "stdafx.h"
#include "MFHSDomain.h"

using namespace mf;
using namespace mfhs;

// TODO: this works only for ortho grids
void t_Block::calc_surf_norm(const t_BlkInd a_ind, t_Vec3Dbl& norm) const{

	t_BlkInd surf_ind = a_ind;
	surf_ind.j = 0;
	t_GeomPoint surf_point = get_rec(surf_ind).get_xyz();
	t_GeomPoint uppr_point = get_rec(t_BlkInd(surf_ind,0,1,0)).get_xyz();
	norm = (uppr_point - surf_point).normalize();
};

void t_Block::calc_surf_norm(const mf::t_GeomPoint a_xyz, t_Vec3Dbl& norm) const{

	t_BlkInd nrst_ind = get_nearest_ind_surf(a_xyz);
	calc_surf_norm(nrst_ind, norm);
};

void t_Block::calc_surf_point(const t_GeomPoint& a_xyz, t_GeomPoint& surf_point, t_Vec3Dbl& norm) const{

	double min_avg_dst = 1.0E+20;//BIG_INFINITE;
	t_BlkInd base_ind;
	t_GeomPoint cur_patch[3], min_patch[3];
	t_BlkInd min_dst_inds[3], cur_indxs[3];
	t_Vec3Dbl dst;

	for (int i=0; i<Nx-1; i++)
	for (int k=0; k<Nz-1; k++)
	for (int p=0; p<2 ; p++)
	{
		cur_indxs[0] = t_BlkInd(i,0,k);
		cur_patch[0].set(get_rec(cur_indxs[0]));

		cur_indxs[1] = (p==0) ? t_BlkInd(i+1, 0, k) : t_BlkInd(i, 0, k+1);
		cur_patch[1].set(get_rec(cur_indxs[1]));

		cur_indxs[2] = t_BlkInd(i+1,0,k+1);
		cur_patch[2].set(get_rec(cur_indxs[2]));

		double dd=0.0;

		for (int m=0; m<3; m++)
		{
			matrix::base::minus<double, double>(a_xyz, cur_patch[m], dst);
			dd+=dst.norm();
		}

		dd/=3.0;

		if (dd<min_avg_dst)
		{
			min_avg_dst = dd;
			for (int m=0; m<3; m++) {
				min_dst_inds[m] = cur_indxs[m];
				min_patch[m] = cur_patch[m];
			}
			
		}


	}

	// construct normal to a found patch
	/*
	t_Vec3Dbl e1, e2;
	for (int m=0; m<3; m++) cur_patch[m] = get_rec(min_dst_inds[m]).get_xyz();

	e1 = (cur_patch[1] - cur_patch[0]).normalize();
	e2 = (cur_patch[2] - cur_patch[0]).normalize();
	norm = vector::cross(e1, e2);
	*/
	calc_surf_norm(a_xyz, norm);
	// get a surface projection point
	// for now just go to nearest surface node

	double min_node_dst = 1.0E+20;  // TODO:BIG_INFINITE
	for (int m=0; m<3; m++){
		matrix::base::minus<double, double>(a_xyz, cur_patch[m], dst);
		if (dst.norm()<min_node_dst){
			surf_point.set(get_rec(min_dst_inds[m]));
		}
	}
};

// IMPORTANT TODO: better way to calculate x scale?
double t_Block::calc_x_scale(const mf::t_GeomPoint& xyz) const{
	t_BlkInd nrst_ind = get_nearest_index_raw(xyz);
	t_BlkInd start_ind = nrst_ind;
	start_ind.i = 0;
	return calc_distance(start_ind, nrst_ind);
}


t_SqMat3Dbl t_Block::calc_jac_to_loc_rf(const t_BlkInd ind) const{

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
	jac.set_col(0, e1);
	jac.set_col(1, e2);
	jac.set_col(2, e3);
	
	return jac;

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
	t_Vec3Dbl rbase = get_rec(base).get_xyz();
	t_Vec3Dbl r = point - rbase;
	t_Vec3Dbl bvecs[3], norms;
	t_Vec3Dbl& bvec_i = bvecs[0]; 
	t_Vec3Dbl& bvec_j = bvecs[1];
	t_Vec3Dbl& bvec_k = bvecs[2];
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
	t_SqMat3Dbl gramm, inv_gramm;
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

	double cur_dst, min_dst = 1.0E+20;


	t_BlkInd cur_ind;
	t_Vec3Dbl cur_rvec;
	t_GeomPoint cur_node_xyz;


	// TODO: is it too slow to iterate over entire block?

	for (int i=0; i<get_Nx(); i++){
	  for (int j=0; j<get_Ny(); j++)
		for (int k=0; k<get_Nz(); k++)
		{

			cur_ind.i = i;
			cur_ind.j = j;
			cur_ind.k = k;

			cur_node_xyz.set(get_rec(cur_ind));
			matrix::base::minus<double, double>(geom_point, cur_node_xyz, cur_rvec);
			cur_dst = cur_rvec.norm();

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

t_BlkInd t_Block::get_nearest_ind_surf(mf::t_GeomPoint geom_point) const{

	t_BlkInd ind_nrst;

	double cur_dst, min_dst = 1.0E+20;


	t_BlkInd cur_ind;
	t_Vec3Dbl cur_rvec;
	t_GeomPoint cur_node_xyz;

	// do not iterate over j

	for (int i=0; i<get_Nx(); i++){
		//for (int j=0; j<get_Ny(); j++)
			for (int k=0; k<get_Nz(); k++)
			{

				cur_ind.i = i;
				cur_ind.j = 0;
				cur_ind.k = k;

				cur_node_xyz.set(get_rec(cur_ind));
				matrix::base::minus<double, double>(geom_point, cur_node_xyz, cur_rvec);
				cur_dst = cur_rvec.norm();

				if (cur_dst<=min_dst)
				{
					ind_nrst = cur_ind;
					min_dst = cur_dst;
				}
			}
	}

	return ind_nrst;

}

t_Rec t_Block::interpolate_to_point(t_GeomPoint point) const{

	t_Rec ret;
	// TODO: good interpolation !!!!!!!!!!!!!!!!!!!!!!!!!
	// for now the worst variant:
	t_BlkInd nrst_surf = get_nearest_ind_surf(point);

	t_GeomPoint p1, p2;
	t_Vec3Dbl dd, r1, r2;
	t_BlkInd bt_ind = nrst_surf;
	t_BlkInd up_ind = nrst_surf;

	bool point_inside = false;

	for (int j=0; j<Ny-1; j++){
		bt_ind.j = j;
		p1.set(get_rec(bt_ind));
		up_ind.j = j+1;
		p2.set(get_rec(up_ind));

		matrix::base::minus<double, double>(p2, p1, dd);
		matrix::base::minus<double, double>(point, p1, r1);
		matrix::base::minus<double, double>(point, p2, r2);

		double cprod1, cprod2;

		cprod1 = vector::dot(dd, r1);
		cprod2 = vector::dot(dd, r2);

		if ((cprod1>=0)&&(cprod2<=0)){

			point_inside = true;

			double d, l1, l2, d1, d2, c;
			d = dd.norm();
			l1 = r1.norm();
			l2 = r2.norm();
			c = (l2*l2-l1*l1+d*d)/d;
			d1 = 0.5*(d-c);
			d2 = 0.5*(d+c);
			double a = 1.0 - d1/d;
			double b = 1.0 - d2/d;
#ifdef _DEBUG
			if (a<0.0 || a>1.0 || b<0. || b>1.)
				wxLogMessage(_T("Two point mf interpolation failed!\n"));
#endif

			const t_Rec& rec1 = get_rec(bt_ind);
			const t_Rec& rec2 = get_rec(up_ind);
			
			#define SET_RET(o) ret.o## = a*rec1.o## + b*rec2.o##;
			SET_RET(u);SET_RET(v);SET_RET(w);
			SET_RET(p);SET_RET(t);SET_RET(r);
			#undef SET_RET
			break;
		}
	}

	if (!point_inside){
		t_GeomPoint p;
		t_Vec3Dbl dr;
		t_BlkInd cur_ind = nrst_surf;
		double min_dst = 1.0E+20;
		t_BlkInd min_dst_ind = nrst_surf;
		double cur_dst;
		for (int j=0; j<Ny; j++){
			cur_ind.j = j;
			p.set(get_rec(cur_ind));
			matrix::base::minus<double, double>(p, point, dr);
			cur_dst = dr.norm();
			if (cur_dst<min_dst){
				min_dst = cur_dst;
				min_dst_ind = cur_ind;
			}

		}
		ret = get_rec(min_dst_ind);
	}

	ret.set_xyz(point);
	return ret;
}

bool t_Block::is_point_inside(const t_GeomPoint& xyz) const{
	// IMPORTANT TODO: implement!!!
	t_BlkInd nrst_ind = get_nearest_index_raw(xyz);
	return (nrst_ind.i!=0) && (nrst_ind.i!=Nx-1) &&
		   (nrst_ind.k!=0) && (nrst_ind.k!=Nz-1);

};
