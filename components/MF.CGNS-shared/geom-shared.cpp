#include "stdafx.h"
#include "cgns_structs.h"

#include <fstream>

using namespace mf::cg;

void TDomain::set_face_iters(int iZone, int iFace, int& is, int& ie, 
							 int& js, int& je, int& ks, int& ke) const{

	const TZone& blk = Zones[iZone-1];

	// first set full blk ranges

	is = blk.is;
	ie = blk.ie;

	js = blk.js;
	je = blk.je;

	get_k_range(iZone, ks, ke);

	switch(iFace){
		case(faceXmin):
			ie = is;
			break;
		case(faceXmax):
			is = ie;
			break;
		case(faceYmin):
			je = js;
			break;
		case(faceYmax):
			js = je;
			break;
		case(faceZmin):
			ke = ks;
			break;
		case(faceZmax):
			ks = ke;
			break;
	}
	
};


t_ZoneNode TDomain::_get_nrst_node_surf(const t_GeomPoint& point) const{

	t_Rec cur_rec;
	t_GeomPoint cur_xyz;
	t_Vec3Dbl dr;
	double cur_dst, min_dst = HUGE_VAL;

	t_ZoneNode min_znode;

	for( int iZone = 1;  iZone <= nZones;  ++iZone )
	{
		TZone& blk = Zones[iZone-1];

		// TODO: 6 is BRICK_FACENUM
		for (int iface = 0; iface<6; iface++){
			if (_is_face_of_bcwall_type(blk.Faces[iface].szBCFamName)){
				// we are on viscous wall
				int i_s, i_e, j_s, j_e, k_s, k_e;

				set_face_iters(iZone, iface, i_s, i_e, j_s, j_e, k_s, k_e);

		
				for (int i=i_s; i<=i_e; i++)
					for (int j=j_s; j<=j_e; j++)
						for (int k=k_s; k<=k_e; k++){

							get_rec(blk, i, j, k, cur_rec);
							cur_xyz.set(cur_rec);
							matrix::base::minus<double, double>(cur_xyz, point, dr);
							cur_dst = dr.norm();

							if (cur_dst<min_dst){
								min_dst = cur_dst;
								min_znode.iZone = iZone;
								min_znode.iNode.set(i,j,k);
								min_znode.iFacePos = iface;
							}

						}

			} 

		}	// ~iterate over blk BCWall faces

	}

	if (min_znode.iZone<0) wxLogError(_T("Failed to find nearest point on viscous wall!"));

	return min_znode;
}

t_ZoneNode TDomain::_get_nrst_node_surf(const t_ZoneNode& src_node) const{
	
	t_Rec rec;
	get_rec(src_node, rec);
	return _get_nrst_node_surf(rec.get_xyz());
}

t_ZoneNode TDomain::_get_nrst_node_raw(const t_GeomPoint& geom_point) const{
	// iterate over all nodes in all zones to find nearest node

	mf::t_Rec cur_rec;
	double cur_dst, min_dst = 1.0E+20;

	t_Vec3Dbl cur_rvec;
	t_GeomPoint cur_node_xyz;
	t_ZoneNode cur_znode, min_dst_znode;

	for (int iZone=1; iZone<=nZones; iZone++){

		TZone& blk = Zones[iZone-1];

		int ks, ke;
		get_k_range(iZone, ks, ke);

		for (int i=blk.is; i<=blk.ie; i++){
			for (int j=blk.js; j<=blk.je; j++)
				for (int k=ks; k<=ke; k++)
				{

					cur_znode.iZone = iZone;
					cur_znode.iNode.set(i,j,k);
					get_rec(cur_znode, cur_rec);
					cur_node_xyz.set(cur_rec);
					matrix::base::minus<double, double>(geom_point, cur_node_xyz, cur_rvec);
					cur_dst = cur_rvec.norm();

					if (cur_dst<=min_dst)
					{
						min_dst_znode = cur_znode;
						min_dst = cur_dst;
					}
				}
		}


	}
	
	return min_dst_znode;
}

void TDomain::calc_surf_point(const t_GeomPoint& a_xyz, t_GeomPoint& surf_point, t_Vec3Dbl& norm) const{

	t_ZoneNode surf_node = _get_nrst_node_surf(a_xyz);
	t_Rec surf_rec;
	get_rec(surf_node, surf_rec);
	surf_point.set(surf_rec);
	calc_surf_norm(surf_node, norm);
	//wxLogError(_T("calc_surf_point not implemented"));
	//ssuTHROW(t_GenException, _T("Not Implemented"));

};

// Simple and relatively fast(~Nz faster than looping through entire Domain) 
// but not generally correct way to find node(i.e. find zone and local index inside Zone) 
// that is the nearest to argument geom point;
// First search for nearest surface point, then go along norm gridline and find min dst
// with 2-point interpolation, if failed, try to rescue with 1-point(zero-order) interpolation

mf::t_Rec TDomain::_interpolate_to_point_surf_raw(const t_GeomPoint& point) const{

	t_Rec ret, cur_rec;

	t_ZoneNode znode_mdsurf = _get_nrst_node_surf(point);
	const t_BlkInd& ind_md = znode_mdsurf.iNode;

	// nearest surface node found, now make interpolation
	// TODO: good interpolation

	//=======================================================

	int i_s = ind_md.i;
	int i_e = i_s;

	int j_s = ind_md.j;
	int j_e = j_s;

	int k_s = ind_md.k;
	int k_e = k_s;

	const TZone& the_zone = Zones[znode_mdsurf.iZone -1];

	bool is_var_x = false;
	bool is_var_y = false;
	bool is_var_z = false;

	switch(znode_mdsurf.iFacePos){
			case(faceXmin):
			case(faceXmax):
				// vary i and keep j,k fixed
				is_var_x = true;
				i_s = the_zone.is;
				i_e = the_zone.ie - 1;
				break;
			case(faceYmin):
			case(faceYmax):
				// vary j and keep i,k fixed
				is_var_y = true;
				j_s = the_zone.js;
				j_e = the_zone.je - 1;
				break;
			case(faceZmin):
			case(faceZmax):
				// vary k and keep i,j fixed
				is_var_z = true;
				k_s = the_zone.ks;
				k_e = the_zone.ke - 1;
				break;
			default:
				wxLogError(
					_T("Wrong face pos in interpolate to point, Zone=%d, facepos=%d"), 
					znode_mdsurf.iZone, znode_mdsurf.iFacePos);
	}

	t_GeomPoint p1, p2;
	t_Vec3Dbl dd, r1, r2;
	t_BlkInd bt_ind = ind_md;
	t_BlkInd up_ind = ind_md;

	bool point_inside = false;

	t_Rec rec1, rec2;

	for (int i=i_s; i<=i_e; i++)
		for (int j=j_s; j<=j_e; j++)
			for (int k=k_s; k<=k_e; k++){

				bt_ind.set(i,j,k);
				TDomain::get_rec(the_zone, bt_ind, cur_rec);
				p1.set(cur_rec);

				up_ind = bt_ind;

				if (is_var_x) {
					up_ind.i = i+1;
					TDomain::get_rec(the_zone, up_ind, cur_rec);
				}
				if (is_var_y){
					up_ind.j = j+1;
					TDomain::get_rec(the_zone, up_ind, cur_rec);
				}
				if (is_var_z){
					up_ind.k = k+1;
					TDomain::get_rec(the_zone, up_ind, cur_rec);

				}

				p2.set(cur_rec);

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
					c = (l2*l2-l1*l1)/d;
					d1 = 0.5*(d-c);
					d2 = 0.5*(d+c);
					double a = 1.0 - d1/d;
					double b = 1.0 - d2/d;
#ifdef _DEBUG
					if (a<0.0 || a>1.0 || b<0. || b>1.)
						wxLogMessage(_T("Two point mf interpolation failed!\n"));
#endif

					//const t_Rec& rec1 = get_rec(bt_ind);
					//const t_Rec& rec2 = get_rec(up_ind);
					TDomain::get_rec(the_zone, bt_ind, rec1);
					TDomain::get_rec(the_zone, up_ind, rec2);

#define SET_RET(o) ret.o = a*rec1.o + b*rec2.o
					SET_RET(u);SET_RET(v);SET_RET(w);
					SET_RET(p);SET_RET(t);SET_RET(r);
#undef SET_RET
					// break looping
					goto stop;
				}
			}

stop:;

			if (!point_inside){
				// now check robustness
				wxLogError(
					_T("Failed to interpolate to point - point is outside. Zoneid = %d\n Rescue with zero order..."), 
					znode_mdsurf.iZone);
				// rescue with zero order

				t_GeomPoint p;
				t_Vec3Dbl dr;
				t_BlkInd cur_ind = ind_md;
				double min_dst = 1.0E+20;
				t_BlkInd min_dst_ind = ind_md;
				double cur_dst;

				for (int i=i_s; i<=i_e; i++)
					for (int j=j_s; j<=j_e; j++)
						for (int k=k_s; k<=k_e; k++){
							get_rec(the_zone, i, j, k, cur_rec);
							p.set(cur_rec);
							matrix::base::minus<double, double>(p, point, dr);
							cur_dst = dr.norm();
							if (cur_dst<min_dst){
								min_dst = cur_dst;
								min_dst_ind.set(i,j,k);
							}
						}
						TDomain::get_rec(the_zone, min_dst_ind, ret);

			}

			ret.set_xyz(point);
			return ret;
};

bool TDomain::_is_point_inside(const t_GeomPoint& xyz) const{

	wxLogError(_T("MF.CGNS-shared: is_point_inside not yet implemented correctly, using bounding box!"));
	return bbox.is_pnt_inside(xyz);

}

bool TDomain::is_point_inside(const t_GeomPoint& xyz) const{

	return _is_point_inside(xyz);

}

// TODO: current version works for ortho near-surface grids only
void TDomain::calc_surf_norm(const t_ZoneNode& surf_node, t_Vec3Dbl& norm) const{

	t_GeomPoint wall_xyz, fld_xyz;
	t_Rec wall_rec, fld_rec;

	const TZone& blk = Zones[surf_node.iZone-1];

	const int i_base = surf_node.iNode.i;
	const int j_base = surf_node.iNode.j;
	const int k_base = surf_node.iNode.k;

	switch(surf_node.iFacePos){
		case(faceXmin):
			get_rec(blk, i_base+1, j_base, k_base, fld_rec);
			break;
		case(faceXmax):
			get_rec(blk, i_base-1, j_base, k_base, fld_rec);
			break;
		case(faceYmin):
			get_rec(blk, i_base, j_base+1, k_base, fld_rec);
			break;
		case(faceYmax):
			get_rec(blk, i_base, j_base-1, k_base, fld_rec);
			break;
		default:
			wxLogError(_T("Mf Domain error: Failed to calc norm, zone=%d, facepos=%d"), surf_node.iZone, surf_node.iFacePos);
			break;
	}

	get_rec(blk, i_base, j_base, k_base, wall_rec);
	wall_xyz.set(wall_rec);

	fld_xyz.set(fld_rec);
	matrix::base::minus<double, double>(fld_xyz, wall_xyz, norm);
	norm.normalize();

};

double TDomain::calc_x_scale(const t_GeomPoint& xyz) const{
	wxLogMessage(_T("MF Domain Warning: LRef scale used as x scale"));
	return get_mf_params().L_ref;
};

// TODO: works only for grids nearly orthogonal to viscous wall 
// TODO: inout bl_thick is for y_ref !!!
void TDomain::_calc_bl_thick_vderiv(const t_GeomPoint& xyz, double& bl_thick, 
					t_ZoneNode& surf_znode, t_ZoneNode& outer_znode) const{

	surf_znode = _get_nrst_node_surf(xyz);

	t_ZoneGrdLine grd_line(*this, surf_znode);

	int i_s = grd_line.ind_s.i;
	int i_e = grd_line.ind_e.i;

	int j_s = grd_line.ind_s.j;
	int j_e = grd_line.ind_e.j;

	int k_s = grd_line.ind_s.k;
	int k_e = grd_line.ind_e.k;

	const TZone& blk = *grd_line.zne;

	int i = i_s;
	int j = j_s;
	int k = k_s;

	int di = grd_line.di;
	int dj = grd_line.dj;
	int dk = grd_line.dk;

	// if du_dy increases search max
	// use du_dy wall otherwise
	double du_dy_max = 0.0;
	double du_dy;

	mf::t_Rec cur_rec, nxt_rec;
	t_GeomPoint wall_xyz, cur_xyz, nxt_xyz, dr;
	t_Vec3Dbl cur_uvw, nxt_uvw, du;

	// wall rec is used to calculate bl_thick
	get_rec(blk, i, j, k, cur_rec);
	wall_xyz.set(cur_rec);

	// first find max deriv
	// then find reference point where dudy = eps*dudy_max : y = dd
	// then find outer record y_out = thick_coef*dd
	bool searching_max = true;
	bool searching_ref = true;

	double y_ref;

	t_ZoneNode ref_znode;

	for (int m=0; m<grd_line.nnodes-1; m++) 
	{

		get_rec(blk, i, j, k, cur_rec);

		cur_xyz.set(cur_rec);

		cur_uvw.set(cur_rec.u, cur_rec.v, cur_rec.w);

		i+=di;
		j+=dj;
		k+=dk;

		get_rec(blk, i, j, k, nxt_rec);

		nxt_xyz.set(nxt_rec);

		nxt_uvw.set(nxt_rec.u, nxt_rec.v, nxt_rec.w);

		matrix::base::minus<double, double>(nxt_xyz, cur_xyz, dr);	
		matrix::base::minus<double, double>(nxt_uvw, cur_uvw, du);

		du_dy = abs(du.norm()/dr.norm());

		if (searching_max){

			if (du_dy<du_dy_max) searching_max = false;

			du_dy_max = du_dy;

			continue;

		}else{

			double eps = _profile_cfg.DerivThreshold;

			if (searching_ref){

				if (du_dy<eps*du_dy_max){
					searching_ref = false;

					ref_znode.iZone = surf_znode.iZone;
					ref_znode.iNode.set(i,j,k);

					matrix::base::minus<double, double>(cur_xyz, wall_xyz, dr);

					y_ref = dr.norm();
					bl_thick = y_ref;	

					continue;
				}
				

			}else{

				matrix::base::minus<double, double>(cur_xyz, wall_xyz, dr);

				double cur_dd = dr.norm();
				// IMPORTANT TODO: passing thick coef !!!
				if (cur_dd>=_profile_cfg.ThickCoefDefault*y_ref){

					outer_znode.iZone = surf_znode.iZone;
					outer_znode.iNode.set(i,j,k);

					return;

				}

			}
/*
			double eps = _profile_cfg.DerivThreshold;
			if (du_dy<eps*du_dy_max){

				// IMPORTANT TODO: 
				// here we found not "outer" znode but znode at dd away from wall
				// outer should be at max_bl_thick_coef*dd !!!
				outer_znode.iZone = surf_znode.iZone;
				outer_znode.iNode.set(i,j,k);

				matrix::base::minus<double, double>(cur_xyz, wall_xyz, dr);

				bl_thick = dr.norm();

				break;

			};
*/

		}

	};

	wxLogError(_T("Error: Zone boundary reached while searching BL bound & outer rec"));

};

// TODO : works only for grids nearly orthogonal to viscous wall
void TDomain::_calc_bl_thick_enthalpy(const t_GeomPoint& xyz, double& bl_thick, 
									t_ZoneNode& surf_znode, t_ZoneNode& outer_znode) const{


	surf_znode = _get_nrst_node_surf(xyz);

	t_ZoneGrdLine grd_line(*this, surf_znode);

	int i_s = grd_line.ind_s.i;
	int i_e = grd_line.ind_e.i;

	int j_s = grd_line.ind_s.j;
	int j_e = grd_line.ind_e.j;

	int k_s = grd_line.ind_s.k;
	int k_e = grd_line.ind_e.k;

	const TZone& blk = *grd_line.zne;

	mf::t_Rec cur_rec, nxt_rec;
	t_GeomPoint wall_xyz, cur_xyz, nxt_xyz, dr;
	double h_cur, h_inf;

	// wall rec is used to calculate bl_thick
	get_rec(blk, i_s, j_s, k_s, cur_rec);
	wall_xyz.set(cur_rec);

	h_inf = calc_enthalpy_freestream();

	// march from outside of bl

	int i = i_e;
	int j = j_e;
	int k = k_e;

	int di = -grd_line.di;
	int dj = -grd_line.dj;
	int dk = -grd_line.dk;

	for (int m=0; m<grd_line.nnodes; m++) 
	{

		get_rec(blk, i, j, k, cur_rec);

		h_cur = calc_enthalpy(cur_rec);

		double tol = abs(h_cur/h_inf -1);

		if (tol>_profile_cfg.DerivThreshold){

			outer_znode.iZone = surf_znode.iZone;
			outer_znode.iNode.set(i,j,k);

			cur_xyz.set(cur_rec);

			matrix::base::minus<double, double>(cur_xyz, wall_xyz, dr);

			bl_thick = dr.norm();

			break;
		}

		i+=di;
		j+=dj;
		k+=dk;


	};

};

void TDomain::_calc_bl_thick_full_gridline(const t_GeomPoint& xyz, double& bl_thick, 
				t_ZoneNode& surf_znode, t_ZoneNode& outer_znode) const{

	surf_znode = _get_nrst_node_surf(xyz);

	t_ZoneGrdLine grd_line(*this, surf_znode);

	outer_znode.iZone = surf_znode.iZone;
	outer_znode.iNode = grd_line.ind_e;

	mf::t_Rec rec1, rec2;

	get_rec(surf_znode, rec1);
	get_rec(outer_znode, rec2);

	t_GeomPoint r1, r2;
	t_Vec3Dbl dr;

	r1.set(rec1);
	r2.set(rec2);

	matrix::base::minus<double, double>(r1, r2, dr);
	bl_thick = dr.norm();
	
}

inline void TDomain::_calc_bl_thick(const t_GeomPoint& xyz, double& bl_thick, 
									t_ZoneNode& surf_znode, t_ZoneNode& outer_znode) const{


	switch (_profile_cfg.BLThickCalcType)
	{
	case t_BLThickCalcType::BLTHICK_BY_VDERIV:
		_calc_bl_thick_vderiv(xyz,bl_thick, surf_znode, outer_znode);
		break;

	case t_BLThickCalcType::BLTHICK_BY_ENTHALPY:
		_calc_bl_thick_enthalpy(xyz,bl_thick, surf_znode, outer_znode);
		break;
	case t_BLThickCalcType::FULL_DATA:
		_calc_bl_thick_full_gridline(xyz,bl_thick, surf_znode, outer_znode);
		break;

	default:
		wxString msg(_T("Provided BL Thick Calc Type not implemented"));
		wxLogError(msg); ssuGENTHROW(msg);
	}


}


double TDomain::calc_bl_thick(const t_GeomPoint& xyz) const{
	
	t_ZoneNode surf_znode, outer_znode;
	double bl_thick;

	_calc_bl_thick(xyz, bl_thick, surf_znode, outer_znode);
	return bl_thick;

};

void TDomain::calc_nearest_surf_rec(const t_GeomPoint& xyz, t_Rec& surf_rec) const{

	t_ZoneNode surf_znode;	

	surf_znode = _get_nrst_node_surf(xyz);

	t_BlkInd& ind = surf_znode.iNode;

	get_rec(Zones[surf_znode.iZone-1], ind.i, ind.j, ind.k, surf_rec);

};

void TDomain::calc_nearest_inviscid_rec(const t_GeomPoint& xyz, t_Rec& outer_rec) const{

	t_ZoneNode surf_znode, outer_znode;
	double bl_thick;

	_calc_bl_thick(xyz,bl_thick, surf_znode, outer_znode);
	t_BlkInd& ind = outer_znode.iNode;

	get_rec(Zones[outer_znode.iZone-1], ind.i, ind.j, ind.k, outer_rec);

};

int TDomain::estim_num_bl_nodes(const t_GeomPoint& xyz) const{

	t_ZoneNode surf_znode, outer_znode;
	double bl_thick;

	_calc_bl_thick(xyz,bl_thick, surf_znode, outer_znode);
	
	// TODO; ensure that only 1 index is different for surf_znode
	// and outer_znode

	int N = abs(outer_znode.iNode.i - surf_znode.iNode.i) +
		    abs(outer_znode.iNode.j - surf_znode.iNode.j) +
			abs(outer_znode.iNode.k - surf_znode.iNode.k);

	return N;
}

t_SqMat3Dbl TDomain::calc_jac_to_loc_rf(const t_GeomPoint& xyz) const{
	
	t_ZoneNode bound_znode, surf_znode;
	double bl_thick;
	_calc_bl_thick(xyz, bl_thick, surf_znode, bound_znode);

	t_SqMat3Dbl jac;

	t_Rec bound_rec, surface_rec;
	get_rec(bound_znode, bound_rec);
	get_rec(surf_znode, surface_rec);

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

// tmp func to print enthalpy profiles for enthalpy criteria
// IMPORTANT TODO: the profile coord Y is along gridline (may be curvilinear)
void TDomain::dump_full_enthalpy_profile(const mf::t_GeomPoint& xyz, int pid) const{

	t_ZoneNode surf_znode = _get_nrst_node_surf(xyz);

	int i_s = surf_znode.iNode.i;
	int i_e = i_s;

	int j_s = surf_znode.iNode.j;
	int j_e = j_s;

	int k_s = surf_znode.iNode.k;
	int k_e = k_s;

	const TZone& blk = Zones[surf_znode.iZone-1];

	int i = i_s;
	int j = j_s;
	int k = k_s;

	int di = 0, dj =0 , dk = 0;

	switch(surf_znode.iFacePos){
		case(faceXmin):
			di = 1;
			i_s = blk.is;
			i_e = blk.ie;
			break;
		case(faceXmax):
			di = -1;
			i_s = blk.ie;
			i_e = blk.is;
			break;
		case(faceYmin):
			dj = 1;
			j_s = blk.js;
			j_e = blk.je;
			break;
		case(faceYmax):
			dj = -1;
			j_s = blk.js;
			j_e = blk.je;
			break;
		case(faceZmin):
			dk = 1;
			get_k_range(surf_znode.iZone, k_s, k_e);
			break;
		case(faceZmax):
			dk=-1;
			get_k_range(surf_znode.iZone, k_e, k_s);
			break;
	}

	char fname[33];
	sprintf(fname, "%s/profile_enthalpy_%d.dat", hsstab::OUTPUT_DIR.ToAscii(), pid);
	std::ofstream ofstr(fname);

	mf::t_Rec mf_rec;
	mf::t_GeomPoint cur_xyz,prev_xyz, surf_xyz, dr;

	get_rec(blk, i, j, k, mf_rec);
	surf_xyz.set(mf_rec);

	// Prints all but 1 points 
	int M = (i_e==i_s ? 0 : abs(i_e - i_s)+1) + 
		    (j_e==j_s ? 0 : abs(j_e - j_s)+1) + 
			(k_e==k_s ? 0 : abs(k_e - k_s)+1);

	double s = 0.0;

	cur_xyz = surf_xyz; prev_xyz = surf_xyz;

	for (int m=0; m<M; m++){

		get_rec(blk, i, j, k, mf_rec);
		cur_xyz.set(mf_rec);
		matrix::base::minus<double, double>(cur_xyz, prev_xyz, dr);
		s+= dr.norm();
		double enth = calc_enthalpy(mf_rec);
		ofstr<<s<<"\t"<<enth<<"\n";

		prev_xyz = cur_xyz;

		i+=di; 
		j+=dj;
		k+=dk;
	}

};
