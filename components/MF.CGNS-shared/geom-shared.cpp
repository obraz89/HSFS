#include "stdafx.h"
#include "cgns_structs.h"

#include <fstream>

using namespace mf;
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
								min_znode.iFacePos = static_cast<TZoneFacePos>(iface);
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

// calculate derivative of velocity vs y like:
// U={u,v,w}
// d|U|/dy
// |dU/dy|
// du/dy
// etc
double TDomain::_calc_specifid_velo_deriv_abs(const std::vector<t_ZoneNode>& data_grdline, 
	int ind, t_VeloDerivType vd_type) const{

	int is = 0;
	int ie = data_grdline.size() - 2;

	if (ind<is || ind>ie) {

		wxLogMessage(_T("Error: _calc_specifid_velo_deriv: index is out of range: i=%d"), ind);

		return -1.0;

	}

	mf::t_Rec cur_rec, nxt_rec;
	t_GeomPoint wall_xyz, cur_xyz, nxt_xyz, dr;
	t_Vec3Dbl cur_uvw, nxt_uvw, du;

	get_rec(data_grdline[ind], cur_rec);

	cur_xyz.set(cur_rec);

	cur_uvw.set(cur_rec.u, cur_rec.v, cur_rec.w);

	double cur_uvw_abs = cur_uvw.norm();

	get_rec(data_grdline[ind + 1], nxt_rec);

	nxt_xyz.set(nxt_rec);

	nxt_uvw.set(nxt_rec.u, nxt_rec.v, nxt_rec.w);

	double nxt_uvw_abs = nxt_uvw.norm();

	matrix::base::minus<double, double>(nxt_xyz, cur_xyz, dr);
	matrix::base::minus<double, double>(nxt_uvw, cur_uvw, du);

	double dd = dr.norm();

	if (vd_type == t_VeloDerivType::VD_ABS) return abs(cur_uvw_abs - nxt_uvw_abs)/dd;

	if (vd_type == t_VeloDerivType::VD_VEC_ABS) return du.norm()/dd;

	if (vd_type == t_VeloDerivType::VD_X_ABS) return abs(du[0]) / dd;

	if (vd_type == t_VeloDerivType::VD_TAU_VEC_ABS){

		t_ZoneNode surf_znode = _get_nrst_node_surf(data_grdline[0]);

		t_Vec3Dbl surf_norm;
		calc_surf_norm(surf_znode, surf_norm);

		double un_cur = vector::dot(cur_uvw, surf_norm);
		double ut_cur = sqrt(cur_uvw_abs*cur_uvw_abs - un_cur*un_cur);

		double un_nxt = vector::dot(nxt_uvw, surf_norm);
		double ut_nxt = sqrt(nxt_uvw_abs*nxt_uvw_abs - un_nxt*un_nxt);

		return abs((ut_nxt - ut_cur) / dr.norm());


	}

	wxLogMessage(_T("Error: velo deriv type not implemented, check TDomain::_calc_specifid_velo_deriv_abs"));

	return -1.0;

	
}

// TODO: works only for grids nearly orthogonal to viscous wall 
// input: data_grdline - for now just an extracted domain gridline
// output: bl_thick & raw_profile (truncated data_grdline)
void TDomain::_calc_bl_thick_vderiv(
		const std::vector<t_ZoneNode>& data_grdline, double& bl_thick, 
		std::vector<t_ZoneNode>& out_raw_profile) const{

	t_ZoneNode surf_znode, outer_znode;

	surf_znode = _get_nrst_node_surf(data_grdline[0]);

	t_GeomPoint wall_xyz, cur_xyz;

	t_Vec3Dbl surf_norm, dr;
	calc_surf_norm(surf_znode, surf_norm);

	wxLogMessage(_T("norm:%s"), surf_norm.to_wxstr());

	// wall rec is used to calculate bl_thick
	mf::t_Rec wall_rec, cur_rec;
	get_rec(surf_znode, wall_rec);
	wall_xyz.set(wall_rec);

	wxLogMessage(_T("surf node xyz:%s"), wall_xyz.to_wxstr());

	double du_dy = 0.0;
	double du_dy_max = 0.0;
	double du_dy_base = 0.0;

	double y_ref;

	t_ZoneNode ref_znode;

	t_VeloDerivType vd_type = t_VeloDerivType::VD_ABS;

	wxLogMessage(_T("_calc_bl_thick_vderiv derivative type: currently VD_ABS"));

	enum {VD_WALL, VD_MAX} vd_ref_type;

	vd_ref_type = VD_WALL;

	wxLogMessage(_T("_calc_bl_thick_vderiv derivative mode: currently WALL value is used"));

	int ind_bound = -1;

	// first find base reference deriv
	// then find reference point where dudy = eps*dudy_base : y = dd
	// then check that derivative is small inside dd and thick_coef*dd
	// then find outer record y_out = thick_coef*dd

	// step1 - find du_dy_base

	for (int m = 0; m < data_grdline.size() - 1; m++)
	{

		du_dy = _calc_specifid_velo_deriv_abs(data_grdline, m, vd_type);

		//wxLogMessage(_T("du_dy=%lf; dutau_dy=%lf"), du_dy, dutau_dy);

		if (vd_ref_type == VD_WALL) {

				if (m == 0) {

					// debug 
					wxLogMessage(_T("_calc_bl_thick_vderiv:: wall deriv = %lf"), du_dy);

					du_dy_base = du_dy;

					break;

				}

			}

		// find first max deriv value moving from the wall to outside
		if (vd_ref_type == VD_MAX) {

			if (du_dy < du_dy_max) {

				// debug 
				wxLogMessage(_T("_calc_bl_thick_vderiv:: max tang-velo deriv = %lf"), du_dy);

				du_dy_base = du_dy;

				break;
			}

		}
	} 

	// step 2 - find bl bound rec

	int ind_ref;

	for (int m = 0; m < data_grdline.size() - 1; m++) {

		du_dy = _calc_specifid_velo_deriv_abs(data_grdline, m, vd_type);

		double eps = _profile_cfg.DerivThreshold;

		if (du_dy<eps*du_dy_base) {

			get_rec(data_grdline[m], cur_rec);

			cur_xyz.set(cur_rec);

			matrix::base::minus<double, double>(cur_xyz, wall_xyz, dr);

			y_ref = dr.norm();

			// additional check required that from y_ref to 
			// y_ref*ThickCoef the derivative is below tolerance

			bool deriv_check = true;

			for (int j = m; j < data_grdline.size() - 1; j++) {

				get_rec(data_grdline[j], cur_rec);

				cur_xyz.set(cur_rec);

				matrix::base::minus<double, double>(cur_xyz, wall_xyz, dr);

				double cur_y = dr.norm();

				// TODO: ThickCoef should be here
				if (cur_y < _profile_cfg.ThickCoefDefault*y_ref) {

					double du_dy_up = _calc_specifid_velo_deriv_abs(data_grdline, j, vd_type);

					if (du_dy_up > eps*du_dy_base) {
						deriv_check = false;
						break;
					}

				}
				else { break; }

			}

			if (deriv_check) {
				ind_ref = m;
				bl_thick = y_ref;
				break;
			}

		}

	}

	// step 3 - extract cropped profile data
	for (int m = 0; m < data_grdline.size() - 1; m++) {

		get_rec(data_grdline[m], cur_rec);

		cur_xyz.set(cur_rec);

		matrix::base::minus<double, double>(cur_xyz, wall_xyz, dr);

		double cur_dd = dr.norm();

		// TODO: here should be ThickCoef, not default...
		if (cur_dd >= _profile_cfg.ThickCoefDefault*bl_thick) {

			outer_znode = data_grdline[m];

			out_raw_profile.resize(m + 1);

			for (int p = 0; p < m + 1; p++) out_raw_profile[p] = data_grdline[p];

			return;

		}

	}

	wxLogError(_T("Error: Zone boundary reached while searching BL bound & outer rec"));

};

// TODO : works only for grids nearly orthogonal to viscous wall
void TDomain::_calc_bl_thick_enthalpy(
	const std::vector<t_ZoneNode>& data_grdline, double& bl_thick, 
	std::vector<t_ZoneNode>& out_profile_data) const{

	t_ZoneNode surf_znode, outer_znode;

	surf_znode = data_grdline[0];

	outer_znode = data_grdline.back();

	mf::t_Rec cur_rec, nxt_rec;
	t_GeomPoint wall_xyz, cur_xyz, nxt_xyz, dr;
	double h_cur, h_inf;

	// wall rec is used to calculate bl_thick
	get_rec(surf_znode, cur_rec);
	wall_xyz.set(cur_rec);

	h_inf = calc_enthalpy_freestream();

	// march from outside of bl

	int m = data_grdline.size()-1;

	for (; m>=0; m--) 
	{

		get_rec(data_grdline[m], cur_rec);

		h_cur = calc_enthalpy(cur_rec);

		double tol = abs(h_cur/h_inf -1);

		if (tol>_profile_cfg.DerivThreshold){

			outer_znode = data_grdline[m];

			cur_xyz.set(cur_rec);

			matrix::base::minus<double, double>(cur_xyz, wall_xyz, dr);

			bl_thick = dr.norm();

			break;
		}

	};

	out_profile_data.resize(m+1);

	for (int p=0; p<m+1; p++) out_profile_data[p] = data_grdline[p];

};

void TDomain::_calc_bl_thick_full_gridline(
	const std::vector<t_ZoneNode>& data_grdline,
	double& bl_thick, std::vector<t_ZoneNode>& out_raw_profile) const{

	t_ZoneNode surf_znode, outer_znode;

	surf_znode = data_grdline[0];

	outer_znode = data_grdline.back();

	mf::t_Rec rec1, rec2;

	get_rec(surf_znode, rec1);
	get_rec(outer_znode, rec2);

	t_GeomPoint r1, r2;
	t_Vec3Dbl dr;

	r1.set(rec1);
	r2.set(rec2);

	matrix::base::minus<double, double>(r1, r2, dr);
	bl_thick = dr.norm();

	out_raw_profile = data_grdline;
	
}

inline void TDomain::_calc_bl_thick(const t_GeomPoint& xyz, double& bl_thick, 
									std::vector<t_ZoneNode>& raw_profile) const{

	std::vector<t_ZoneNode> data_grdline;

	_extract_profile_data_grdline(xyz, data_grdline);

	switch (_profile_cfg.BLThickCalcType)
	{
	case t_BLThickCalcType::BLTHICK_BY_VDERIV:
		_calc_bl_thick_vderiv(data_grdline, bl_thick, raw_profile);
		break;

	case t_BLThickCalcType::BLTHICK_BY_ENTHALPY:
		_calc_bl_thick_enthalpy(data_grdline, bl_thick, raw_profile);
		break;
	case t_BLThickCalcType::BLTHICK_FULL_GRIDLINE:
		_calc_bl_thick_full_gridline(data_grdline, bl_thick, raw_profile);
		break;

	default:
		wxString msg(_T("Provided BL Thick Calc Type not implemented"));
		wxLogError(msg); ssuGENTHROW(msg);
	}


}


double TDomain::calc_bl_thick(const t_GeomPoint& xyz) const{
	
	std::vector<t_ZoneNode> raw_profile;
	double bl_thick;

	_calc_bl_thick(xyz, bl_thick, raw_profile);
	return bl_thick;

};

void TDomain::calc_nearest_surf_rec(const t_GeomPoint& xyz, t_Rec& surf_rec) const{

	t_ZoneNode surf_znode;	

	surf_znode = _get_nrst_node_surf(xyz);

	t_BlkInd& ind = surf_znode.iNode;

	get_rec(Zones[surf_znode.iZone-1], ind.i, ind.j, ind.k, surf_rec);

};

void TDomain::calc_nearest_inviscid_rec(const t_GeomPoint& xyz, t_Rec& outer_rec) const{

	std::vector<t_ZoneNode> raw_profile;
	double bl_thick;

	_calc_bl_thick(xyz, bl_thick, raw_profile);

	get_rec(raw_profile.back(), outer_rec);

};

int TDomain::estim_num_bl_nodes(const t_GeomPoint& xyz) const{

	std::vector<t_ZoneNode> raw_profile;
	double bl_thick;

	_calc_bl_thick(xyz, bl_thick, raw_profile);

	return raw_profile.size();
}

t_SqMat3Dbl TDomain::calc_jac_to_loc_rf(const t_GeomPoint& xyz) const{
	
	std::vector<t_ZoneNode> raw_profile;
	double bl_thick;
	_calc_bl_thick(xyz, bl_thick, raw_profile);

	t_SqMat3Dbl jac;

	t_Rec bound_rec, surface_rec;
	get_rec(raw_profile.back(), bound_rec);
	get_rec(raw_profile[0], surface_rec);

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

	wxLogMessage(_T("bound rec xyz:%s"), bound_rec.get_xyz().to_wxstr());
	wxLogMessage(_T("bound rec UVW:%s"), bound_rec.get_uvw().to_wxstr());
	wxLogMessage(_T("Jac to loc rf:%s"), jac.to_wxstr());
	
	return jac;

}

// tmp func to print enthalpy profiles for enthalpy criteria
// IMPORTANT TODO: the profile coord Y is along gridline (may be curvilinear)
void TDomain::dump_full_enthalpy_profile(const mf::t_GeomPoint& xyz, int pid) const{

	t_ZoneNode surf_znode = _get_nrst_node_surf(xyz);

	t_DomainGrdLine dom_grdline(*this);

	dom_grdline.init(surf_znode);

	char fname[33];
	sprintf(fname, "%s/profile_enthalpy_%d.dat", hsstab::OUTPUT_DIR.ToAscii(), pid);
	std::ofstream ofstr(fname);

	mf::t_Rec mf_rec;
	mf::t_GeomPoint cur_xyz,prev_xyz, surf_xyz, dr;

	get_rec(dom_grdline.znodes[0], mf_rec);
	surf_xyz.set(mf_rec);

	// Prints all but 1 points 
	int M = dom_grdline.znodes.size();

	double s = 0.0;

	cur_xyz = surf_xyz; prev_xyz = surf_xyz;

	for (int m=0; m<M; m++){

		get_rec(dom_grdline.znodes[m], mf_rec);
		cur_xyz.set(mf_rec);
		matrix::base::minus<double, double>(cur_xyz, prev_xyz, dr);
		s+= dr.norm();
		double enth = calc_enthalpy(mf_rec);
		ofstr<<s<<"\t"<<enth<<"\n";

		prev_xyz = cur_xyz;
	}

};
