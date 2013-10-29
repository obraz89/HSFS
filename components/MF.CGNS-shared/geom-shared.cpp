#include "stdafx.h"
#include "cgns_structs.h"

static const double BL_BOUND_VELO_TOL = 0.01;

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

		for (int iface = 0; iface<4; iface++){
			if (_is_face_of_bcwall_type(blk.Faces[iface].szName)){
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

void TDomain::_calc_bl_thick(const t_GeomPoint& xyz, double& bl_thick, 
					t_ZoneNode& surf_znode, t_ZoneNode& outer_znode) const{

	surf_znode = _get_nrst_node_surf(xyz);

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

// get du_dn wall
	mf::t_Rec cur_rec, nxt_rec;
	t_GeomPoint cur_xyz, nxt_xyz, dr;
	t_GeomPoint wall_xyz, out_xyz;
	t_Vec3Dbl cur_uvw, nxt_uvw, du;
	double du_dn_wall, du_dn_cur;

	get_rec(blk, i_s, j_s, k_s, cur_rec);
	get_rec(blk, i_s + di, j_s + dj, k_s + dk, nxt_rec);

	cur_xyz.set(cur_rec);
	nxt_xyz.set(nxt_rec);

	cur_uvw.set(cur_rec.u, cur_rec.v, cur_rec.w);
	nxt_uvw.set(nxt_rec.u, nxt_rec.v, nxt_rec.w);

	matrix::base::minus<double, double>(nxt_xyz, cur_xyz, dr);	
	matrix::base::minus<double, double>(nxt_uvw, cur_uvw, du);

	du_dn_wall = du.norm()/dr.norm();
	wall_xyz.set(cur_rec);

	// 1 step done
	i+=di;
	j+=dj;
	k+=dk;
	
	do 
	{
		do 
		{
			do 
			{
				
				get_rec(blk, i, j, k, cur_rec);
				get_rec(blk, i + di, j + dj, k + dk, nxt_rec);

				cur_xyz.set(cur_rec);
				nxt_xyz.set(nxt_rec);

				cur_uvw.set(cur_rec.u, cur_rec.v, cur_rec.w);
				nxt_uvw.set(nxt_rec.u, nxt_rec.v, nxt_rec.w);

				matrix::base::minus<double, double>(nxt_xyz, cur_xyz, dr);	
				matrix::base::minus<double, double>(nxt_uvw, cur_uvw, du);

				du_dn_cur = du.norm()/dr.norm();

				if (du_dn_cur<=du_dn_wall*BL_BOUND_VELO_TOL)
				{

					outer_znode.iZone = surf_znode.iZone;
					outer_znode.iNode.set(i,j,k);

					out_xyz = cur_xyz;
					matrix::base::minus<double, double>(cur_xyz, wall_xyz, dr);
					bl_thick = dr.norm();
					return;
				}
				

				k+=dk;
			} while ((k_e-k)*dk>0);
			j+=dj;
		} while ((j_e-j)*dj>0);
		i+=di;
	} while ((i_e-i)*di>0);

	wxLogError(_T("Failed to calculate bl thickness, iZone=%d"), surf_znode.iZone);
	return;
	

};


double TDomain::calc_bl_thick(const t_GeomPoint& xyz) const{
	
	t_ZoneNode surf_znode, outer_znode;
	double bl_thick;

	_calc_bl_thick(xyz,bl_thick, surf_znode, outer_znode);
	return bl_thick;

};

int TDomain::estim_num_bl_nodes(t_GeomPoint xyz) const{

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
