#include "stdafx.h"

#include "MF_CGNS2D.h"
#include "common_data.h"

using namespace hsstab;
using namespace mf;
using namespace mf::cg;
//

void t_MFCGNS2D::init(const TPlugin& g_plug){

	const TPluginParamsGroup& g = g_plug.get_settings_grp_const("");

	_fld_bin_path = g.get_string_param("FldBinPath");
	_grd_bin_path = g.get_string_param("GrdBinPath");

	if (!wxFileName::FileExists(_fld_bin_path) || 
		!wxFileName::FileExists(_grd_bin_path)){
			wxLogError(_T("Fatal Error: MF.CGNS2D: Provided fld or grd file doesn't exist"));
			ssuTHROW(t_GenException, _T("Error: Failed to initialize MF.CGNS2D, see errlog"));
	}

	mf::cg::hsf2d::_init_fld_base_params(_base_params, g);

	_base_params.MFSym = g.get_int_param("AxeSym")==mf::t_AxeSym::AxeSym ? 
		mf::t_AxeSym::AxeSym : mf::t_AxeSym::Plane;

	_base_params.ZSpan = g.get_real_param("ZSpan");

	_base_params.ThetaSpan = g.get_real_param("ThetaSpan");

	_base_params.Nz = g.get_int_param("Nz");

	if (_base_params.Nz<1){
		wxLogError(_T("Error: Bad Nz value during CGNS2D fld initialization, check ini files"));
		ssuTHROW(t_GenException, _T("Error: in CGNS2D init, see err log"));
	}

	// TODO: read from config ?
	nu = 4;  // unknown functions number
	nDim = 2;

	//G_strFunctionNames = "u\nv\np\nT";

	G_vecCGNSFuncNames.clear();
	G_vecCGNSFuncNames.push_back("VelocityX");
	G_vecCGNSFuncNames.push_back("VelocityY");
	G_vecCGNSFuncNames.push_back("Pressure");
	G_vecCGNSFuncNames.push_back("Temperature");

	// here come all viscous wall bc identifiers
	// TODO: do i need to keep this 33 size to compare in _is_face_of_bcwall_type
	// TODO: make entry in config file for wall BCs. With Regexps. o_O

	char viscBCWallName[33];

	sprintf(viscBCWallName, "blk1-Ymin");
	_vecBCWallNames.push_back(std::string(viscBCWallName));

	sprintf(viscBCWallName, "blk2-Ymin");
	_vecBCWallNames.push_back(std::string(viscBCWallName));

	sprintf(viscBCWallName, "blk3-Ymin");
	_vecBCWallNames.push_back(std::string(viscBCWallName));

	sprintf(viscBCWallName, "blk4-Ymin");
	_vecBCWallNames.push_back(std::string(viscBCWallName));

	_init();	// allocate space and read grd and fld

}

void t_MFCGNS2D::_init(){

	loadGrid(_grd_bin_path);
	initField(_fld_bin_path);

}

t_Rec t_MFCGNS2D::get_rec(const t_GeomPoint& xyz) const{

	return interpolate_to_point(xyz);

};

void t_MFCGNS2D::get_rec(const TZone& blk, int i, int j, int k, mf::t_Rec& rec) const{

// fld
	// in fact we have a 2D field
	int k_fld = 1;

	int glob_ind = blk.absIdx(i,j,k_fld);
	// cgns-shared.cpp ln.167
	int glob_offset = nu*( blk.absIdx(i,j,k_fld) - 1 );
	// TODO: avoid hardcoding?
	double u_base = blk.U[glob_offset + 0];
	double v_base = blk.U[glob_offset + 1];
	double p_base = blk.U[glob_offset + 2];
	double t_base = blk.U[glob_offset + 3];

	const t_CGNS2DParams& params = _base_params;
	double gmama = params.Gamma*params.Mach*params.Mach;
	double r_base = gmama*p_base/t_base;

	const TgridCell2D& cell = blk.cell(i, j);
	double x_base = cell.ij.x;
	double y_base = cell.ij.y;

	if (params.MFSym==mf::t_AxeSym::AxeSym){

		double psi = (-0.5+double(k-1)/double(params.Nz-1))*params.ThetaSpan;
		//double psi = M_PI/double(params.Nz-1)*(k-1);
		rec.x = x_base;
		rec.y = y_base*cos(psi);
		rec.z = y_base*sin(psi);

		rec.u = u_base;
		rec.v = v_base*cos(psi);
		rec.w = v_base*sin(psi);
		rec.p = p_base;
		rec.t = t_base;
		rec.r = r_base;
	}else{
		rec.x = x_base;
		rec.y = y_base;
		rec.z = (-0.5+double(k-1)/double(params.Nz-1))*params.ZSpan; // z from -0.5*ZSpan to 0.5*ZSpan;

		rec.u = u_base;
		rec.v = v_base;
		rec.w = 0.0;
		rec.p = p_base;
		rec.t = t_base;
		rec.r = r_base;

	};

}

//t_ZoneNode t_MFCGNS2D


t_Rec t_MFCGNS2D::interpolate_to_point(const t_GeomPoint& point) const{

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

#define SET_RET(o) ret.o## = a*rec1.o## + b*rec2.o##;
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
				_T("Failed to interpolate to point. Zoneid = %d"), 
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


	}  // ~interpolate_to_point


