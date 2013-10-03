#include "stdafx.h"
#include "MFHSDomain.h"

using namespace mf;
using namespace mfhs;

//---------------------------------------------------------------------t_Domain

void t_Domain::calc_surf_point(const mf::t_GeomPoint& a_xyz, mf::t_GeomPoint& surf_point, t_Vec3Dbl& norm) const{
	return get_blk().calc_surf_point(a_xyz, surf_point, norm);
}

t_Rec t_Domain::get_rec(const t_GeomPoint& xyz) const{

	return get_blk().interpolate_to_point(xyz);

};

t_Rec t_Domain::interpolate_to_point(const t_GeomPoint& xyz) const{

	return get_blk().interpolate_to_point(xyz);
};

t_SqMat3Dbl t_Domain::calc_jac_to_loc_rf(const t_GeomPoint& xyz) const{

	const t_Block& blk = get_blk();
	// TODO: get nearest index exact!
	t_BlkInd ind = blk.get_nearest_index_raw(xyz);

	return blk.calc_jac_to_loc_rf(ind);
};

double t_Domain::calc_enthalpy(const t_GeomPoint& xyz) const{

	const t_Block& blk = get_blk();
	// TODO: get nearest index exact!
	t_BlkInd ind = blk.get_nearest_index_raw(xyz);

	return blk.calc_enthalpy(ind);

};

double t_Domain::calc_enthalpy_freestream() const{

	return get_blk().calc_enthalpy_freestream();
};

double t_Domain::calc_viscosity(const double t) const{
	return get_blk().calc_viscosity(t);
}

double t_Domain::calc_viscosity(const t_GeomPoint& xyz) const{

	const t_Block& blk = get_blk();
	// TODO: get nearest index exact!
	t_BlkInd ind = blk.get_nearest_index_raw(xyz);

	return blk.calc_viscosity(ind);

};

double t_Domain::calc_cin_visc_dim(const t_GeomPoint& xyz)  const{

	const t_Block& blk = get_blk();
	// TODO: get nearest index exact!
	t_BlkInd ind = blk.get_nearest_index_raw(xyz);

	return get_blk().calc_cin_visc_dim(ind);
};

double t_Domain::calc_cin_visc_inf_dim() const{
	return get_blk().calc_cin_visc_inf_dim();
};

double t_Domain::calc_u_inf() const{
	return get_blk().calc_u_inf();
};

double t_Domain::calc_c_dim(const double t) const{

	return get_blk().calc_c_dim(t);

}

double t_Domain::calc_c_dim(const t_GeomPoint& xyz) const{

	const t_Block& blk = get_blk();
	// TODO: get nearest index exact!
	t_BlkInd ind = blk.get_nearest_index_raw(xyz);

	return blk.calc_c_dim(ind);
};


double t_Domain::calc_mach(const t_GeomPoint& xyz) const{

	const t_Block& blk = get_blk();
	// TODO: get nearest index exact!
	t_BlkInd ind = blk.get_nearest_index_raw(xyz);

	return blk.calc_mach(ind);

};

double t_Domain::calc_bl_thick(const t_GeomPoint& xyz) const{
	return get_blk().calc_bl_thick(xyz);
}

double t_Domain::calc_x_scale(const mf::t_GeomPoint& xyz) const{
	return get_blk().calc_x_scale(xyz);
}

bool t_Domain::is_point_inside(const t_GeomPoint& xyz) const{
	return get_blk().is_point_inside(xyz);
}

int  t_Domain::estim_num_bl_nodes(mf::t_GeomPoint xyz) const{
	const t_Block& blk = get_blk();
	t_BlkInd nrst_ind = blk.get_nearest_index_raw(xyz);
	return blk.get_bound_index(nrst_ind);
};

// tmp