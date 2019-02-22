#include "stdafx.h"
#include "mf_shared.h"

#include "io_helpers.h"

using namespace mf;

//-------------------------------------------------------------------------t_Rec

void t_Rec::set_xyz(const t_GeomPoint& point){
	x = point.x();
	y = point.y();
	z = point.z();
};

void t_Rec::set_uvw(const t_Vec3Dbl& vec){
	u = vec[0];
	v = vec[1];
	w = vec[2];
}

t_GeomPoint t_Rec::get_xyz() const{return t_GeomPoint(x,y,z);};
t_Vec3Dbl t_Rec::get_uvw() const{return t_Vec3Dbl(u,v,w);};

std::wostream& operator<<(std::wostream& os, const t_Rec& rec){
	os<<_T("x:")<<std_manip::std_format_fixed<double>(rec.x)<<
		_T("y:")<<std_manip::std_format_fixed<double>(rec.y)<<
		_T("z:")<<std_manip::std_format_fixed<double>(rec.z)<<std::endl
		<<_T("u:")<<std_manip::std_format_fixed<double>(rec.u)<<
		_T("v:")<<std_manip::std_format_fixed<double>(rec.v)<<
		_T("w:")<<std_manip::std_format_fixed<double>(rec.w)<<std::endl
		<<_T("p:")<<std_manip::std_format_fixed<double>(rec.p)<<
		_T("t:")<<std_manip::std_format_fixed<double>(rec.t)<<
		_T("r:")<<std_manip::std_format_fixed<double>(rec.r)<<std::endl;
	return os;
};

//---------------------------------------------------------------------operators

t_Rec t_Rec::lin_comb(double c1, const t_Rec& rec1, double c2, const t_Rec& rec2) {

	t_Rec res;
	res.x = c1*rec1.x + c2*rec2.x;
	res.y = c1*rec1.y + c2*rec2.y;
	res.z = c1*rec1.z + c2*rec2.z;
	res.u = c1*rec1.u + c2*rec2.u;
	res.v = c1*rec1.v + c2*rec2.v;
	res.w = c1*rec1.w + c2*rec2.w;
	res.p = c1*rec1.p + c2*rec2.p;
	res.t = c1*rec1.t + c2*rec2.t;
	res.r = c1*rec1.r + c2*rec2.r;
	return res;

}

t_Rec mf::operator-(const t_Rec& rec1, const t_Rec& rec2){
	t_Rec res;
	res.x = rec1.x - rec2.x;
	res.y = rec1.y - rec2.y;
	res.z = rec1.z - rec2.z;
	res.u = rec1.u - rec2.u;
	res.v = rec1.v - rec2.v;
	res.w = rec1.w - rec2.w;
	res.p = rec1.p - rec2.p;
	res.t = rec1.t - rec2.t;
	res.r = rec1.r - rec2.r;
	return res;
};