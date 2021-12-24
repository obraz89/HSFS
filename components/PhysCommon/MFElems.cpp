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

double t_Rec::get_val(char name) {
	if (name == 'x')
		return x;
	if (name == 'y')
		return y;
	if (name == 'z')
		return z;
	if (name == 'u')
		return u;
	if (name == 'v')
		return v;
	if (name == 'w')
		return w;
	if (name == 'p')
		return p;
	if (name == 't')
		return t;
	if (name == 'r')
		return r;

	wxLogError(_T("mf::t_Rec::get_val: unknown name"));
	return 0.0;
}

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

// RecGradKed

t_Vec3Dbl& t_RecGradKed::get_vec(char name) {

	if (name == 'x')
		return xked;
	if (name == 'y')
		return yked;
	if (name == 'z')
		return zked;
	if (name == 'u')
		return uked;
	if (name == 'v')
		return vked;
	if (name == 'w')
		return wked;
	if (name == 'p')
		return pked;
	if (name == 't')
		return tked;
	if (name == 'r')
		return rked;

	wxLogError(_T("mf::t_RecGradKed::get_val: unknown name"));
	return xked;

}

t_RecGrad t_RecGrad::lin_comb(double c1, const t_RecGrad& r1, double c2, const t_RecGrad& r2) {
	t_RecGrad ret;

	ret.ug = c1*r1.ug + c2*r2.ug;
	ret.vg = c1*r1.vg + c2*r2.vg;
	ret.wg = c1*r1.wg + c2*r2.wg;
	ret.pg = c1*r1.pg + c2*r2.pg;
	ret.tg = c1*r1.tg + c2*r2.tg;

	return ret;
}