#include "stdafx.h"
#include "MFBlockBase.h"

using namespace mf;

//-------------------------------------------------------------------------t_Rec

void t_Rec::set_xyz(t_GeomPoint point){
	x = point.x();
	y = point.y();
	z = point.z();
};

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

//----------------------------------------------------------------------t_BlkInd
t_BlkInd::t_BlkInd():i(0),j(0),k(0){};
t_BlkInd::t_BlkInd(int _i, int _j, int _k):i(_i), j(_j), k(_k){};
t_BlkInd::t_BlkInd(const t_BlkInd& _ind, int di, int dj, int dk)
{
	i = _ind.i + di;
	j = _ind.j + dj;
	k = _ind.k + dk;
};

bool mf::operator==(const t_BlkInd a, const t_BlkInd b)
{
	return ((a.i==b.i)&&(a.j==b.j)&&(a.k==b.k));
};
bool mf::operator!=(const t_BlkInd a, const t_BlkInd b)
{
	return !(mf::operator==(a,b));
}

std::wostream& mf::operator<<(std::wostream& str, const t_BlkInd& ind){
	return str<<_T("[")
		<<ind.i<<_T(";")
		<<ind.j<<_T(";")
		<<ind.k<<_T("]");
};

//---------------------------------------------------------------------operators

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
	return res;
};