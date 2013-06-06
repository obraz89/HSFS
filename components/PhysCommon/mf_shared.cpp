#include "stdafx.h"
#include "mf_shared.h"

using namespace mf;

t_GeomPoint::t_GeomPoint(double x/*=0*/, double y/*=0*/, double z/*=0*/)
	:t_Vec3Dbl(x,y,z){};

t_GeomPoint::t_GeomPoint(const t_Vec3Dbl& v):t_Vec3Dbl(v){};

double t_GeomPoint::x() const{return this->operator[](0);};
double& t_GeomPoint::x(){return this->operator[](0);};

double t_GeomPoint::y() const{return this->operator[](1);};
double& t_GeomPoint::y(){return this->operator[](1);};

double t_GeomPoint::z() const{return this->operator[](2);};
double& t_GeomPoint::z(){return this->operator[](2);}; 

//------------------------------------------------------------- Field parameters

const int t_ViscType::ViscPower=0;
const int t_ViscType::ViscSuther=1;