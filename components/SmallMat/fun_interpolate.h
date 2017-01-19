#ifndef __FUN_INTERPOLATE
#define __FUN_INTERPOLATE

#include "dll_impexp_smat.h"

namespace smat{

	// instantiated in lib.cpp
	template<typename T> T interpolate_parab(double x1, T f1, double x2, T f2, 
		double x3, T f3, double x){

			return f1*(x-x2)*(x-x3)/((x1-x2)*(x1-x3))+
				f2*(x-x3)*(x-x1)/((x2-x1)*(x2-x3))+
				f3*(x-x1)*(x-x2)/((x3-x1)*(x3-x2));
	};

	template<typename T> T interpolate_parab(double x_v[3], T f_v[3], double x){
		return interpolate_parab<T>(x_v[0], f_v[0], x_v[1], f_v[1], x_v[2], f_v[2], x);
	}

	IMPEXP_SMALLMAT void interpolate_profile_sm_deriv_cubic(double* x, double* y, const int nnodes, double* y1, double* y2);

}

#endif  // __FUN_INTERPOLATE