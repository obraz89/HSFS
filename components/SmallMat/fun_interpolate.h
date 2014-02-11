#ifndef __FUN_INTERPOLATE
#define __FUN_INTERPOLATE

namespace smat{

	// instantiated in lib.cpp
	template<typename T> T interpolate_parab(double x1, T f1, double x2, T f2, 
		double x3, T f3, double x){

			return f1*(x-x2)*(x-x3)/((x1-x2)*(x1-x3))+
				f2*(x-x3)*(x-x1)/((x2-x1)*(x2-x3))+
				f3*(x-x1)*(x-x2)/((x3-x1)*(x3-x2));
	};

}

#endif  // __FUN_INTERPOLATE