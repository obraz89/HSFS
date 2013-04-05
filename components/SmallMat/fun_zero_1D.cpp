#include "stdafx.h"
#include "fun_zero_1D.h"

using namespace smat;

void t_RootDicho::set_tol(double tol){_tol=tol;}

double t_RootDicho::search_root(double init_lft_arg, double init_rgt_arg){

	double lft = init_lft_arg;
	double rgt = init_rgt_arg;
	double mid, fun_l, fun_r, fun_mid;

	do
	{
		fun_l = _calc_fun(lft);
		fun_r = _calc_fun(rgt);
		mid = 0.5*(lft+rgt);
		fun_mid = _calc_fun(mid);
		if (fun_l*fun_mid<0.0){
			rgt = mid;
		}else{
			lft = mid;
		}
	}while(fun_l*fun_r<0.0&&abs(lft-rgt)>_tol);

	return mid;

}