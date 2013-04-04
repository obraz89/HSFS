#include "stdafx.h"

#include "io_helpers.h"

#include "conj_minmax.h"

#include "tests.h"

using namespace smat;

double test_fun1(double x, double y){return 24.0-5.0*(x-1.0)*(x-1.0)-(y-1.0)*(y-1.0);}

t_VecDbl test_fun1_grad(double x, double y){

	t_VecDbl ret(2);
	ret[0]=-10.0*(x-1);
	ret[1]=-2.0*(y-1);
	return ret;
}

void test::test_conj_grad_min2D(){

	t_Conj2D conj_srch(test_fun1, test_fun1_grad);

	t_VecDbl init_guess(2);
	init_guess[0]=7.0;
	init_guess[1]=3.0;

	conj_srch.search_max(init_guess);
	int bla=1;
}

