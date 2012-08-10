#include "math_operands.h"

// tests
void abc(){
	t_MatCmplx m1(3,4), m2(4,3);
	m1 = m1*m2;
	t_VecCmplx v(3);
	m1 = v*m1;
//	m1 = m2 + v;
	t_Complex norm = v.norm();
	v = v-v;
	m1[0]=v;
};


