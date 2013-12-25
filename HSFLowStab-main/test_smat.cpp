#include "stdafx.h"

#include "math_operands.h"
#include "tests.h"

#include "time.h"
#include <sys/types.h>
#include <sys/timeb.h>


void test::test_smat(){
	/*
	t_MatCmplx m1(3,4), m2(4,3);
	m1 = m1*m2;
	t_VecCmplx v(3);
	m1 = v*m1;
	//	m1 = m2 + v;
	t_Complex norm = v.norm();
	v = v-v;
	m1[0]=v;
	*/

	/*
	int N = 1.0E+08;
	_timeb time_start_b, time_end_b;
	double time_spent;

	//std::vector<double> vec1(3, 1.), vec2(3, 2.), vec3(3, 3.);

	t_Vec3Dbl vec1, vec2, vec3;

	for (int i=0;i<3; i++) {
		vec1[i] = vec2[i] = vec3[i] = 1.;
	}
	//t_Vec3Dbl vec1, vec2, vec3;
	_ftime_s(&time_start_b);

	for (int n=0; n<N; n++){

		//double *vec1 = new double[3];
		//double *vec2 = new double[3];
		//double *vec3 = new double[3];

		//double **new_mat = new double*[3];

		//double new_mat[3][3];

		//for (int i=0; i<3; i++){
			
			//new_mat[i] = new double[3];

		//	mat1[0][i] = 1.;
		//    vec1[i] = vec2[i] = vec3[i] = 1.;
		//	vec3[i] = vec1[i] + vec2[i];
		//}
		
		matrix::base::plus<double, double>(vec1, vec2, vec3);
		//vec3 = vec1 + vec2;

	}

	_ftime_s(&time_end_b);

	time_spent = (time_end_b.time - time_start_b.time) + ((double)(time_end_b.millitm - time_start_b.millitm))/1000.;
	std::cout<<"time spent:"<<time_spent;
	*/

	t_Vec3Dbl r;
	t_SqMat3Dbl l;
	t_Vec3Dbl ret;
	r[0] = 0.;
	r[1] = 1.;
	r[2] = 0.;
	l.setToUnity();
	ret = l*r;
	std::wcout<<ret;
	matrix::base::mat_mul<double, double>(l,r,ret);
	std::wcout<<ret;
	getchar();
};