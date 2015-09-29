#include "stdafx.h"
#include "fun_integrate.h"
#include "gen_exception.h"

IMPEXP_SMALLMAT double smat::fun_integrate(
	const std::vector<double>& x, const std::vector<double>& y, int a_nsteps/* =-1 */){

	if (x.size()!=y.size() || a_nsteps>=x.size()){
		wxLogError(_T("smat: Fun Integrate Error")); return 0.0;
	}

	int n = a_nsteps>0 ? a_nsteps : x.size()-1;

	double ret=0.0;

	for (int i=0; i<n; i++){

		ret+=0.5*(y[i]+y[i+1])*(x[i+1]-x[i]);

	}

	return ret;

}

IMPEXP_SMALLMAT t_Complex smat::fun_integrate(
	const std::vector<double>& x, const std::vector<t_Complex>& y, int a_nsteps/* =-1 */){

		int i1 = x.size();
		int i2 = y.size();
		if ((i1!=i2)||(a_nsteps>=i1)){
			wxLogError(_T("smat: Fun Integrate Error")); return 0.0;
		}

		int n = a_nsteps>0 ? a_nsteps : i1-1;

		t_Complex ret=0.0;

		for (int i=0; i<n; i++){

			ret+=0.5*(y[i]+y[i+1])*(x[i+1]-x[i]);

		}

		return ret;

}

template<typename T> inline T _integrate_simp4_uniform(const std::vector<double>& x, const std::vector<T>& y){
	if (x.size()!=y.size()) wxLogError(_T("Error in integration - size mismatch"));

	int n;

	if (x.size() % 2 ==0) {
		wxLogError(_T("Smat: Fun Integrate: Number of points N should be odd, using N-1"));
		n = x.size()-1;
	}else{
		n = x.size();
	}

	double d = x[1] - x[0];

	T ret = d/3.0*(y[0]+y[n-1]);

	for (int i=1; i<n-1; i+=2) ret+=4.0/3.0*d*y[i];
	for (int i=2; i<n-2; i+=2) ret+=2.0/3.0*d*y[i];

	return ret;
};

IMPEXP_SMALLMAT double smat::fun_integrate_simp4_uniform(const std::vector<double>& x, const std::vector<double>& y, int nsteps){
	//return 0.0;
	return _integrate_simp4_uniform<double>(x,y);
}

IMPEXP_SMALLMAT t_Complex smat::fun_integrate_simp4_uniform(const std::vector<double>& x, const std::vector<t_Complex>& y, int nsteps){
	//return 0.0;
	return _integrate_simp4_uniform<t_Complex>(x,y);
}

// find ff - the antiderivative of ff over x
// note that ff is calculated on a different grid - in halfnodes i+1/2
// the overbound value 0-1/2 is set to zero
IMPEXP_SMALLMAT void smat::integrate_over_range(
	const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& ff){

	ff[0] = 0;

	for (int i=1; i<x.size(); i++){
		ff[i] = smat::fun_integrate(x, y, i);
	}

}