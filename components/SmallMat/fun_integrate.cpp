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