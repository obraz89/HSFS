#include "WavePackLine.h"

#include "log.h"

using namespace mf;
using namespace pf;
using namespace stab;

// calculate derivative of amplitude function along wpline
// (usually in 2d task it is just deriv along x)
// i - index of record in wpline
// IMPORTANT TODO: this is slightly incorrect way to compute deriv
// because the physical grids of adjacent profiles are different in normal direction (correct this !)
// IMPORTANT: here x must be non-dim by stab scale
void t_WavePackLine::_calc_amp_fun_deriv_dx(int i, stab::t_LSBase& loc_solver, std::vector<t_VecCmplx>& fun_l,
	std::vector<t_VecCmplx>& fun_r, std::vector<t_VecCmplx>& amp_funcs_deriv) {

	// TODO: correct expressions for first and last points
	// for now using second-order approx from second and last but one derivs

	if (i == 0) i = 1;
	if (i == _line.size() - 1) i = _line.size() - 2;

	t_GeomPoint p1, p2;

	p1 = _line[i - 1].mean_flow.get_xyz();
	p2 = _line[i + 1].mean_flow.get_xyz();

	double coef = _line[i].wchars_loc.scales().Dels/_rFldMF.get_mf_params().L_ref;
	//t_StabScales scales = loc
	double dx = (p2 - p1).norm()/coef;

	get_amp_funcs(i - 1, loc_solver, fun_l);
	get_amp_funcs(i + 1, loc_solver, fun_r);

	loc_solver.normalizeAmpFuncsByPressureAtWall(fun_l);
	loc_solver.normalizeAmpFuncsByPressureAtWall(fun_r);

	int nnodes_stab = fun_l.size();

	const int stab_matrix_dim = 8;

	for (int j = 0; j < nnodes_stab; j++) {
		for (int i = 0; i < stab_matrix_dim; i++)
			amp_funcs_deriv[j][i] = 0.5*(fun_r[j][i] - fun_l[j][i]) / dx;
	}
	

}
// weakly non-parallel addition to eigenvalue
// see Fedorov AIAA-2002-2846, Nayfeh AIAA-J ARTICLE NO. 79-0262R
// solution in non-par approximation:
// Z = C(x)*dze(y;x)*exp(iS)
// dC/dx = WC
// W = -(<H1*dze_dx, ksi> + <H2*dze, ksi>)/<H1*dze, ksi>
// C = C_0*exp(Wx) => da = -i*W
// -i*W is thus addition to eigenvalue a
// NB : dze_dx is computed with amp func normalization, so dze must be normalized here too!
// NB : ksi has arbitrary normalization
void t_WavePackLine::_calc_nonpar_sigma_additions(stab::t_LSBase& loc_solver) {

	std::vector<t_VecCmplx> fun_l(loc_solver.getNNodes(), t_VecCmplx(8));
	std::vector<t_VecCmplx> fun_r(loc_solver.getNNodes(), t_VecCmplx(8));
	std::vector<t_VecCmplx> dze_ddx(loc_solver.getNNodes(), t_VecCmplx(8));
	// direct amplitude func
	std::vector<t_VecCmplx> dze(loc_solver.getNNodes(), t_VecCmplx(8));
	// conjugate amplitude func 
	std::vector<t_VecCmplx> ksi(loc_solver.getNNodes(), t_VecCmplx(8));

	t_GeomPoint xyz;

	t_WCharsLoc wchars;

	t_Complex v1, v2, v3;

	t_Complex W;

	const t_Complex E = t_Complex(0.0, 1.0);

	t_Complex da_ratio;

	for (int i = 0; i < _line.size(); i++) {

		xyz = _line[i].mean_flow.get_xyz();

		loc_solver.setContext(xyz);

		wchars = _line[i].wchars_loc;

		// get conjugate amp fun
		loc_solver.setLSMode(stab::t_LSMode(stab::t_LSMode::CONJUGATE | stab::t_LSMode::ASYM_HOMOGEN));

		loc_solver.searchWave(wchars, stab::t_LSCond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED),
			stab::t_TaskTreat::SPAT);

		loc_solver.getAmpFuncs(ksi);

		// get direct amp fun

		loc_solver.setLSMode(stab::t_LSMode(stab::t_LSMode::DIRECT | stab::t_LSMode::ASYM_HOMOGEN));

		wchars = _line[i].wchars_loc;
		loc_solver.searchWave(wchars, stab::t_LSCond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED),
			stab::t_TaskTreat::SPAT);

		loc_solver.getAmpFuncs(dze);
		loc_solver.normalizeAmpFuncsByPressureAtWall(dze);

		// get deriv of amp_fun along x
		_calc_amp_fun_deriv_dx(i, loc_solver, fun_l, fun_r, dze_ddx);

		// restore context at this point after d_dx computations
		loc_solver.setContext(xyz);

		wchars = _line[i].wchars_loc;

		loc_solver.searchWave(wchars, stab::t_LSCond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED),
			stab::t_TaskTreat::SPAT);

		// calc <H1*dze_ddx, ksi>
		v1 = loc_solver.calcScalarProd_H1(dze_ddx, ksi);

		// calc <H2*dze, ksi>
		v2 = 0.0;// loc_solver.calcScalarProd_H2(dze, ksi);

		// calc <H1*dze, ksi>
		v3 = loc_solver.calcScalarProd_H1(dze, ksi);

		wxLogMessage(_T("v1=(%lf, %lf)"), v1.real(), v1.imag());
		wxLogMessage(_T("v2=(%lf, %lf)"), v2.real(), v2.imag());
		wxLogMessage(_T("v3=(%lf, %lf)"), v3.real(), v3.imag());

		W = -(v1 + v2) / v3;

		_line[i].da_nonpar = -1.0*E*W;

		// debug 
		wxLogMessage(_T("da_nonpar = (%lf, %lf)"), 
			_line[i].da_nonpar.real(), _line[i].da_nonpar.imag());

		da_ratio = _line[i].da_nonpar / wchars.a;

		wxLogMessage(_T("da_ratio = da/a = (%lf, %lf)"), 
			da_ratio.real(), da_ratio.imag());

		// modify wave chars
		_line[i].wchars_loc.a += _line[i].da_nonpar;

		t_WCharsGlob wchars_glob(_line[i].wchars_loc, _rFldMF.calc_jac_to_loc_rf(xyz),
			loc_solver.get_stab_scales());

		_line[i].wave_chars = wchars_glob;	

	}

}

