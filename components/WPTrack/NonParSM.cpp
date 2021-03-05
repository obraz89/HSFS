#include "WavePackLine.h"

#include "log.h"

using namespace mf;
using namespace pf;
using namespace stab;

// calculate derivative of amplitude function along wpline
// (usually in 2d task it is just deriv along x)
// i - index of record in wpline
// IMPORTANT: grids of y distribution in physical space must exactly coincide
// for points in the stencil (x+dx and x-dx)
// IMPORTANT: here x must be non-dim by stab scale
void t_WavePackLine::_calc_amp_fun_deriv_dx(int i, stab::t_LSBase& loc_solver, std::vector<t_VecCmplx>& fun_l,
	std::vector<t_VecCmplx>& fun_r, std::vector<t_VecCmplx>& amp_funcs_deriv) {

	// TODO: correct expressions for first and last points
	// for now using second-order approx from second and last but one derivs

	if (i == 0) i = 1;
	if (i == _line.size() - 1) i = _line.size() - 2;

	t_GeomPoint pp, p1, p2;

	pp = _line[i].mean_flow.get_xyz();
	p1 = _line[i - 1].mean_flow.get_xyz();
	p2 = _line[i + 1].mean_flow.get_xyz();

	double coef = _line[i].wchars_loc.scales().Dels/_rFldMF.get_mf_params().L_ref;
	//t_StabScales scales = loc
	const double dx = (p2 - p1).norm()/coef;
	const double c = 0.5 / dx;

	// set context to get thickness to be used as thick for +- calcs
	loc_solver.setContext(pp);

	double y_thick = loc_solver.getThickMF();

	mf::t_ProfDataCfg prof_cfg;

	prof_cfg.ThickFixed = y_thick;

	// get amp funcs with fixed thickness of prof stab
	get_amp_funcs(i - 1, loc_solver, fun_l, &prof_cfg);
	get_amp_funcs(i + 1, loc_solver, fun_r, &prof_cfg);

	loc_solver.normalizeAmpFuncsByPressureAtWall(fun_l);
	loc_solver.normalizeAmpFuncsByPressureAtWall(fun_r);

	int nnodes_stab = fun_l.size();

	const int stab_matrix_dim = 8;

	for (int j = 0; j < nnodes_stab; j++) {
		for (int i = 0; i < stab_matrix_dim; i++)
			amp_funcs_deriv[j][i] = c*(fun_r[j][i] - fun_l[j][i]);
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

	t_Complex v1_a, v2_a, v3_a;
	t_Complex v1_w, v2_w, v3_w;

	t_Complex W;

	const t_Complex E = t_Complex(0.0, 1.0);

	t_Complex da_ratio;

	// tmp, debug group velo
	std::ofstream ofstr("output/wplines_group_velo.dat");
	ofstr << "da_dw_fd_real, da_dw_fd_imag, da_dw_mat_real, da_dw_mat_imag\n";

	for (int i = 0; i < _line.size(); i++) {

		xyz = _line[i].mean_flow.get_xyz();

		loc_solver.setContext(xyz);

		wchars = _line[i].wchars_loc;

		// get conjugate amp fun
		loc_solver.setLSMode(stab::t_LSMode(stab::t_LSMode::CONJUGATE | stab::t_LSMode::ASYM_HOMOGEN));

		loc_solver.searchWave(wchars, stab::t_LSCond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED),
			stab::t_TaskTreat::SPAT);

		loc_solver.getAmpFuncs(ksi);

		// x=0.923
		if (i==157)
			stab::dumpEigenFuncs("output/amp_funcs_conjug.dat", loc_solver.getNNodes(), loc_solver.get_y_distrib(), ksi);
			//loc_solver.dumpEigenFuctions("output/amp_funcs_conjug.dat");

		// get direct amp fun

		loc_solver.setLSMode(stab::t_LSMode(stab::t_LSMode::DIRECT | stab::t_LSMode::ASYM_HOMOGEN));

		wchars = _line[i].wchars_loc;
		loc_solver.searchWave(wchars, stab::t_LSCond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED),
			stab::t_TaskTreat::SPAT);

		loc_solver.getAmpFuncs(dze);

		// x=0.923
		if (i == 157) {
			wxLogMessage(_T("Current point x=%lf"), xyz.x());
			loc_solver.normalizeAmpFuncsByPressureAtWall(dze);
			stab::dumpEigenFuncs("output/amp_funcs_direct.dat", loc_solver.getNNodes(), loc_solver.get_y_distrib(), dze);
			getchar();
		}

		t_Complex da_dw, sp_hw, sp_ha;
		// test group velo calcs
		da_dw = loc_solver.calcDaDwSpat(wchars);
		wxLogMessage(_T("da_dw_fd = (%lf, %lf)"), da_dw.real(), da_dw.imag());
		ofstr << xyz.x()<<"\t"<< da_dw.real() << "\t" << da_dw.imag() << "\t";

		loc_solver.calcGroupVelocity_ScalProd(wchars);

		t_Complex vga_inv = 1.0 / wchars.vga;
		wxLogMessage(_T("da_dw_mat = (%lf, %lf)"), vga_inv.real(), vga_inv.imag());
		loc_solver.calcScalarProd_H1_HW(dze, ksi, sp_ha, sp_hw);
		t_Complex da_dw_mat = sp_hw/sp_ha;
		wxLogMessage(_T("da_dw_mat = (%lf, %lf)"), da_dw_mat.real(), da_dw_mat.imag());
		ofstr << da_dw_mat.real() << "\t" << da_dw_mat.imag() << "\n";

		getchar();

		loc_solver.normalizeAmpFuncsByPressureAtWall(dze);

		// get deriv of amp_fun along x
		_calc_amp_fun_deriv_dx(i, loc_solver, fun_l, fun_r, dze_ddx);

		// restore context at this point after d_dx computations
		loc_solver.setContext(xyz);

		wchars = _line[i].wchars_loc;

		loc_solver.searchWave(wchars, stab::t_LSCond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED),
			stab::t_TaskTreat::SPAT);

		// calc <H1*dze_ddx, ksi>
		loc_solver.calcScalarProd_H1_HW(dze_ddx, ksi, v1_a, v1_w);

		// calc <H2*dze, ksi>
		v2_a = 0.0;// loc_solver.calcScalarProd_H2(dze, ksi);

		// calc <H1*dze, ksi>
		loc_solver.calcScalarProd_H1_HW(dze, ksi, v3_a, v3_w);

		wxLogMessage(_T("v1=(%lf, %lf)"), v1_a.real(), v1_a.imag());
		wxLogMessage(_T("v2=(%lf, %lf)"), v2_a.real(), v2_a.imag());
		wxLogMessage(_T("v3=(%lf, %lf)"), v3_a.real(), v3_a.imag());

		W = -(v1_a + v2_a) / v3_a;

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

