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
	std::vector<t_VecCmplx>& fun_m, std::vector<t_VecCmplx>& fun_r, std::vector<t_VecCmplx>& amp_funcs_deriv, t_QmData& QmData) {

	// TODO: correct expressions for first and last points
	// for now using second-order approx from second and last but one derivs

	if (_line.size() <= 1) wxLogMessage(_T("Error: can't calculate amp fun derivs with only one point in wpline!"));

	if (i == 0) i = 1;
	if (i == _line.size() - 1) i = _line.size() - 2;

	t_GeomPoint pp, p1, p2;

	pp = _line[i].mean_flow.get_xyz();
	p1 = _line[i - 1].mean_flow.get_xyz();
	p2 = _line[i + 1].mean_flow.get_xyz();

	double coef = _line[i].wchars_loc.scales().Dels/_rFldMF.get_mf_params().L_ref;
	//t_StabScales scales = loc
	const double dx = (p2 - p1).norm()/coef;
	const double c = 1.0 / dx;

	// set context to get thickness to be used as thick for +- calcs
	wxLogMessage(_T("PP:x=%lf, y=%lf,z=%lf"), pp.x(), pp.y(), pp.z());
	loc_solver.setContext(pp);

	double y_thick = loc_solver.getThickMF();

	mf::t_ProfDataCfg prof_cfg;

	prof_cfg.ThickFixed = y_thick;

	// get p wall to non dim amp funcs in left and right points
	get_amp_funcs(i, loc_solver, fun_m, &prof_cfg);

	int nnodes_stab = fun_m.size();
	t_Complex pwall_mid = fun_m[nnodes_stab - 1][3];
	loc_solver.normalizeAmpFuncsByFixedVal(fun_m, pwall_mid);

	loc_solver.calcQmAmpFun(fun_m, QmData.qm_m);

	// find index where qm=max
	int j_qm_max = 0;
	QmData.qmabs_max = 0.0;
	for (int j = 0; j < QmData.qm_m.size(); j++) {
		if (smat::norm(QmData.qm_m[j]) > QmData.qmabs_max) {
			QmData.qmabs_max = smat::norm(QmData.qm_m[j]);
			j_qm_max = j;
		}
	}

	// debug
	wxLogMessage(_T("J where Qm is max:%d"), j_qm_max);

	// debug
	std::vector<double> yy_l(nnodes_stab);
	std::vector<double> yy_r(nnodes_stab);

	// get amp funcs with fixed thickness of prof stab
	get_amp_funcs(i - 1, loc_solver, fun_l, &prof_cfg);
	yy_l = loc_solver.get_y_distrib();
	loc_solver.normalizeAmpFuncsByFixedVal(fun_l, pwall_mid);
	loc_solver.calcQmAmpFun(fun_l, QmData.qm_l);


	get_amp_funcs(i + 1, loc_solver, fun_r, &prof_cfg);
	yy_r = loc_solver.get_y_distrib();
	loc_solver.normalizeAmpFuncsByFixedVal(fun_r, pwall_mid);
	loc_solver.calcQmAmpFun(fun_r, QmData.qm_r);

	// debug
	if (i == 158) {
		std::ofstream ofstr("output/qm_ddx.dat");
		for (int j = 0; j < nnodes_stab; j++) {
			double dqmabs_dx = c*(smat::norm(QmData.qm_r[j]) - smat::norm(QmData.qm_l[j]));
			ofstr << yy_l[j] << "\t" << dqmabs_dx << "\n";
		}
	}

	// calculate qm deriv at point where qm=max 
	{

		QmData.dqmabs_dx = c*(smat::norm(QmData.qm_r[j_qm_max]) - smat::norm(QmData.qm_l[j_qm_max]));

	}

	const int stab_matrix_dim = 8;

	for (int j = 0; j < nnodes_stab; j++) {
		//if (yy_l[j] != yy_r[j]) { wxLogMessage(_T("Error:y_l=%lf, y_r=%lf"), yy_l[j], yy_r[j]); getchar(); }
		for (int i = 0; i < stab_matrix_dim; i++)
			amp_funcs_deriv[j][i] = c*(fun_r[j][i] - fun_l[j][i]);
	}

	static bool do_write_ddx = true;

	if (pp.x() > 0.44 && do_write_ddx) {
		do_write_ddx = false;

		dumpEigenFuncs("output/dze_l.dat", nnodes_stab, yy_l, fun_l);
		dumpEigenFuncs("output/dze_r.dat", nnodes_stab, yy_r, fun_r);
		dumpEigenFuncs("output/ddze_dx.dat", nnodes_stab, yy_l, amp_funcs_deriv);

		// dump qm_mid distribution (from wall to outer rec)
		std::ofstream ofstr("output/dze_qm.dat");
		for (int j = 0; j < QmData.qm_m.size(); j++) 
			ofstr << yy_l[QmData.qm_m.size() - 1 - j] << "\t" << smat::norm(QmData.qm_m[QmData.qm_m.size() - 1 - j]) << "\n";
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
	std::vector<t_VecCmplx> fun_m(loc_solver.getNNodes(), t_VecCmplx(8));
	std::vector<t_VecCmplx> fun_r(loc_solver.getNNodes(), t_VecCmplx(8));
	std::vector<t_VecCmplx> dze_ddx(loc_solver.getNNodes(), t_VecCmplx(8));
	// direct amplitude func
	std::vector<t_VecCmplx> dze(loc_solver.getNNodes(), t_VecCmplx(8));
	// conjugate amplitude func 
	std::vector<t_VecCmplx> ksi(loc_solver.getNNodes(), t_VecCmplx(8));

	t_QmData qm_data;
	qm_data.resize(loc_solver.getNNodes());

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
//		if (i == 157) {
//			wxLogMessage(_T("Current point x=%lf"), xyz.x());
//			loc_solver.normalizeAmpFuncsByPressureAtWall(dze);
//			stab::dumpEigenFuncs("output/amp_funcs_direct.dat", loc_solver.getNNodes(), loc_solver.get_y_distrib(), dze);
//			getchar();
//		}

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

		//getchar();

		loc_solver.normalizeAmpFuncsByPressureAtWall(dze);

		// get deriv of amp_fun along x
		_calc_amp_fun_deriv_dx(i, loc_solver, fun_l, fun_m, fun_r, dze_ddx, qm_data);

		// restore context at this point after d_dx computations
		loc_solver.setContext(xyz);

		wchars = _line[i].wchars_loc;

		loc_solver.searchWave(wchars, stab::t_LSCond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED),
			stab::t_TaskTreat::SPAT);

		// calc <H1*dze_ddx, ksi>
		loc_solver.calcScalarProd_H1_HW(dze_ddx, ksi, v1_a, v1_w);

		//v1_a = 0.0;

		// calc <H2*dze, ksi>
		v2_a = loc_solver.calcScalarProd_H2(dze, ksi);

		//v2_a = 0.0;

		// calc <H1*dze, ksi>
		loc_solver.calcScalarProd_H1_HW(dze, ksi, v3_a, v3_w);

		wxLogMessage(_T("v1=(%lf, %lf)"), v1_a.real(), v1_a.imag());
		wxLogMessage(_T("v2=(%lf, %lf)"), v2_a.real(), v2_a.imag());
		wxLogMessage(_T("v3=(%lf, %lf)"), v3_a.real(), v3_a.imag());

		W = -(v1_a + v2_a) / v3_a;

		wxLogMessage(_T("Uniform non par addition W=%lf"), W.real());

		double Da_Qm = 1.0 / qm_data.qmabs_max*qm_data.dqmabs_dx;

		wxLogMessage(_T("Da_Qm=%lf"), Da_Qm);

		// dump qm at station x=0.926 <=> i=158
		if (i == 158) {

			const std::vector<double>& y_vec = loc_solver.get_y_distrib();
			// debug, dump qm
			std::ofstream ofstr("output/qm.dat");
			for (int j = 0; j < qm_data.qm_m.size(); j++) {

				ofstr << y_vec[j] << "\t" << smat::norm(qm_data.qm_m[j]) << "\n";
			}

		}

		if (_params.CalcNonParEffectsAtQmax) {
			_line[i].da_nonpar = -1.0*E*(W + Da_Qm);
		}
		else {
			_line[i].da_nonpar = -1.0*E*W;
		}

		// debug 
		wxLogMessage(_T("da_nonpar = (%lf, %lf)"), 
			_line[i].da_nonpar.real(), _line[i].da_nonpar.imag());

		da_ratio = _line[i].da_nonpar / wchars.a;

		wxLogMessage(_T("Sigma addition (nondim local):%lf"), -1.0*_line[i].da_nonpar.imag());
		wxLogMessage(_T("Debug : current x=%lf"), xyz.x());
		getchar();

		// modify wave chars glob and DO NOT modify wchars loc!
		t_WCharsLoc wchars_new;
		wchars_new = _line[i].wchars_loc;
		wchars_new.a = _line[i].wchars_loc.a + _line[i].da_nonpar;

		t_WCharsGlob wchars_glob(wchars_new, _rFldMF.calc_jac_to_loc_rf(xyz),
			loc_solver.get_stab_scales());

		_line[i].wave_chars = wchars_glob;	

	}

}

