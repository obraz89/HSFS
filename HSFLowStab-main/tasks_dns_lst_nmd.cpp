#include "stdafx.h"

#include "tasks.h"

#include "ProfileStab.h"

#include "solvers_glob.h"

#include "fstream"

using namespace hsstab;

void task::calc_scal_prod_self() {

	wxLogMessage(_T("Scal prod LST test: both modes from LST"));

	const mf::t_GeomPoint& xyz = g_pStabDB->get_pave_pt(0).xyz;

	g_pStabSolver->setContext(xyz);

	g_pGSSolverSpat->setContext(xyz);

	g_pGSSolverTime->setContext(xyz);

	std::vector<t_WCharsLoc> init_waves_raw;
	std::vector<t_WCharsLoc> init_waves_filtered;

	t_WCharsLoc init_wave, ret_wave;

	//init_wave.b = g_taskParams.b_ndim_min;
	//init_wave.w = g_taskParams.w_ndim_min;

	//init_wave.w = t_Complex(0.562, 0.0);
	//init_wave.b = 0.91;

	//init_waves_raw = g_pGSSolverSpat->getInstabModes(init_wave);

	//init_waves_filtered = g_pStabSolver->filter_gs_waves_spat(init_waves_raw, 
	//	stab::t_LSCond(stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED));

	//ret_wave = t_WCharsLoc::find_max_instab_spat(init_waves_filtered); 

	//wxLogMessage(_T("wave chars for direct problem:\n%s"), &(ret_wave.to_wstr()[0]) );
	//return;

	t_WCharsLoc wchars_A, wchars_B;

	// mode "A"

	//wchars_A.a = t_Complex(0.898, -0.0121);
	//wchars_A.b = 1.22;
	//wchars_A.w = t_Complex(0.762, 0.0);

	wchars_A.a = t_Complex(0.682, -0.0105);;
	wchars_A.b = 0.91;
	wchars_A.w = t_Complex(0.562, 0.0);


	wchars_A.set_treat(stab::t_TaskTreat::SPAT);

	wchars_A.resid = 1.0E+06;

	// 4 debugging
	//g_pStabSolver->setLSMode(stab::t_LSMode(stab::t_LSMode::DIRECT|stab::t_LSMode::ASYM_HOMOGEN));

	g_pStabSolver->searchWave(wchars_A, stab::t_LSCond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED),
		stab::t_TaskTreat::SPAT);

	g_pStabSolver->dumpEigenFuctions("output/mode_A.dat");

	//wxLogMessage(_T("Wave chars A (direct):%s"), &(wchars_A.to_wstr()[0]));


	// mode "B"

	//	wchars_B.a = t_Complex(0.5157, 0.000645);
	//	wchars_B.b = 1.22;
	//	wchars_B.w = t_Complex(0.762, 0.0);

	wchars_B.a = t_Complex(0.378, 0.00025);
	wchars_B.b = 0.91;
	wchars_B.w = t_Complex(0.562, 0.0);

	wchars_B.set_treat(stab::t_TaskTreat::SPAT);

	wchars_B.resid = 1.0E+06;

	/*
	double min_resi = 1.0e+06;
	double a_re, a_im;

	for (int i=0; i<100; i++){
	for (int j=0; j<100; j++)
	{
	wchars_B.a = t_Complex(0.4+i*0.002, -0.1+j*0.002);

	wchars_B.resid = 1.0E+06;

	t_Complex resi = g_pStabSolver->solve(wchars_B);

	double rr = sqrt(resi.real()*resi.real()+resi.imag()*resi.imag());

	wxLogMessage(_T("a=(%lf,%lf):res = %lf"),wchars_B.a.real(), wchars_B.a.imag(),
	rr);

	if (rr<min_resi){
	min_resi = rr;
	a_re = wchars_B.a.real();
	a_im = wchars_B.a.imag();
	}
	}

	}

	wxLogMessage(_T("Minimum: a=(%lf,%lf):res = %lf"),a_re, a_im, min_resi);
	return;*/

	g_pStabSolver->searchWave(wchars_B, stab::t_LSCond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED),
		stab::t_TaskTreat::SPAT);

	g_pStabSolver->dumpEigenFuctions("output/mode_B.dat");

	t_Complex val = g_pStabSolver->calcScalarProd(wchars_A, wchars_B);

	wxLogMessage(_("Scalar product result:(%f,%f)"), val.real(), val.imag());

	int vvvv = 1;

}

// read bl data required to initialize DNS disturbance profile

void _get_bl_info_from_file(std::string fname_bl_data, t_ProfStabCfgDNSDisturb& pstab_cfg) {

	std::ifstream ifstr(&fname_bl_data[0]);
	std::stringstream istr;

	const int max_lsize = 256;
	char line[max_lsize];
	char ch;

	// first line is header
	ifstr.get(line, max_lsize, '\n');
	ifstr.get(ch);

	ifstr.get(line, max_lsize, '\n');

	istr << line;

	istr >> pstab_cfg.Ymax;

	istr >> pstab_cfg.NNodes;

	istr >> pstab_cfg.Dels_gndim;

	istr >> pstab_cfg.Ue;

	istr >> pstab_cfg.Rhoe;

	istr >> pstab_cfg.Te;

	ifstr.close();

}

t_Complex _search_max_pressure_disturb(const std::vector<t_VecCmplx>& amp_funcs){

	double max_abs = 0.0;

	t_Complex max_val;

	double cur_abs;

	for (int i = 0; i < amp_funcs.size(); i++) {

		cur_abs = smat::norm(amp_funcs[i][3]);

		if (cur_abs > max_abs) {
			max_abs = cur_abs;
			max_val = amp_funcs[i][3];
		}

	}

	return max_val;

}

// helper function to read DNS amplitude vectors from file

void _get_DNS_amp_vecs_from_file(std::string fname_re, std::string fname_im, std::string fname_bl_data,
	std::vector<t_VecCmplx>& amp_vec_out, std::vector<double>& yy_dns) {

	t_ProfStabCfgDNSDisturb pstab_cfg;

	_get_bl_info_from_file(fname_bl_data, pstab_cfg);

	mf::t_ProfDataCfg data_cfg;

	data_cfg.ThickCoef = pstab_cfg.Ymax / pstab_cfg.Dels_gndim;

	t_ProfileNS prof_ns_re(*g_pMFDomain);
	prof_ns_re.initialize_extract(fname_re, data_cfg);

	t_ProfileNS prof_ns_im(*g_pMFDomain);
	prof_ns_im.initialize_extract(fname_im, data_cfg);

	t_ProfileStab ps_re, ps_im;

	ps_re.initialize_dist_DNS(prof_ns_re, pstab_cfg);
	ps_im.initialize_dist_DNS(prof_ns_im, pstab_cfg);

	int nnodes_stab = pstab_cfg.NNodes;

	if (nnodes_stab != amp_vec_out.size())
		wxLogError(_T("Error while processing DNS amp funcs from file: size mismatch"));

	// amp vecs in stability code are from outer to wall !!!
	for (int i = 0; i<nnodes_stab; i++)
	{

		int i_stab = nnodes_stab - 1 - i;

		// y
		yy_dns[i] = ps_re.get_rec(i).y;

		// (u_re, u_im)
		amp_vec_out[i_stab][0] = t_Complex(ps_re.get_rec(i).u, ps_im.get_rec(i).u);

		// (u'_re, u'_im)
		amp_vec_out[i_stab][1] = t_Complex(ps_re.get_rec(i).u1, ps_im.get_rec(i).u1);

		// (v_re, v_im)
		amp_vec_out[i_stab][2] = t_Complex(ps_re.get_rec(i).v, ps_im.get_rec(i).v);

		// (p_re, p_im)
		amp_vec_out[i_stab][3] = t_Complex(ps_re.get_rec(i).p, ps_im.get_rec(i).p);

		// (t_re, t_im)
		amp_vec_out[i_stab][4] = t_Complex(ps_re.get_rec(i).t, ps_im.get_rec(i).t);

		// (t'_re, t'_im)
		amp_vec_out[i_stab][5] = t_Complex(ps_re.get_rec(i).t1, ps_im.get_rec(i).t1);

		// (w_re, w_im)
		amp_vec_out[i_stab][6] = t_Complex(ps_re.get_rec(i).w, ps_im.get_rec(i).w);

		// (w'_re, w'_im)
		amp_vec_out[i_stab][7] = t_Complex(ps_re.get_rec(i).w1, ps_im.get_rec(i).w1);

	}



	// debug, check amp_vec_out
	
	std::ofstream ofstr("output/amp_vec_DNS.dat");

	for (int i = 0; i<nnodes_stab; i++) {
		ofstr << ps_re.get_rec(i).y << "\t";
		for (int j = 0; j<8; j++) {
			ofstr << amp_vec_out[nnodes_stab - 1 - i][j].real() << "\t";
			ofstr << amp_vec_out[nnodes_stab - 1 - i][j].imag() << "\t";
		}
		ofstr << "\n";
	}
	ofstr.close();
	
}

void task::calc_scal_prod_particle_test() {

	std::string fname_bl_data("spectrum_input/amp_funcs_bl_data.dat");

	std::string fname_re("spectrum_input/amp_funcs_dnsx0.2_om-24_b80_re.dat");
	std::string fname_im("spectrum_input/amp_funcs_dnsx0.2_om-24_b80_im.dat");

	const int STAB_MATRIX_DIM = 8;

	const int nnodes_stab = 501;

	std::vector<double> yy_DNS(nnodes_stab);
	std::vector<t_VecCmplx> amp_fun_DNS(nnodes_stab, STAB_MATRIX_DIM);

	// x~0.2
	// DNS wave packet initiated by particle near x~0.13:
	// packet hump at Om_glob_ndim = 150, B_glob_ndim = 500
	// get lst amp funcs for specified Om and Beta,
	// normalize them
	// and calculate scal prod coefficient
	// A - direct task, B - conjugate, Q - DNS (direct)
	// <Q_dns, B>/<A, B>
	// Om_glob = 150 <=> w = 0.181
	// B_glob = 500 <=> b = 0.55

	_get_DNS_amp_vecs_from_file(fname_re, fname_im, fname_bl_data, amp_fun_DNS, yy_DNS);

	t_Complex pmax_DNS = _search_max_pressure_disturb(amp_fun_DNS);

	const mf::t_GeomPoint& xyz = g_pStabDB->get_pave_pt(0).xyz;

	g_pStabSolver->setContext(xyz);

	t_WCharsLoc wchars_A;

	wchars_A.a = t_Complex(0.3029, -0.01141);
	wchars_A.b = 0.55;
	wchars_A.w = t_Complex(0.181, 0.0);

	wchars_A.set_treat(stab::t_TaskTreat::SPAT);

	// direct problem

	g_pStabSolver->setLSMode(stab::t_LSMode(stab::t_LSMode::DIRECT | stab::t_LSMode::ASYM_HOMOGEN));

	g_pStabSolver->searchWave(wchars_A, stab::t_LSCond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED),
		stab::t_TaskTreat::SPAT);
	
	wxLogMessage(_T("Wave chars A (direct):%s"), &(wchars_A.to_wstr()[0]));

	//g_pStabSolver->dumpEigenFuctions("output/amp_funcs_lst");

	std::vector<t_VecCmplx> amp_fun_LST(nnodes_stab, STAB_MATRIX_DIM);

	g_pStabSolver->getAmpFuncs(amp_fun_LST);

	t_Complex pmax_LST = _search_max_pressure_disturb(amp_fun_LST);

	t_Complex coef = pmax_DNS / pmax_LST;

	wxLogMessage(_T("pmax_DNS/pmax_LST=(%lf,%lf)"), coef.real(), coef.imag());
	wxLogMessage(_T("abs(pmax_DNS/pmax_LST)=%lf"), smat::norm(coef));

	// scale lst amp_funcs to dns amp funcs

	std::ofstream ofstr("output/amp_vec_lst_scaled2dns.dat");

	for (int i = 0; i<nnodes_stab; i++) {

		ofstr << yy_DNS[i] << "\t";
		for (int j = 0; j<8; j++) {
			t_Complex val = amp_fun_LST[nnodes_stab - 1 - i][j] * coef;
			ofstr << val.real() << "\t";
			ofstr << val.imag() << "\t";
		}
		ofstr << "\n";
	}
	ofstr.close();

	t_Complex val_lst = coef*g_pStabSolver->calcScalarProd(wchars_A, wchars_A, NULL);

	wxLogMessage(_T("LST <A_lst_scaled2dns, A_lst_conj>=%lf"), smat::norm(val_lst));

	t_Complex val_dns = g_pStabSolver->calcScalarProd(wchars_A, wchars_A, &amp_fun_DNS);

	wxLogMessage(_T("DNS <Q_dns, A_lst_conj>=%lf"), smat::norm(val_dns));

}


void task::calc_scal_prod_dns_test() {

	const mf::t_GeomPoint& xyz = g_pStabDB->get_pave_pt(0).xyz;

	g_pStabSolver->setContext(xyz);

	g_pGSSolverSpat->setContext(xyz);

	g_pGSSolverTime->setContext(xyz);

	// "A" - direct problem, taking vector from DNS
	// "B" - solving conjugate problem, taking mode from LST solver
	// scal prod matrix is set with eigenvalues from "B"

	t_WCharsLoc wchars_A, wchars_B;

	int nnodes_stab = g_pStabSolver->getNNodes();


	wchars_A.a = t_Complex(0.682, -0.0105);;
	wchars_A.b = 0.91;
	wchars_A.w = t_Complex(0.562, 0.0);

	wchars_A.set_treat(stab::t_TaskTreat::SPAT);


	wchars_B.a = t_Complex(0.378, 0.00025);
	wchars_B.b = 0.91;
	wchars_B.w = t_Complex(0.562, 0.0);

	//wchars_B = wchars_A;

	wchars_B.set_treat(stab::t_TaskTreat::SPAT);

	// get DNS amplitude funcs

	const int STAB_MATRIX_DIM = 8;
	std::vector<t_VecCmplx> sol_dir_A(nnodes_stab, STAB_MATRIX_DIM);

	//_get_DNS_amp_vecs_from_file(sol_dir_A, nnodes_stab);

	// setting scal prod matrix 

	t_Complex val;

	val = g_pStabSolver->calcScalarProd(wchars_A, wchars_B, &sol_dir_A);

	wxLogMessage(_("Scalar product result:(%f,%f)"), val.real(), val.imag());

	int vvvv = 1;

}