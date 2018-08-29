#include "stdafx.h"

#include "tasks.h"

#include "ProfileStab.h"

#include "solvers_glob.h"

#include "fstream"

using namespace hsstab;

struct t_SpctrInfo {

	int Nw, Nb, Ny;

	double Om_fft_min, Om_fft_max;

	double B_fft_min, B_fft_max;

};

// iinformation required from Mean-Flow profiles
// to correctly initialize disturbance profiles
// Dels gndim is Dels/L_ref from StabProfile
// Ue, Ve, We - direction of velocity at bl edge
// Eta_max is profile thickness in gndim
//
struct t_MFProfInfo {

	double Dels_gndim;

	int NNodesStab;

	double Ue_dir, Ve_dir, We_dir;

	double Eta_max;

	double Ue_abs, Rho_e, T_e;

};

void _get_spctr_info(t_SpctrInfo& spd, std::string fname_spctr_info) {

	std::ifstream ifstr(fname_spctr_info);
	std::stringstream istr;

	const int max_line_size = 256;
	char line[max_line_size];

	// header
	ifstr.getline(line, max_line_size, '\n');
	//data
	ifstr.getline(line, max_line_size, '\n');

	istr << line;

	istr >> spd.Ny >> spd.Nw >> spd.Nb;

	istr >> spd.Om_fft_min >> spd.Om_fft_max;

	istr >> spd.B_fft_min >> spd.B_fft_max;

	ifstr.close();

}

// TODO: all info from 1 file
void _get_mfprof_info(t_MFProfInfo& mfd, std::string fname_mf_info) {

	std::ifstream ifstr(fname_mf_info);
	std::stringstream istr;

	const int max_line_size = 256;
	char line[max_line_size];

	// header
	ifstr.getline(line, max_line_size, '\n');
	//data
	ifstr.getline(line, max_line_size, '\n');

	istr << line;

	istr >> mfd.Eta_max >> mfd.NNodesStab >> mfd.Dels_gndim;

	istr >> mfd.Ue_abs >> mfd.Rho_e >> mfd.T_e;

	ifstr.getline(line, max_line_size, '\n');

	istr.clear();
	istr << line;

	istr >> mfd.Ue_dir >> mfd.Ve_dir >> mfd.We_dir;

	ifstr.close();

}

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
	
}

void _get_lst_wpacket_from_dns();

// initial testing for a point near spectral maximum for 
// a wwave packet near x=0.2
void _do_calc_scal_prod_spctr_max_point() {

	std::string fname_bl_data("spectrum_input/_amp_funcs_bl_data.dat");

	std::string fname_re("spectrum_input/_amp_funcs_tmp_re.dat");
	std::string fname_im("spectrum_input/_amp_funcs_tmp_im.dat");

	const int STAB_MATRIX_DIM = 8;

	const int nnodes_stab = 501;

	std::vector<double> yy_DNS(nnodes_stab);
	std::vector<t_VecCmplx> amp_fun_DNS(nnodes_stab, STAB_MATRIX_DIM);

	// calculate scal prod coefficient
	// A - direct task, B - conjugate, Q - DNS (direct)
	// <Q_dns, B>/<A, B>

	_get_DNS_amp_vecs_from_file(fname_re, fname_im, fname_bl_data, amp_fun_DNS, yy_DNS);

	// debug, check amp_vec_out

	std::ofstream ofstr_dns("output/amp_vec_DNS.dat");

	for (int i = 0; i<nnodes_stab; i++) {
		ofstr_dns << yy_DNS[i] << "\t";
		for (int j = 0; j<8; j++) {
			ofstr_dns << amp_fun_DNS[nnodes_stab - 1 - i][j].real() << "\t";
			ofstr_dns << amp_fun_DNS[nnodes_stab - 1 - i][j].imag() << "\t";
		}
		ofstr_dns << "\n";
	}
	ofstr_dns.close();

	t_Complex pmax_DNS = _search_max_pressure_disturb(amp_fun_DNS);

	const mf::t_GeomPoint& xyz = g_pStabDB->get_pave_pt(0).xyz;

	g_pStabSolver->setContext(xyz);

	t_WCharsLoc wchars_A;

	bool read_ok = read_max_wave_pid(0, _T("wchars_max_loc.dat"), wchars_A);

	if (!read_ok) wxLogMessage(_T("cannot read wchars from wchars_max_loc.dat"));

	wchars_A.set_treat(stab::t_TaskTreat::SPAT);

	// direct problem

	g_pStabSolver->setLSMode(stab::t_LSMode(stab::t_LSMode::DIRECT | stab::t_LSMode::ASYM_HOMOGEN));

	g_pStabSolver->searchWave(wchars_A, stab::t_LSCond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED),
		stab::t_TaskTreat::SPAT);

	wxLogMessage(_T("Wave chars A (direct):%s"), &(wchars_A.to_wstr()[0]));

	//g_pStabSolver->dumpEigenFuctions("output/amp_funcs_lst");

	std::vector<t_VecCmplx> amp_fun_LST(nnodes_stab, STAB_MATRIX_DIM);

	g_pStabSolver->getAmpFuncs(amp_fun_LST);

	const std::vector<double>& y_lst = g_pStabSolver->get_y_distrib();

	t_Complex pmax_LST = _search_max_pressure_disturb(amp_fun_LST);

	t_Complex coef_p = pmax_DNS / pmax_LST;

	wxLogMessage(_T("pmax_DNS/pmax_LST=(%lf,%lf)"), coef_p.real(), coef_p.imag());
	wxLogMessage(_T("abs(pmax_DNS/pmax_LST)=%lf"), smat::norm(coef_p));

	// scale lst amp_funcs to dns amp funcs by pressure coef

	std::ofstream ofstr("output/amp_vec_lst_scaled2dns_by_pressure.dat");

	for (int i = 0; i<nnodes_stab; i++) {

		ofstr << y_lst[nnodes_stab - 1 - i] << "\t";
		for (int j = 0; j<8; j++) {
			t_Complex val = amp_fun_LST[nnodes_stab - 1 - i][j] * coef_p;
			ofstr << val.real() << "\t";
			ofstr << val.imag() << "\t";
		}
		ofstr << "\n";
	}
	ofstr.close();

	t_Complex val_lst = g_pStabSolver->calcScalarProd(wchars_A, wchars_A, NULL);

	wxLogMessage(_T("LST <A_lst_scaled2dns_p, A_lst_conj>=%lf"), smat::norm(coef_p*val_lst));

	t_Complex val_dns = g_pStabSolver->calcScalarProd(wchars_A, wchars_A, &amp_fun_DNS);

	wxLogMessage(_T("DNS <Q_dns, A_lst_conj>=%lf"), smat::norm(val_dns));

	// calculate (<Q_DNS, A_conj>/<A, A_conj>)*A

	t_Complex coef_scal = val_dns / val_lst;

	std::ofstream ofstr1("output/amp_vec_lst_scaled2dns_by_sprod.dat");

	for (int i = 0; i<nnodes_stab; i++) {

		ofstr1 << y_lst[nnodes_stab - 1 - i] << "\t";
		for (int j = 0; j<8; j++) {
			t_Complex val = amp_fun_LST[nnodes_stab - 1 - i][j] * coef_scal;
			ofstr1 << val.real() << "\t";
			ofstr1 << val.imag() << "\t";
		}
		ofstr1 << "\n";
	}
	ofstr1.close();

}

// get amp funcs from bunch of spectrum files
// output is a pair of intermediate files _amp_funcs_tmp_*.dat *={re,im}
// which are then used to initialize amp functions in LST solver basis
void _get_amp_funcs_from_full_spectrum(int i_w, int j_b, const t_SpctrInfo spd, const t_MFProfInfo mfd) {

	char fname_sp_j_fix[128];

	const int max_lsize = 1024;
	char line[max_lsize];
	char ch;

	char fname_sp_base[] = "spectrum_input/amp_funcs_dns_ombeta";

	std::ifstream ifstr;

	std::stringstream istr;

	const int glob_ind = i_w*spd.Nb + j_b;

	int iw_read, jb_read;

	double om_fft_read, b_fft_read;

	double x, y, z, ur, ui, vr, vi, wr, wi, pr, pi, tr, ti;

	std::ifstream ifstr_xyz("spectrum_input/_xyz_data.dat");

	std::ofstream ofstr_re("spectrum_input/_amp_funcs_tmp_re.dat");
	std::ofstream ofstr_im("spectrum_input/_amp_funcs_tmp_im.dat");

	ofstr_re << spd.Ny << "\t" << mfd.Dels_gndim << "\n";
	ofstr_re << mfd.Ue_dir << "\t" << mfd.Ve_dir << "\t" << mfd.We_dir << "\n";

	ofstr_im << spd.Ny << "\t" << mfd.Dels_gndim << "\n";
	ofstr_im << mfd.Ue_dir << "\t" << mfd.Ve_dir << "\t" << mfd.We_dir << "\n";

	for (int j = 0; j < spd.Ny; j++) {

		ifstr_xyz.getline(line, max_lsize, '\n');
		istr.clear();
		istr << line;

		istr >> x >> y >> z;

		sprintf(fname_sp_j_fix, "%s_j_%d.dat", fname_sp_base, j);

		//wxLogMessage(_T("reading spectrum %s"), wxString::FromAscii(fname_sp_j_fix));

		ifstr.open(fname_sp_j_fix);

		// read Nw, Nb line - just skipping
		ifstr.getline(line, max_lsize, '\n');

		// skip all data up to a required line
		for (int ln = 0; ln < glob_ind; ln++)
			ifstr.getline(line, max_lsize, '\n');

		ifstr.getline(line, max_lsize, '\n');

		istr.clear();
		istr << line;

		istr >> iw_read >> jb_read;
		istr >> om_fft_read >> b_fft_read;

		if (i_w != iw_read || j_b != jb_read) {
			wxLogMessage(_T("Error while parsing spectrum data file: iw=%d, iw_read=%d, jb=%d, jb_read=%d"), i_w, iw_read, j_b, jb_read);
		}

		istr >> ur >> ui;
		istr >> vr >> vi;
		istr >> wr >> wi;
		istr >> pr >> pi;
		istr >> tr >> ti;

		ofstr_re << x << "\t" << y << "\t" << z << "\t";
		ofstr_re << ur << "\t" << vr << "\t" << wr << "\t" << pr << "\t" << tr << "\n";

		ofstr_im << x << "\t" << y << "\t" << z << "\t";
		ofstr_im << ui << "\t" << vi << "\t" << wi << "\t" << pi << "\t" << ti << "\n";

		ifstr.close();
	}

	ifstr_xyz.close();
}

void _do_calc_dns2lst_wall_disturb(t_WCharsLoc wchars_A, std::vector<t_Complex>& a_wall_data) {

	std::string fname_bl_data("spectrum_input/_amp_funcs_bl_data.dat");

	std::string fname_re("spectrum_input/_amp_funcs_tmp_re.dat");
	std::string fname_im("spectrum_input/_amp_funcs_tmp_im.dat");

	const int STAB_MATRIX_DIM = 8;

	const int nnodes_stab = 501;

	std::vector<double> yy_DNS(nnodes_stab);
	std::vector<t_VecCmplx> amp_fun_DNS(nnodes_stab, STAB_MATRIX_DIM);

	_get_DNS_amp_vecs_from_file(fname_re, fname_im, fname_bl_data, amp_fun_DNS, yy_DNS);

	// direct problem

	g_pStabSolver->setLSMode(stab::t_LSMode(stab::t_LSMode::DIRECT | stab::t_LSMode::ASYM_HOMOGEN));

	g_pStabSolver->searchWave(wchars_A, stab::t_LSCond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED),
		stab::t_TaskTreat::SPAT);

	wxLogMessage(_T("Wave chars A (direct):%s"), &(wchars_A.to_wstr()[0]));

	std::vector<t_VecCmplx> amp_fun_LST(nnodes_stab, STAB_MATRIX_DIM);

	g_pStabSolver->getAmpFuncs(amp_fun_LST);

	if (true) {

		// debug - write amplitude functions

		t_Complex pmax_DNS = _search_max_pressure_disturb(amp_fun_DNS);

		t_Complex pmax_LST = _search_max_pressure_disturb(amp_fun_LST);

		t_Complex coef_p = pmax_DNS / pmax_LST;

		wxLogMessage(_T("pmax_DNS/pmax_LST=(%lf,%lf)"), coef_p.real(), coef_p.imag());
		wxLogMessage(_T("abs(pmax_DNS/pmax_LST)=%lf"), smat::norm(coef_p));

		std::ofstream ofstr_dns("output/amp_vec_DNS.dat");

		for (int i = 0; i<nnodes_stab; i++) {
			ofstr_dns << yy_DNS[i] << "\t";
			for (int j = 0; j<8; j++) {
				ofstr_dns << amp_fun_DNS[nnodes_stab - 1 - i][j].real() << "\t";
				ofstr_dns << amp_fun_DNS[nnodes_stab - 1 - i][j].imag() << "\t";
			}
			ofstr_dns << "\n";
		}
		ofstr_dns.close();

		std::ofstream ofstr("output/amp_vec_lst_scaled2dns_by_pressure.dat");

		const std::vector<double>& y_lst = g_pStabSolver->get_y_distrib();

		for (int i = 0; i<nnodes_stab; i++) {

			ofstr << y_lst[nnodes_stab - 1 - i] << "\t";
			for (int j = 0; j<8; j++) {
				t_Complex val = amp_fun_LST[nnodes_stab - 1 - i][j] * coef_p;
				ofstr << val.real() << "\t";
				ofstr << val.imag() << "\t";
			}
			ofstr << "\n";
		}
		ofstr.close();	

	}

	if(true){

		// version 2 - calculate (<Q_DNS, A_conj>/<A, A_conj>)*A

		t_Complex val_lst = g_pStabSolver->calcScalarProd(wchars_A, wchars_A, NULL);

		t_Complex val_dns = g_pStabSolver->calcScalarProd(wchars_A, wchars_A, &amp_fun_DNS);

		//wxLogMessage(_T("LST <A_lst_scaled2dns_p, A_lst_conj>=%lf"), smat::norm(coef_p*val_lst));

		t_Complex coef_scal = val_dns / val_lst;

		wxLogMessage(_T("<Q_dns, A_lst_conj>/<A_lst, A_lst_conj>=(%lf, %lf)"), coef_scal.real(), coef_scal.imag());
		wxLogMessage(_T("abs(coef_scal)=%lf"), smat::norm(coef_scal));

		// du
		a_wall_data[0] = amp_fun_LST[nnodes_stab - 1][0] * coef_scal;
		// dv
		a_wall_data[1] = amp_fun_LST[nnodes_stab - 1][2] * coef_scal;
		// dw
		a_wall_data[2] = amp_fun_LST[nnodes_stab - 1][6] * coef_scal;
		// dp
		a_wall_data[3] = amp_fun_LST[nnodes_stab - 1][3] * coef_scal;
		// dt
		a_wall_data[4] = amp_fun_LST[nnodes_stab - 1][4] * coef_scal;
	
	}

}

t_WCharsLoc _get_wc_s_default(const t_WCharsLoc& wc_dest) {

	t_WCharsLoc wchars_A;

	bool read_ok = read_max_wave_pid(0, _T("wchars_max_loc.dat"), wchars_A);

	if (!read_ok) wxLogMessage(_T("cannot read wchars from wchars_max_loc.dat"));

	wchars_A.set_treat(stab::t_TaskTreat::SPAT);

	wchars_A.set_scales(g_pStabSolver->get_stab_scales());

	if (wc_dest.b.real() < 0.0) wchars_A.b = -1.0*wchars_A.b;

	if (wc_dest.w.real() < 0.0) {

		wchars_A.w = -1.0*wchars_A.w;
		wchars_A.a.real(-1.0*wchars_A.a.real());

	}

	return wchars_A;

}

// when a maximum of wave packet spectrum was found
// in form ws, bs, ar_s, ai_s
// get (ar, ai) eigen for a given values of w and b
// 
void _get_wchars(t_WCharsLoc& wc_dest, const t_WCharsLoc& wc_s) {

	// direct problem

	wxLogMessage(_T("_get_wchars() in process, this may take a while..."));

	wxLogMessage(_T("_get_wave:wchars input:%s"), &(wc_dest.to_wstr()[0]));
	wxLogMessage(_T("_get_wave:wchars starting point:%s"), &(wc_s.to_wstr()[0]));

	g_pStabSolver->setLSMode(stab::t_LSMode(stab::t_LSMode::DIRECT | stab::t_LSMode::ASYM_HOMOGEN));

	t_WCharsLoc wc_cur = wc_s;

	t_Complex a_arr[3];

	// move along w line to final destination

	const double dw = wc_dest.w.real() >= wc_s.w.real() ? 0.0025 : -0.0025;

	const int Nw = int(abs((wc_dest.w.real() - wc_s.w.real()) / dw));

	for (int i = 0; i < Nw; i++) {

		g_pStabSolver->searchWave(wc_cur, stab::t_LSCond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED),
			stab::t_TaskTreat::SPAT);

		//wxLogMessage(_T("wc_cur:%s"), &(wc_cur.to_wstr()[0]));

		if (i < 3) {

			a_arr[i] = wc_cur.a;

		}
		else {

			a_arr[0] = a_arr[1];
			a_arr[1] = a_arr[2];
			a_arr[2] = wc_cur.a;

			wc_cur.a = smat::interpolate_parab(0.0, a_arr[0], dw, a_arr[1], 2 * dw, a_arr[2], 3 * dw);

		}

		wc_cur.w += dw;
	}

	// move along b line with wixed w

	const double db = wc_dest.b.real() >= wc_s.b.real() ? 0.0025 : -0.0025;

	const int Nb = int(abs((wc_dest.b.real() - wc_s.b.real()) / db));

	for (int i = 0; i < Nb; i++) {

		g_pStabSolver->searchWave(wc_cur, stab::t_LSCond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED),
			stab::t_TaskTreat::SPAT);

		//wxLogMessage(_T("wc_cur:%s"), &(wc_cur.to_wstr()[0]));

		if (i < 3) {

			a_arr[i] = wc_cur.a;

		}
		else {

			a_arr[0] = a_arr[1];
			a_arr[1] = a_arr[2];
			a_arr[2] = wc_cur.a;

			wc_cur.a = smat::interpolate_parab(0.0, a_arr[0], db, a_arr[1], 2 * db, a_arr[2], 3 * db);

		}

		wc_cur.b += db;
	}

	wxLogMessage(_T("_get_wave:wchars after w&b iters:%s"), &(wc_cur.to_wstr()[0]));

	wc_dest = wc_cur;

	// TODO:final adjustment can easily break convergence...
	//g_pStabSolver->searchWave(wc_dest, stab::t_LSCond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED),
	//	stab::t_TaskTreat::SPAT);

	wxLogMessage(_T("_get_wave:wchars output:%s"), &(wc_dest.to_wstr()[0]));

}

void _get_lst_wpacket_from_dns() {

	t_SpctrInfo spd;

	t_MFProfInfo mfd;

	_get_spctr_info(spd, "spectrum_input/_fft_info.dat");

	_get_mfprof_info(mfd, "spectrum_input/_amp_funcs_bl_data.dat");

	double w_fft, b_fft;
	double w_lst, b_lst;

	const double dOm_fft = (spd.Om_fft_max - spd.Om_fft_min) / (spd.Nw - 1);
	const double dB_fft = (spd.B_fft_max - spd.B_fft_min) / (spd.Nb - 1);

	const double pi = acos(-1.0);

	const mf::t_GeomPoint& xyz = g_pStabDB->get_pave_pt(0).xyz;

	g_pStabSolver->setContext(xyz);

	t_WCharsLoc wc_cur;
	t_WCharsLoc wc_prev;

	wc_cur.set_treat(stab::t_TaskTreat::SPAT);

	wc_cur.set_scales(g_pStabSolver->get_stab_scales());

	std::vector<t_Complex> wall_data(5);

	std::ofstream ofstr("output/_wall_spectrum_dns2lst.dat", std::ios::app);

	ofstr << spd.Nw << "\t" << spd.Nb << "\n";

	double w_fft_min_crop = 4.77;
	double w_fft_max_crop = 238.0;

	double b_fft_min_crop = 0.0;
	double b_fft_max_crop = 300.0;

	for (int i = 0; i < spd.Nw; i++){

		bool need_full_search = true;

		for (int j = 0; j < spd.Nb; j++) {

			w_fft = spd.Om_fft_min + i*dOm_fft;
			b_fft = spd.B_fft_min + j*dB_fft;

			w_lst = -2.0*pi*mfd.Dels_gndim / mfd.Ue_abs*w_fft;
			b_lst = 2.0*pi*mfd.Dels_gndim*b_fft;

			// no need to calc products where spectrum is zero (speed up only)
			// this is defined by hand and depends on how spectrum looks!
			// to avoid this set all min crops by 0 and max crops by HUGE_VAL
			bool crop_by_w = (abs(w_fft) < w_fft_min_crop) || (abs(w_fft) > w_fft_max_crop);
			bool crop_by_b = (abs(b_fft) < b_fft_min_crop) || (abs(b_fft) > b_fft_max_crop);
			if (crop_by_w || crop_by_b) {

				wxLogMessage(_T("Spectrum cropped for Om_fft=%lf, B_fft=%lf"), w_fft, b_fft);

				for (int k = 0; k < 5; k++) wall_data[k] = t_Complex(0.0, 0.0);

				need_full_search = true;

			}
			else {

				wxLogMessage(_T("Calculating lst2dns for Om_fft=%lf, B_fft=%lf"), w_fft, b_fft);

				try {

					wc_cur.w = w_lst;
					wc_cur.b = b_lst;

					if (need_full_search) {
						_get_wchars(wc_cur, _get_wc_s_default(wc_cur));
						wc_prev = wc_cur;
					}
					else {
						_get_wchars(wc_cur, wc_prev);
						wc_prev = wc_cur;
					}

					need_full_search = false;

					// wc_cur is set, calc scaled amp function

					_get_amp_funcs_from_full_spectrum(i, j, spd, mfd);

					_do_calc_dns2lst_wall_disturb(wc_cur, wall_data);

				}
				catch (...) {

					for (int k = 0; k < 5; k++) wall_data[k] = t_Complex(0.0, 0.0);

					need_full_search = true;

				}

			}

			ofstr << i << "\t" << j << "\t" << w_fft << "\t" << b_fft << "\t";

			// output should be in global nondim rf
			double coef2gndim[5];

			// u,v,w was ndim by ue
			coef2gndim[0] = mfd.Ue_abs;
			coef2gndim[1] = mfd.Ue_abs;
			coef2gndim[2] = mfd.Ue_abs;
			// p was ndim by rhoe*ue*ue
			coef2gndim[3] = mfd.Rho_e*mfd.Ue_abs*mfd.Ue_abs;
			// t was ndim by Te
			coef2gndim[4] = mfd.T_e;

			for (int k = 0; k < 5; k++) {
				ofstr << coef2gndim[k]*wall_data[k].real() << "\t" 
					  << coef2gndim[k]*wall_data[k].imag() << "\t";
			}

			ofstr << "\n";
			ofstr.flush();

		}
	}
	ofstr.close();

}

void task::calc_scal_prod_particle_test() {

	// get lst part of dns spectrum
	// using biortho scalar product
	if (true){

		_get_lst_wpacket_from_dns();

	}
	// test _get_wchars for a particular spectral point
	if (false) {

		const mf::t_GeomPoint& xyz = g_pStabDB->get_pave_pt(0).xyz;

		g_pStabSolver->setContext(xyz);

		t_WCharsLoc wchars_A;

		wchars_A.a = t_Complex(0.0, 0.0);
		wchars_A.b = -2.00;
		wchars_A.w = t_Complex(0.05, 0.0);

		wchars_A.set_treat(stab::t_TaskTreat::SPAT);

		wchars_A.set_scales(g_pStabSolver->get_stab_scales());

		_get_wchars(wchars_A, _get_wc_s_default(wchars_A));
	}
	// test dns vs lst amplitude functions for a particular spectral point
	if (false){
		t_SpctrInfo spd;

		t_MFProfInfo mfd;

		_get_spctr_info(spd, "spectrum_input/_fft_info.dat");

		_get_mfprof_info(mfd, "spectrum_input/_amp_funcs_bl_data.dat");

		_get_amp_funcs_from_full_spectrum(43,63,spd, mfd);

		_do_calc_scal_prod_spctr_max_point();
	}

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