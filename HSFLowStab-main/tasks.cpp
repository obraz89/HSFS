#include "stdafx.h"

#include "tasks.h"

#include "PluginsManager.h"

#include "common_data.h"

#include "ProfileStab.h"

#include "Log.h"

#include "io_helpers.h"

#include "solvers_glob.h"

using namespace hsstab;

// global references

mf::t_DomainBase* g_pMFDomain;

stab::t_LSBase* g_pStabSolver;

stab::t_GSBase* g_pGSSolverTime;

stab::t_GSBase* g_pGSSolverSpat;

stab::t_StabDBase* g_pStabDB;

// configure values

wxString task::TaskNames[TaskNum];
wxString task::SpatTimeNames[SpatTimeNum];


void task::init_glob_solvers(){

	TCapsMF& caps_mf = G_Plugins.get_caps_mf();
	TCapsLS& caps_ls = G_Plugins.get_caps_ls();
	TCapsGS& caps_gs = G_Plugins.get_caps_gs();
	TCapsWPTrack& caps_wp = G_Plugins.get_caps_wp();
	
	g_pMFDomain = caps_mf.create_domain();
	try
	{
		g_pMFDomain->init(G_Plugins.get_plugin(hsstab::plgMF));
	}
	catch (t_GenException e)
	{
		wxLogMessage(e.what());
		return;
	}

	g_pStabSolver = caps_ls.create_ls_solver(*g_pMFDomain);
	g_pStabSolver->init(G_Plugins.get_plugin(plgLS));

	g_pGSSolverTime = caps_gs.create_gs_solver(*g_pMFDomain, stab::t_TaskTreat::TIME);
	g_pGSSolverTime->init(G_Plugins.get_plugin(plgGS));

	g_pGSSolverSpat = caps_gs.create_gs_solver(*g_pMFDomain, stab::t_TaskTreat::SPAT);
	g_pGSSolverSpat->init(G_Plugins.get_plugin(plgGS));

}

void task::init_stab_db(){

	int NPavePnts;

	std::string points_fname = wx_to_stdstr(g_taskParams.pave_grd_fname);

	std::ifstream f_cin(&points_fname[0]);

	f_cin>>NPavePnts;

	std::vector<mf::t_GeomPoint> PaveGrd(NPavePnts);
	const int BufSize = 256;
	char line[BufSize];

	double x,y,z;

	for (int i=0; i<NPavePnts; i++){

		if (!f_cin.good()){
			wxString msg(_T("Error during reading of pave points"));
			wxLogMessage(msg);ssuGENTHROW(msg);
		}

		f_cin>>x>>y>>z;
	
		PaveGrd[i] = mf::t_GeomPoint(x,y,z);
		
	}

	g_pStabDB = new stab::t_StabDBase();
	g_pStabDB->init_pave_pts(PaveGrd);

}

void task::destroy_glob_solvers(){

delete g_pStabSolver;

delete g_pGSSolverTime;

delete g_pGSSolverSpat;

delete g_pStabDB;

delete g_pMFDomain;

}

bool search_global_initial_wr_fixed(const task::TTaskParams& task_params, double a_wr, t_WCharsLoc& ret_wave){

	// solvers Contexts must be already set here!

	stab::t_LSBase* const stab_solver = g_pStabSolver;

	stab::t_GSBase* const gs_solver = g_pGSSolverSpat;

	const t_StabScales& stab_scales = stab_solver->get_stab_scales();

	const int n_bt=g_taskParams.N_b;

	double bt_min = g_taskParams.b_ndim_min;
	double bt_max = g_taskParams.b_ndim_max;

	double dbt = (bt_max - bt_min)/double(n_bt);

	const double w = a_wr;

	std::vector<t_WCharsLoc> waves_spat;

	for (int j=0; j<n_bt; j++){

		std::vector<t_WCharsLoc> init_waves_raw;
		std::vector<t_WCharsLoc> init_waves_filtered;

		std::cout<<"J="<<j<<"\n";

		t_WCharsLoc init_wave;

		init_wave.b = bt_min + dbt*j;
		init_wave.w = w;

		init_waves_raw = gs_solver->getInstabModes(init_wave);

		init_waves_filtered = g_pStabSolver->filter_gs_waves_spat(init_waves_raw, 
			stab::t_LSCond(stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED));

		for (int k=0; k<init_waves_filtered.size(); k++)
			waves_spat.push_back(init_waves_filtered[k]);	

	}	// ~loop over betas

	if (waves_spat.size()>0){
		ret_wave = t_WCharsLoc::find_max_instab_spat(waves_spat); 
		return true;
	}else 
		return false;
};


void task::search_max_instab_fixed_point_spat(const task::TTaskParams& task_params){

	double w_min = g_taskParams.w_ndim_min;
	double w_max = g_taskParams.w_ndim_max;

	int Nw = g_taskParams.N_w;

	double dw = (w_max - w_min)/double(Nw);

	std::vector<t_WCharsLoc> max_waves_spat;

	t_WCharsLoc cur_max_wave;

	int pave_pnt_ind = g_taskParams.pave_point_id;
	if ( pave_pnt_ind>=g_pStabDB->get_npoints())
	{
		wxLogError(_T("StabDb: point_id is out of range: point_id=%d"), pave_pnt_ind);
		return;
	}

	int pid_s, pid_e;
	if (pave_pnt_ind<0){
		pid_s = 0;
		pid_e = g_pStabDB->get_npoints()-1;
	}else{
		pid_s = pave_pnt_ind;
		pid_e = pave_pnt_ind;
	}

	char fname[33];

	sprintf(fname, "%s/wchars_all_loc.dat", hsstab::OUTPUT_DIR.ToAscii());
	std::ofstream fostr_all(fname);

	sprintf(fname, "%s/wchars_max_loc.dat", hsstab::OUTPUT_DIR.ToAscii());
	std::ofstream fostr_max(fname);

	for (int pid=pid_s; pid<=pid_e; pid++){

		if (pid_s!=pid_e){

			int task_perc_done = double(pid-pid_s)/double(pid_e-pid_s)*100;
			wxLogMessage(_T("\n=============Task : %d perc done =============\n"), task_perc_done);

		}

		max_waves_spat.resize(0);max_waves_spat.clear();

		mf::t_GeomPoint xyz = g_pStabDB->get_pave_pt(pid).xyz;

		g_pStabSolver->setContext(xyz);

		g_pGSSolverSpat->setContext(xyz);

		for (int j=0; j<Nw; j++){

			int perc_done = double(j)/double(Nw)*100;
			wxLogMessage(_T("\n=============SearchMax Loc : %d perc done =============\n"), perc_done);

			double cur_w = w_min + dw*j;

			bool ok = search_global_initial_wr_fixed(task_params, cur_w, cur_max_wave);

			if (ok) {
				max_waves_spat.push_back(cur_max_wave);
			} else
				{continue;}

			// debug
			const t_WCharsLoc& lw = cur_max_wave;
			fostr_all<<xyz.x()<<"\t"<<xyz.y()<<"\t"<<xyz.z()<<"\t"
				<<lw.a.real()<<"\t"<<lw.a.imag()<<"\t"<<lw.b.real()<<"\t"<<lw.b.imag()<<"\t"<<lw.w.real()<<"\t"<<lw.w.imag()
				<<"\n";
			fostr_all.flush();

		}

		if (max_waves_spat.size()>0){
			cur_max_wave = t_WCharsLoc::find_max_instab_spat(max_waves_spat);

			t_WCharsLocDim ret_wave = cur_max_wave.make_dim();

			const t_WCharsLoc& lw = cur_max_wave;
			const t_WCharsLocDim& ld = ret_wave;

			fostr_max<<pid<<"\t"<<xyz.x()<<"\t"<<xyz.y()<<"\t"<<xyz.z()<<"\t"
				<<lw.a.real()<<"\t"<<lw.a.imag()<<"\t"<<lw.b.real()<<"\t"<<lw.b.imag()<<"\t"<<lw.w.real()<<"\t"<<lw.w.imag()<<"\t"
				<<ld.a.real()<<"\t"<<ld.a.imag()<<"\t"<<ld.b.real()<<"\t"<<ld.b.imag()<<"\t"<<ld.w.real()<<"\t"<<ld.w.imag()
				<<"\n";
		}else{
			fostr_max<<pid<<"\t"<<xyz.x()<<"\t"<<xyz.y()<<"\t"<<xyz.z()<<" --- Failed to find wave"<<"\n";

		}

		fostr_max.flush();


	}

	return;

};

void task::analyze_wchars(const std::string& fname){

	std::vector<t_WaveChars> wchars_loc_v;
	std::vector<t_WaveChars> wchars_loc_dim_v;

	std::ifstream fistr("wchars_max_loc.dat");

	std::string line;line.resize(256);

	mf::t_GeomPoint xyz; t_WaveChars lw; t_WaveChars ld;
	double ar, ai, br, bi, wr, wi;
	double adr, adi, bdr, bdi, wdr, wdi;
	double x, y, z;
	while(true){
		if (fistr.getline(&line[0], 256)){

			std::istringstream istr(line);
			istr>>x>>y>>z
				>>ar>>ai>>br>>bi>>wr>>wi
				>>adr>>adi>>bdr>>bdi>>wdr>>wdi;

			//xyz.set(x,y,z);

			lw.a.real(ar);lw.a.imag(ai);
			lw.b.real(br);lw.b.imag(bi);
			lw.w.real(wr);lw.w.imag(wi);

			ld.a.real(adr);ld.a.imag(adi);
			ld.b.real(bdr);ld.b.imag(bdi);
			ld.w.real(wdr);ld.w.imag(wdi);

			wchars_loc_v.push_back(lw);
			wchars_loc_dim_v.push_back(ld);

		} else break;
	}

	double an_min, an_max, bn_min, bn_max, wn_min, wn_max, cn_min, cn_max;
	double ad_min, ad_max, bd_min, bd_max, wd_min, wd_max, cd_min, cd_max;

	an_min = bn_min = wn_min = cn_min = ad_min = bd_min = wd_min = cd_min = 1.0E+12;
	an_max = bn_max = wn_max = cn_max = ad_max = bd_max = wd_max = cd_max = -1.0E+12;

	for (int i=0; i<wchars_loc_v.size(); i++){

		const t_WCharsLoc& wl = wchars_loc_v[i];

		if (wl.a.real()<an_min) an_min = wl.a.real();
		if (wl.a.real()>an_max) an_max = wl.a.real();

		if (wl.b.real()<bn_min) bn_min = wl.b.real();
		if (wl.b.real()>bn_max) bn_max = wl.b.real();

		if (wl.w.real()<wn_min) wn_min = wl.w.real();
		if (wl.w.real()>wn_max) wn_max = wl.w.real();

	}

	for (int i=0; i<wchars_loc_dim_v.size(); i++){

		const t_WaveChars& wd = wchars_loc_dim_v[i];

		if (wd.a.real()<ad_min) ad_min = wd.a.real();
		if (wd.a.real()>ad_max) ad_max = wd.a.real();

		if (wd.b.real()<bd_min) bd_min = wd.b.real();
		if (wd.b.real()>bd_max) bd_max = wd.b.real();

		if (wd.w.real()<wd_min) wd_min = wd.w.real();
		if (wd.w.real()>wd_max) wd_max = wd.w.real();

	}

	std::ofstream out_str("wchars_range.dat");

	out_str<<"an=["<<an_min<<";"<<an_max<<"]\n"
		   <<"bn=["<<bn_min<<";"<<bn_max<<"]\n"
		   <<"wn=["<<wn_min<<";"<<wn_max<<"]\n\n"
		   <<"ad=["<<ad_min<<";"<<ad_max<<"]\n"
		   <<"bd=["<<bd_min<<";"<<bd_max<<"]\n"
		   <<"wd=["<<wd_min<<";"<<wd_max<<"]\n";


}

void task::search_max_instab_fixed_point_time(const task::TTaskParams& task_params){
	wxLogError(_T("Migrate from tests!"));
};


void do_retrace_wplines_wfixed_bfree(const task::TTaskParams& params){
	wxLogMessage(_T("Migrate do_retrace_wplines_wfixed_bfree from tests"));
}

bool read_max_wave_pid(int pid, const std::wstring& fname_max_waves, t_WCharsLoc& wave){

	std::wifstream ifstr(&fname_max_waves[0]);

	bool line_read_ok = true; const int MaxBufSize = 256;
	wchar_t line[MaxBufSize];

	int cur_pid; double x,y,z,ar,ai,br,bi,wr,wi;

	bool rec_found = false;
	do 
	{
		line_read_ok = ifstr.getline(&line[0], MaxBufSize); 

		std::wistringstream istr(line);
		istr>>cur_pid>>x>>y>>z
			>>ar>>ai>>br>>bi>>wr>>wi;
		if (pid==cur_pid){

			wave.a = t_Complex(ar, ai);
			wave.b = t_Complex(br, bi);
			wave.w = t_Complex(wr, wi);
			rec_found = true;
			break;

		}

	} while (line_read_ok);

	return rec_found;
};


void do_retrace_wplines_wfixed_bfixed(const task::TTaskParams& params){

	wxChar szFname[64];

	swprintf(szFname, _T("%s/Wave_pack_lines_wbfixed.dat.dat"),hsstab::OUTPUT_DIR.c_str());
	std::wstring fout_wplines_path(szFname);

	swprintf(szFname, _T("%s/max_N_wbfixed.dat"),hsstab::OUTPUT_DIR.c_str());
	std::wstring fout_maxnfactor_path(szFname);


	int npave_pts = g_pStabDB->get_npoints();

	TCapsWPTrack& caps_wp = G_Plugins.get_caps_wp();
	stab::t_WPTrackBase* wp_line = caps_wp.create_wp_track(*g_pMFDomain);
	wp_line->init(G_Plugins.get_plugin(plgWPTrack));

	//double bs = g_taskParams.b_ndim_min;
	//double be = g_taskParams.b_ndim_max;
	//int Nb = g_taskParams.N_b;

	//double ws = g_taskParams.w_ndim_min;
	//double we = g_taskParams.w_ndim_max;
	//int Nw = g_taskParams.N_w;

	for (int pid=0; pid<npave_pts; pid++){

		//for (int i_b=0; i_b<g_taskParams.N_b; i_b++)
		//for (int i_w=0; i_w<g_taskParams.N_w; i_w++)

		try{

			const mf::t_GeomPoint& test_xyz = g_pStabDB->get_pave_pt(pid).xyz;

			g_pGSSolverSpat->setContext(test_xyz);

			g_pStabSolver->setContext(test_xyz);

			t_WCharsLoc w_init;w_init.set_treat(stab::t_TaskTreat::SPAT);

			bool read_ok = read_max_wave_pid(pid, _T("wchars_max_loc.dat"), w_init);

			if (!read_ok) 
				ssuGENTHROW(_T("Failed to read max wave chars, skipping wpline"));
			else
				wxLogMessage(_T("Max Wave pid=%d read from file: ok"), pid);

			stab::t_LSCond cond(stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED);
			g_pStabSolver->searchWave(w_init, cond, stab::t_TaskTreat::SPAT);

			w_init.set_scales(g_pStabSolver->get_stab_scales());
/*
			double cur_b = bs + (be-bs)/double(Nb)*i_b;
			double cur_w = ws + (we-ws)/double(Nw)*i_w;

			w_init.b = cur_b;
			w_init.w = cur_w;

			std::vector<t_WCharsLoc> gs_waves = g_pGSSolverSpat->getInstabModes(w_init);
			std::vector<t_WCharsLoc> waves_filtered = g_pStabSolver->filter_gs_waves_spat(gs_waves, 
				stab::t_LSCond(stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED));

			w_init = t_WCharsLoc::find_max_instab_spat(waves_filtered);
*/

			// TODO: WPLine ids
			if (w_init.a.imag()<0.0){

				int perc_complete = double(pid)/double(npave_pts)*100.;
				wxLogMessage(_T("start retrace WPLine for point = %d, completed %d perc"), pid, perc_complete);

				std::wcout<<_T("Init for wpline:")<<w_init; 

				wp_line->retrace(test_xyz, w_init, *g_pStabSolver, *g_pGSSolverSpat, stab::t_WPRetraceMode::WB_FIXED);
				wp_line->print_to_file(fout_wplines_path, std::ios::app);

				g_pStabDB->update(*wp_line);


			}else{

				wxLogMessage(_T("No instabs found, skipping wpline"));

			}

			}catch(const t_GenException& ex){
				wxLogMessage(ex.what());
			}
			catch(...){

				wxLogMessage(_T("Failed to retrace WPLIneId :Point_id=%d"), pid);
			}
		}

		//StabDB.to_cone_ref_frame(HALF_CONE_ANGLE);
		g_pStabDB->export(fout_maxnfactor_path);

		delete wp_line;
		return;

}

void task::retrace_wplines_wfixed(){

	switch (g_taskParams.retrace_mode)
	{
	case 0:
		do_retrace_wplines_wfixed_bfree(g_taskParams);
		break;
	case 1:
		do_retrace_wplines_wfixed_bfixed(g_taskParams);
		break;
	default:
		wxLogError(_T("Error: Wrong mode for retrace"));
	}

};

void task::get_profiles(){

	int npts = g_pStabDB->get_npoints();

	if (npts<=0) wxLogError(_T("Err: No points provided to extract profiles from..."));

	for (int j=0; j<npts; j++){

		mf::t_GeomPoint xyz = g_pStabDB->get_pave_pt(j).xyz;

		g_pStabSolver->setContext(xyz);

		wchar_t szFname[33];

		t_ProfileNS prof_NS(*g_pMFDomain);

		mf::t_ProfDataCfg data_cfg;
		data_cfg.ThickCoef = g_pMFDomain->get_prof_extr_cfg().ThickCoefDefault;
		prof_NS.initialize(xyz, data_cfg, blp::t_NSInit::EXTRACT);

		swprintf(szFname, _T("%s/ProfileNS_%d.dat"),hsstab::OUTPUT_DIR.c_str(), j);
		prof_NS.dump(szFname);

		swprintf(szFname, _T("%s/ProfileStab_%d.dat"),hsstab::OUTPUT_DIR.c_str(), j);
		g_pStabSolver->dumpProfileStab(szFname);

		g_pMFDomain->dump_full_enthalpy_profile(xyz, j);

	}

}