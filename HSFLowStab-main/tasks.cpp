#include "stdafx.h"

#include "tasks.h"

#include "PluginsManager.h"

#include "common_data.h"

#include "ProfileStab.h"

//#include "Log.h"

#include "io_helpers.h"

#include "solvers_glob.h"

#include "mpi.h"

#include "fstream"

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


static const int MAX_FNAME_LEN = 100;

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

void task::init_stab_dbs(){

	int NPntsGlob;

	std::wifstream f_cin(g_taskParams.pave_grd_fname.ToAscii());

	f_cin>>NPntsGlob;

	if (NPntsGlob<=0 || NPntsGlob>1000000)
		wxLogMessage(_T("Error reading pave points : bad global number of points : %d"), NPntsGlob);

	std::vector<mf::t_GeomPoint> PaveGrdGlob(NPntsGlob);

	const int BufSize = 256;
	wxChar line[BufSize];

	double x,y,z;

	for (int i=0; i<NPntsGlob; i++){

		if (!f_cin.good()){
			wxString msg(_T("Error during reading of pave points"));
			wxLogMessage(msg);ssuGENTHROW(msg);
		}

		f_cin>>x>>y>>z;
	
		PaveGrdGlob[i] = mf::t_GeomPoint(x,y,z);
		
	}

	//=========================================

	// rank of worker, total number of workers 
	int mpi_rank, mpi_size;

	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

	const int nAvg   = NPntsGlob / mpi_size;
	const int nResdl = NPntsGlob % mpi_size;

	int r = mpi_rank;

	// interval of global pave points
	// to work with on current worker

	int bs = r * nAvg + ((r<nResdl) ?r :nResdl);
	int be = bs + nAvg + ((r<nResdl) ?1 :0) - 1;

	wxLogMessage( _("* MPI rank %d owns range : %d-%d"), mpi_rank, bs, be );

	g_pStabDB = new stab::t_StabDBase();
	g_pStabDB->init_pave_pts(PaveGrdGlob, bs, be);

}

void task::destroy_glob_solvers(){

delete g_pStabSolver;

delete g_pGSSolverTime;

delete g_pGSSolverSpat;

delete g_pStabDB;

delete g_pMFDomain;

}


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

bool read_max_wave_pid(int pid, const std::wstring& fname_max_waves, t_WCharsLoc& wave){


	wxChar szFname[64];

	swprintf(szFname, MAX_FNAME_LEN, _T("%s/%s"),hsstab::OUTPUT_DIR.c_str(), &fname_max_waves[0]);
	
	wxString szFname_wx(szFname);

	std::wifstream ifstr(szFname_wx.ToAscii());

	if (!ifstr.is_open())
	{
		// Dump the contents of the file to cout.
		wxLogError(_("Error: Can't find file with max wchars!"));
		return false;
	}


	bool line_read_ok = true; const int MaxBufSize = 256;
	wchar_t line[MaxBufSize];

	int cur_pid; double x,y,z,ar,ai,br,bi,wr,wi;

	do 
	{
		line_read_ok = ifstr.getline(&line[0], MaxBufSize); 

		std::wistringstream istr(line);

		int good_rec;

		istr>>cur_pid>>x>>y>>z>>good_rec
			>>ar>>ai>>br>>bi>>wr>>wi;
		if (pid==cur_pid){

			if (good_rec){

				wave.a = t_Complex(ar, ai);
				wave.b = t_Complex(br, bi);
				wave.w = t_Complex(wr, wi);
				return true;

			}else{
				return false;
			}

		}

	} while (line_read_ok);

	return false;
};

void retrace_wplines_cond(stab::t_WPRetraceMode a_mode_retrace){

	wxChar szFname[64];

	swprintf(szFname, MAX_FNAME_LEN, _T("%s/Wave_pack_lines_mode%d.dat.dat"),
		hsstab::OUTPUT_DIR.c_str(), g_taskParams.retrace_mode);

	std::wstring fout_wplines_path(szFname);

	swprintf(szFname, MAX_FNAME_LEN, _T("%s/max_N_mode%d.dat"),
		hsstab::OUTPUT_DIR.c_str(), g_taskParams.retrace_mode);

	std::wstring fout_maxnfactor_path(szFname);

	TCapsWPTrack& caps_wp = G_Plugins.get_caps_wp();
	stab::t_WPTrackBase* wp_line = caps_wp.create_wp_track(*g_pMFDomain);
	wp_line->init(G_Plugins.get_plugin(plgWPTrack));

	int pid_s, pid_e;

	if (g_taskParams.pave_point_id>=0){

		pid_s=g_taskParams.pave_point_id;

		pid_e=pid_s;

	}else{

		pid_s=0;

		pid_e=g_pStabDB->get_npoints()-1;

	};

	for (int pid=pid_s; pid<=pid_e; pid++){

		try{

			const mf::t_GeomPoint& test_xyz = g_pStabDB->get_pave_pt(pid).xyz;

			g_pGSSolverSpat->setContext(test_xyz);

			g_pStabSolver->setContext(test_xyz);

			// TODO: for now GS is always spat, should be read from gs file later
			t_WCharsLoc w_init;w_init.set_treat(stab::t_TaskTreat::SPAT);

			bool read_ok = read_max_wave_pid(pid, _T("wchars_max_loc.dat"), w_init);

			if (!read_ok) 
				ssuGENTHROW(_T("Failed to read max wave chars, skipping wpline"));
			else
				wxLogMessage(_T("Max Wave pid=%d read from file: ok"), pid);


			// TODO: pass w_init as is to retrace!
			//stab::t_LSCond cond(stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED)
			//g_pStabSolver->searchWave(w_init, cond, stab::t_TaskTreat::SPAT);
			//g_pStabSolver->calcGroupVelocity(w_init);

			w_init.set_scales(g_pStabSolver->get_stab_scales());

			// TODO: WPLine ids
			if (w_init.a.imag()<0.0){

				int perc_complete = double(pid)/double(pid_e-pid_s+1)*100.;
				wxLogMessage(_T("start retrace WPLine for point = %d, completed %d perc"), pid, perc_complete);

				std::wcout<<_T("Init for wpline:")<<w_init; 

				wp_line->retrace(test_xyz, w_init, *g_pStabSolver, *g_pGSSolverSpat, a_mode_retrace);
				
				wxString fout_wplines_path_wx(fout_wplines_path);
				wp_line->print_to_file(std::string(fout_wplines_path_wx.ToAscii()), std::ios::app);

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
		
		wxString fout_maxnfactor_path_wx(fout_maxnfactor_path);
		g_pStabDB->write_to_file(std::string(fout_maxnfactor_path_wx.ToAscii()));

		delete wp_line;
		return;


}


void task::retrace_wplines_wfixed_bfixed(){
	
	//stab::t_LSCond cond(stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED);
	//retrace_wplines_cond(cond, stab::t_WPRetraceMode::WB_FIXED);
	retrace_wplines_cond(stab::t_WPRetraceMode::WB_FIXED);
}

void task::retrace_wplines_wfixed_bfree(){

	//stab::t_LSCond cond(stab::t_LSCond::W_FIXED);
	//retrace_wplines_cond(cond, stab::t_WPRetraceMode::W_FIXED);
	retrace_wplines_cond(stab::t_WPRetraceMode::W_FIXED);
}

void task::get_amplitude_funcs(){

	int npave_pts = g_pStabDB->get_npoints();

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

	for (int pid=pid_s; pid<=pid_e; pid++){

		try{
			const mf::t_GeomPoint& test_xyz = g_pStabDB->get_pave_pt(pid).xyz;

			g_pStabSolver->setContext(test_xyz);

			t_WCharsLoc w_init;w_init.set_treat(stab::t_TaskTreat::SPAT);

			bool read_ok = read_max_wave_pid(pid, _T("wchars_max_loc.dat"), w_init);

			if (!read_ok) 
				ssuGENTHROW(_T("Failed to read max wave chars, skipping wpline"));
			else
				wxLogMessage(_T("Max Wave pid=%d read from file: ok"), pid);

			wxLogMessage(_T("Warning: using default mode for amplitude funcs reconstruction"));

			stab::t_LSCond cond;
			stab::t_TaskTreat stab_treat;
			switch (g_taskParams.spattime)
			{
			case task::Spat:

				cond.set(stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED, w_init);

				stab_treat = stab::t_TaskTreat::SPAT;

				break;
			case task::Time:

				cond.set(stab::t_LSCond::A_FIXED|stab::t_LSCond::B_FIXED, w_init);

				stab_treat = stab::t_TaskTreat::TIME;

				break;
			default:
				return;
			}

			char strFname[33];

			sprintf(strFname, "output/amp_funcs_%d.dat", pid);

			g_pStabSolver->searchWave(w_init, cond, stab_treat);

			g_pStabSolver->dumpEigenFuctions(strFname);

		}
		catch(t_GenException e){
			wxLogMessage(e.what());
		}
		catch(...){
			wxLogMessage(_T("Failed to reconstruct amp funcs for pid=%d"), pid);
		}

	};

}

void task::get_profiles(){

	int npts = g_pStabDB->get_npoints();

	if (npts<=0) wxLogError(_T("Err: No points provided to extract profiles from..."));

	for (int j=0; j<npts; j++){

		mf::t_GeomPoint xyz = g_pStabDB->get_pave_pt(j).xyz;

		g_pStabSolver->setContext(xyz);

		wchar_t szFname[33];

		// generate profile NS
		t_ProfileNS prof_NS(*g_pMFDomain);

		mf::t_ProfDataCfg data_cfg;
		data_cfg.ThickCoef = g_pMFDomain->get_prof_extr_cfg().ThickCoefDefault;
		prof_NS.initialize(xyz, data_cfg, blp::t_NSInit::EXTRACT);

		swprintf(szFname, MAX_FNAME_LEN, _T("%s/ProfileNS_%d.dat"),hsstab::OUTPUT_DIR.c_str(), j);
		
		wxString szFnameNS_wx(szFname);
		prof_NS.dump(std::string(szFnameNS_wx.ToAscii()));

		// generate profile MF (Glob RF)
		t_ProfMFGlob prof_MF(*g_pMFDomain);
		prof_MF.initialize(xyz, data_cfg, blp::t_NSInit::EXTRACT);

		swprintf(szFname, MAX_FNAME_LEN, _T("%s/ProfileMFGlob_%d.dat"),hsstab::OUTPUT_DIR.c_str(), j);
		
		wxString szFnameMF_wx(szFname);
		prof_MF.dump(std::string(szFnameMF_wx.ToAscii()));

		// generate profile Stab
		swprintf(szFname, MAX_FNAME_LEN, _T("%s/ProfileStab_%d.dat"),hsstab::OUTPUT_DIR.c_str(), j);
		
		wxString szFnameST_wx(szFname);
		g_pStabSolver->dumpProfileStab(std::string(szFnameST_wx.ToAscii()));

		g_pMFDomain->dump_full_enthalpy_profile(xyz, j);

	}

}