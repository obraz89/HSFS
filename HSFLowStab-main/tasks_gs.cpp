#include "stdafx.h"

#include "tasks.h"

#include "PluginsManager.h"

#include "common_data.h"

#include "ProfileStab.h"

//#include "Log.h"

#include "wx/log.h"

#include "io_helpers.h"

#include "solvers_glob.h"

#include "mpi.h"

#include "msg_pack.h"
// tmp
#include "tests.h"
#include <time.h>

using namespace hsstab;
using namespace task;

//******************************************************************************
// spatial global search subroutines
//******************************************************************************

bool search_global_initial_wr_fixed_spat(double a_wr, t_WCharsLoc& ret_wave){

	// solvers Contexts must be already set here!

	stab::t_LSBase* const stab_solver = g_pStabSolver;

	stab::t_GSBase* const gs_solver = g_pGSSolverSpat;

	const t_StabScales& stab_scales = stab_solver->get_stab_scales();

	const int n_bt=g_taskParams.N_b;

	double bt_min = g_taskParams.b_ndim_min;
	double bt_max = g_taskParams.b_ndim_max;

	double dbt = (n_bt>1) ? (bt_max - bt_min)/double(n_bt-1) : 0.0;

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

		// tmp, dump ai vs br
		std::ofstream ofstr("output/ai_vs_br.dat");
		for (int i = 0; i < waves_spat.size(); i++) {
			ofstr << waves_spat[i].b.real() <<"\t"<< abs(waves_spat[i].a.imag())<<"\n";
		}
		ofstr.close();

		ret_wave = t_WCharsLoc::find_max_instab_spat(waves_spat); 
		return true;
	}else 
		return false;
};


void task::search_max_instab_fixed_point_spat(int pid, t_GSWaveInfo& winfo){

	double w_min = g_taskParams.w_ndim_min;
	double w_max = g_taskParams.w_ndim_max;

	int Nw = g_taskParams.N_w;

	double dw = (Nw>1) ? (w_max - w_min)/double(Nw-1) : 0.0;

	std::vector<t_WCharsLoc> max_waves_spat;

	t_WCharsLoc cur_max_wave;

	max_waves_spat.resize(0);max_waves_spat.clear();

	mf::t_GeomPoint xyz = g_pStabDB->get_pave_pt(pid).xyz;

	g_pStabSolver->setContext(xyz);

	g_pGSSolverSpat->setContext(xyz);

	for (int j=0; j<Nw; j++){

		int perc_done = double(j)/double(Nw)*100;
		wxLogMessage(_T("\n=============SearchMax Loc : %d perc done =============\n"), perc_done);

		double cur_w = w_min + dw*j;

		bool ok = search_global_initial_wr_fixed_spat(cur_w, cur_max_wave);

		if (ok) {
			max_waves_spat.push_back(cur_max_wave);
		} else
		{continue;}

		// debug
		/*
		const t_WCharsLoc& lw = cur_max_wave;
		fostr_all<<xyz.x()<<"\t"<<xyz.y()<<"\t"<<xyz.z()<<"\t"
			<<lw.a.real()<<"\t"<<lw.a.imag()<<"\t"<<lw.b.real()<<"\t"<<lw.b.imag()<<"\t"<<lw.w.real()<<"\t"<<lw.w.imag()
			<<"\n";
		fostr_all.flush();
		*/

	}

	if (max_waves_spat.size()>0){
		cur_max_wave = t_WCharsLoc::find_max_instab_spat(max_waves_spat);

		winfo.wave = cur_max_wave;
		t_WCharsLocDim wave_dim = cur_max_wave.make_dim();
		winfo.ok = true;

		const t_WCharsLoc& lw = winfo.wave;
		const t_WCharsLocDim& ld = wave_dim;

		/*
		// after xyz write flag - is point ok or not
		// to simplify reading for retrace

		fostr_max<<pid<<"\t"<<xyz.x()<<"\t"<<xyz.y()<<"\t"<<xyz.z()<<"\t"<<ok<<"\t"
			<<lw.a.real()<<"\t"<<lw.a.imag()<<"\t"<<lw.b.real()<<"\t"<<lw.b.imag()<<"\t"<<lw.w.real()<<"\t"<<lw.w.imag()<<"\t"
			<<ld.a.real()<<"\t"<<ld.a.imag()<<"\t"<<ld.b.real()<<"\t"<<ld.b.imag()<<"\t"<<ld.w.real()<<"\t"<<ld.w.imag()
			<<"\n";
		*/
	}else{
		winfo.ok = 0;
		/*fostr_max<<pid<<"\t"<<xyz.x()<<"\t"<<xyz.y()<<"\t"<<xyz.z()<<ok<<"\t"<<" --- Failed to find wave"<<"\n";*/
	}

	return;

};

void task::search_wchars_loc_wvec_fix() {
	
	int npave_pts = g_pStabDB->get_npoints();

	int pave_pnt_ind = g_taskParams.pave_point_id;
	if (pave_pnt_ind >= g_pStabDB->get_npoints())
	{
		wxLogError(_T("StabDb: point_id is out of range: point_id=%d"), pave_pnt_ind);
		return;
	}

	int pid_s, pid_e;
	if (pave_pnt_ind<0) {
		pid_s = 0;
		pid_e = g_pStabDB->get_npoints() - 1;
	}
	else {
		pid_s = pave_pnt_ind;
		pid_e = pave_pnt_ind;
	}

	wxLogMessage(_T("SearchWCharsLocWVecFix: point_id=%d"), pid_s);

	// for now, working only with first pid = pid_s
	for (int pid = pid_s; pid <= pid_s; pid++) {

		try {
			const mf::t_GeomPoint& test_xyz = g_pStabDB->get_pave_pt(pid).xyz;

			g_pStabSolver->setContext(test_xyz);

			t_WCharsLoc w_init; w_init.set_treat(stab::t_TaskTreat::SPAT);

			bool read_ok = read_max_wave_pid(pid, _T("wchars_max_loc.dat"), w_init);

			if (!read_ok)
				ssuGENTHROW(_T("Failed to read max wave chars, skipping wpline"));
			else
				wxLogMessage(_T("Max Wave pid=%d read from file: ok"), pid);

			t_WCharsLoc w_dest = w_init;

			wxLogMessage(_T("SearchWCharsLocWVecFix: using w_ndim_min as target freq value"));
			w_dest.w = g_taskParams.w_ndim_min;

			g_pStabSolver->searchWaveFixWVecDirSpat(w_init, w_dest);

		}
		catch (t_GenException e) {
			wxLogMessage(e.what());
		}
		catch (...) {
			wxLogMessage(_T("Failed to find wave, pid=%d"), pid);
		}

	};

}

void task::search_wchars_loc_wb_shift() {

	int npave_pts = g_pStabDB->get_npoints();

	int pave_pnt_ind = g_taskParams.pave_point_id;
	if (pave_pnt_ind >= g_pStabDB->get_npoints())
	{
		wxLogError(_T("StabDb: point_id is out of range: point_id=%d"), pave_pnt_ind);
		return;
	}

	int pid_s, pid_e;
	if (pave_pnt_ind<0) {
		pid_s = 0;
		pid_e = g_pStabDB->get_npoints() - 1;
	}
	else {
		pid_s = pave_pnt_ind;
		pid_e = pave_pnt_ind;
	}

	wxLogMessage(_T("SearchWCharsLocWBShift: point_id=%d"), pid_s);

	// for now, working only with first pid = pid_s
	for (int pid = pid_s; pid <= pid_s; pid++) {

		try {
			const mf::t_GeomPoint& test_xyz = g_pStabDB->get_pave_pt(pid).xyz;

			g_pStabSolver->setContext(test_xyz);

			t_WCharsLoc w_init; w_init.set_treat(stab::t_TaskTreat::SPAT);

			bool read_ok = read_max_wave_pid(pid, _T("wchars_max_loc.dat"), w_init);

			if (!read_ok)
				ssuGENTHROW(_T("Failed to read max wave chars, skipping wpline"));
			else
				wxLogMessage(_T("Max Wave pid=%d read from file: ok"), pid);

			t_WCharsLoc w_dest = w_init;

			wxLogMessage(_T("SearchWCharsLocWBShift: using w_ndim_min as target freq value"));
			w_dest.w = g_taskParams.w_ndim_min;

			wxLogMessage(_T("SearchWCharsLocWBShift: using b_ndim_min as target beta value"));
			w_dest.b = g_taskParams.b_ndim_min;

			wxLogMessage(_T("SearchWCharsLocWBShift: using Nw as number of iterations"));

			g_pStabSolver->searchWaveWBShift(w_init, w_dest, g_taskParams.N_w);

		}
		catch (t_GenException e) {
			wxLogMessage(e.what());
		}
		catch (...) {
			wxLogMessage(_T("Failed to find wave, pid=%d"), pid);
		}

	};

}

void write_wave_to_db(const t_GSWaveInfo& info){

	char fname[33];

	sprintf(fname, "%s/wchars_all_loc.dat", hsstab::OUTPUT_DIR.ToAscii());
	std::ofstream fostr_all(fname);

	sprintf(fname, "%s/wchars_max_loc.dat", hsstab::OUTPUT_DIR.ToAscii());
	std::ofstream fostr_max(fname);

}


//******************************************************************************
// temporal global search subroutines
//******************************************************************************

bool search_global_initial_ar_fixed_time(double a_ar, t_WCharsLoc& ret_wave){

	// solvers Contexts must be already set here!

	stab::t_LSBase* const stab_solver = g_pStabSolver;

	stab::t_GSBase* const gs_solver = g_pGSSolverTime;

	const t_StabScales& stab_scales = stab_solver->get_stab_scales();

	const int n_bt=g_taskParams.N_b;

	double bt_min = g_taskParams.b_ndim_min;
	double bt_max = g_taskParams.b_ndim_max;

	double dbt = (n_bt>1) ? (bt_max - bt_min)/double(n_bt-1) : 0.0;

	const double a = a_ar;

	std::vector<t_WCharsLoc> waves_time;

	for (int j=0; j<n_bt; j++){

		std::vector<t_WCharsLoc> init_waves_raw;
		std::vector<t_WCharsLoc> init_waves_filtered;

		std::cout<<"J="<<j<<"\n";

		t_WCharsLoc init_wave;

		init_wave.b = bt_min + dbt*j;
		init_wave.a = a;

		init_waves_raw = gs_solver->getInstabModes(init_wave);

		init_waves_filtered = g_pStabSolver->filter_gs_waves_time(init_waves_raw, 
			stab::t_LSCond(stab::t_LSCond::B_FIXED|stab::t_LSCond::A_FIXED));

		for (int k=0; k<init_waves_filtered.size(); k++)
			waves_time.push_back(init_waves_filtered[k]);	

	}	// ~loop over betas

	if (waves_time.size()>0){
		ret_wave = t_WCharsLoc::find_max_instab_time(waves_time); 
		return true;
	}else 
		return false;
};



void task::search_max_instab_fixed_point_time(int pid, t_GSWaveInfo& winfo){

	double a_min = g_taskParams.a_ndim_min;
	double a_max = g_taskParams.a_ndim_max;

	int Na = g_taskParams.N_a;

	double da = (Na>1) ? (a_max - a_min)/double(Na) : 0.0;

	std::vector<t_WCharsLoc> max_waves_time;

	t_WCharsLoc cur_max_wave;

	max_waves_time.resize(0);max_waves_time.clear();

	mf::t_GeomPoint xyz = g_pStabDB->get_pave_pt(pid).xyz;

	g_pStabSolver->setContext(xyz);

	g_pGSSolverTime->setContext(xyz);

	for (int j=0; j<Na; j++){

		int perc_done = double(j)/double(Na)*100;
		wxLogMessage(_T("\n=============SearchMax Loc : %d perc done =============\n"), perc_done);

		double cur_a = a_min + da*j;

		bool ok = search_global_initial_ar_fixed_time(cur_a, cur_max_wave);

		if (ok) {
			max_waves_time.push_back(cur_max_wave);
		} else
		{continue;}

	}

	if (max_waves_time.size()>0){

		cur_max_wave = t_WCharsLoc::find_max_instab_time(max_waves_time);

		winfo.wave = cur_max_wave;
		winfo.ok = 1;

	}else{
		winfo.ok = 0;
	}

	return;
};


// serializer for a gs record

class t_MsgGSRec : public t_MsgBase{

public:

static const int SerSize;
    
	int getSerSize() const{return SerSize;};

	t_MsgGSRec(double* cont);
	t_MsgGSRec(const t_MsgBase& );
	void write2file(std::wofstream& ofstr) const;

};

t_MsgGSRec::t_MsgGSRec(double *cont):t_MsgBase(cont){}

t_MsgGSRec::t_MsgGSRec(const t_MsgBase& bb):t_MsgBase(bb){}

const int t_MsgGSRec::SerSize = 17;

void t_MsgGSRec::write2file(std::wofstream& ofstr) const{

	int pid = _pCont[0];
	double x = _pCont[1];
	double y = _pCont[2];
	double z = _pCont[3];
	int ok = _pCont[4];

	double w_ar = _pCont[5];
	double w_ai = _pCont[6];
	double w_br = _pCont[7];
	double w_bi = _pCont[8];
	double w_wr = _pCont[9];
	double w_wi = _pCont[10];

	double d_ar = _pCont[11];
	double d_ai = _pCont[12];
	double d_br = _pCont[13];
	double d_bi = _pCont[14];
	double d_wr = _pCont[15];
	double d_wi = _pCont[16];


	ofstr<<pid<<"\t"<<x<<"\t"<<y<<"\t"<<z<<"\t"<<ok<<"\t"
		<<w_ar<<"\t"<<w_ai<<"\t"<<w_br<<"\t"<<w_bi<<"\t"<<w_wr<<"\t"<<w_wi<<"\t"
		<<d_ar<<"\t"<<d_ai<<"\t"<<d_br<<"\t"<<d_bi<<"\t"<<d_wr<<"\t"<<d_wi
		<<"\n";


};

// generate gs message for a given pid and insert it into msg pack

void search_max_instab(int pid, t_GSWaveInfo& winfo){

	switch (g_taskParams.spattime)
	{
	case task::Spat:
		search_max_instab_fixed_point_spat(pid, winfo);
		break;
	case task::Time:
		search_max_instab_fixed_point_time(pid, winfo);
		break;
	default:
		wxLogError(_T("Error: unknown spat-time option in search max instab"));
		break;
	}

};

void update_msg_pack(t_MsgPack& msg_pack, int pid){

	t_MsgGSRec a_msg(msg_pack.get_msg(pid));

	mf::t_GeomPoint xyz = g_pStabDB->get_pave_pt(pid).xyz;

	t_GSWaveInfo winfo;

	search_max_instab(pid, winfo);

	if(winfo.ok==1){

		const t_WCharsLoc& lw = winfo.wave;
		const t_WCharsLocDim ld = lw.make_dim();

		// fill msg_flat, init msg and push to msg pack
		a_msg[0] = g_pStabDB->get_global_pid(pid);
		a_msg[1] = xyz.x();
		a_msg[2] = xyz.y();
		a_msg[3] = xyz.z();
		a_msg[4] = winfo.ok;
		a_msg[5] = lw.a.real();
		a_msg[6] = lw.a.imag();
		a_msg[7] = lw.b.real();
		a_msg[8] = lw.b.imag();
		a_msg[9] = lw.w.real();
		a_msg[10]= lw.w.imag();
		a_msg[11]= ld.a.real();
		a_msg[12]= ld.a.imag();
		a_msg[13]= ld.b.real();
		a_msg[14]= ld.b.imag();
		a_msg[15]= ld.w.real();
		a_msg[16]= ld.w.imag();

	}else{
		double false_val = -1.0;
		a_msg[0] = g_pStabDB->get_global_pid(pid);
		a_msg[1] = xyz.x();
		a_msg[2] = xyz.y();
		a_msg[3] = xyz.z();
		a_msg[4] = winfo.ok;
		a_msg[5] = false_val;
		a_msg[6] = false_val;
		a_msg[7] = false_val;
		a_msg[8] = false_val;
		a_msg[9] = false_val;
		a_msg[10]= false_val;
		a_msg[11]= false_val;
		a_msg[12]= false_val;
		a_msg[13]= false_val;
		a_msg[14]= false_val;
		a_msg[15]= false_val;
		a_msg[16]= false_val;

	}

}

// testing mpi send-receive, remove when done

void generate_msg_pack_gs(t_MsgPack & msg_pack){

	int pave_pnt_ind = g_taskParams.pave_point_id;
	if ( pave_pnt_ind>=g_pStabDB->get_npoints())
	{
		wxLogError(_T("StabDb: point_id is out of range: point_id=%d"), pave_pnt_ind);
		return;
	};

	int pid_s, pid_e;
	if (pave_pnt_ind<0){
		pid_s = 0;
		pid_e = g_pStabDB->get_npoints()-1;
	}else{
		wxLogError(_T("Fix gs calc for single point - multiproc version"));
		pid_s = pave_pnt_ind;
		pid_e = pave_pnt_ind;
	}

	for (int pid=pid_s; pid<=pid_e; pid++){

		if (pid_s!=pid_e){

			int task_perc_done = double(pid-pid_s)/double(pid_e-pid_s)*100;
			wxLogMessage(_T("\n=============Task : %d perc done =============\n"), task_perc_done);

		}

		update_msg_pack(msg_pack, pid);

	}

}

void write_msg(double *msg, int size, std::ofstream& ofstr){

	for (int i=0; i<size; i++)
		ofstr<<msg[i]<<"\t";

	ofstr<<"\n";

}

void write_msg_pack_gs(const t_MsgPack& msg_pack, std::wofstream& ofstr){

	for (int i=0; i<msg_pack.get_n_msgs(); i++){
		t_MsgBase msg_base = msg_pack.get_msg(i);
		t_MsgGSRec msg_gs = msg_base;
		msg_gs.write2file(ofstr);
	}

}

void task::do_global_search(){

	int  numtasks, rank, len, rc; 
	char hostname[MPI_MAX_PROCESSOR_NAME];

	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Get_processor_name(hostname, &len);
	printf ("Number of tasks= %d My rank= %d Running on %s\n", numtasks,rank,hostname);

	const int MASTER = 0;

	// maximum number of points on a worker
	// all message packs should be that long
	int max_pnts_num;

	int npoints_loc = g_pStabDB->get_npoints();

	wxLogMessage(_T("loc msg_pack size:%d"), npoints_loc);

	MPI_Allreduce(&npoints_loc, &max_pnts_num, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	wxLogMessage(_T("max msg_pack size:%d"), max_pnts_num);

	t_MsgPack msg_pack_snd(t_MsgGSRec::SerSize, npoints_loc);
	int msg_pack_len_max = msg_pack_snd.get_flat_len();

	generate_msg_pack_gs(msg_pack_snd);

	wxLogMessage(_T("msg generated"));

	if (rank!=MASTER){
		MPI_Send(msg_pack_snd.get_cont(), msg_pack_len_max, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}else{

		t_MsgPack msg_pack_rcv(t_MsgGSRec::SerSize, max_pnts_num);

		MPI_Status status;

		std::wofstream ofstr("output/wchars_max_loc.dat");


//		calcs are done on master too, write them
//		write_msg(msg_snd, BUF_SIZE, ofstr);
		write_msg_pack_gs(msg_pack_snd, ofstr);ofstr.flush();

		
		for (int i=1; i<numtasks; i++){
			MPI_Recv(msg_pack_rcv.get_cont(), msg_pack_rcv.get_flat_len(), MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			wxLogMessage(_T("master received msg from worker %d, writing to output file now...\n"), i);

			write_msg_pack_gs(msg_pack_rcv, ofstr);

		}
		
	}

	MPI_Barrier(MPI_COMM_WORLD);

}

// single point test
// to compare with AVF data

// lst closure of saddle point
// dai/dbr=0
// implies that bi=0, i.e.
// growth of instability along inviscid streamline (see definition of local rf)
// use spatial approach to find zero of f(br) = dai/dbr
bool task::do_global_search_find_max(const int pid){

	const mf::t_GeomPoint& test_xyz = g_pStabDB->get_pave_pt(pid).xyz;

	g_pGSSolverSpat->setContext(test_xyz);

	g_pStabSolver->setContext(test_xyz);

	// TODO: for now GS is always spat, should be read from gs file later
	t_WCharsLoc w_init; w_init.set_treat(stab::t_TaskTreat::SPAT);

	bool read_ok = read_max_wave_pid(pid, _T("wchars_max_loc.dat"), w_init);

	bool ok = true;
	const int max_iters = 100;

	// TODO: emprics with d_rg, move to config when tested
	// not working with small d_arg (like 1.0e-07)

	wxLogMessage(_T("Warning: Check darg in do_global_search_find_max !"));

	const double d_arg = 1.0e-03;
	const double tol = 1.0e-06;

	t_WCharsLoc w_max = w_init;

	double fun, fun_deriv;

	for (int i = 0; i<max_iters; i++) {

		// do gs when making large newton iteration step
		g_pStabSolver->calcAiDbDerivs(w_max, fun, fun_deriv, d_arg);

		// check if we are converged
		if (abs(fun)<tol) {
			std::wcout<<_T("search dai_db=0 Converged\n")<<w_max<<std::endl;
			return true;
		};

		t_WCharsLoc base_wave = w_max;

		w_max.b = base_wave.b - fun / fun_deriv;
	};

	wxLogMessage(_T("Error: search max ai vs br - no convergence\n"));

	w_max = w_init;
	return false;


}

// tmp, debugging...

void test::gs_lapack_vs_petsc(){

	stab::t_LSBase* const stab_solver = g_pStabSolver;

	stab::t_GSBase* const gs_solver = g_pGSSolverSpat;

	mf::t_GeomPoint xyz(0.6309, 0.0, 0.0553);

	stab_solver->setContext(xyz);

	gs_solver->setContext(xyz);

	t_WCharsLoc wave;

	wave.a = 0.3;
	wave.b = 0.7;
	wave.w = 0.196;

	int n_iters = 1;

	time_t start_t, end_t;

	time(&start_t);

	for (int i=0; i<n_iters; i++){

			gs_solver->getInstabModes(wave);

	}

	time(&end_t); double dt = difftime(end_t, start_t);

	std::cout<<n_iters<<"iters: elapsed"<<dt<<"\n";

	gs_solver->writeSpectrum("output/spectrum.dat");

}
