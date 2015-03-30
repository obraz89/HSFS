#include "stdafx.h"

#include "tasks.h"

#include "PluginsManager.h"

#include "common_data.h"

#include "ProfileStab.h"

#include "Log.h"

#include "io_helpers.h"

#include "solvers_glob.h"

#include "mpi.h"

#include "msg_pack.h"
// tmp
#include "tests.h"
#include <time.h>

using namespace hsstab;

//******************************************************************************
// spatial global search subroutines
//******************************************************************************

bool search_global_initial_wr_fixed_spat(const task::TTaskParams& task_params, double a_wr, t_WCharsLoc& ret_wave){

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

			bool ok = search_global_initial_wr_fixed_spat(task_params, cur_w, cur_max_wave);

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

			// after xyz write flag - is point ok or not
			// to simplify reading for retrace

			int ok = 1;

			fostr_max<<pid<<"\t"<<xyz.x()<<"\t"<<xyz.y()<<"\t"<<xyz.z()<<"\t"<<ok<<"\t"
				<<lw.a.real()<<"\t"<<lw.a.imag()<<"\t"<<lw.b.real()<<"\t"<<lw.b.imag()<<"\t"<<lw.w.real()<<"\t"<<lw.w.imag()<<"\t"
				<<ld.a.real()<<"\t"<<ld.a.imag()<<"\t"<<ld.b.real()<<"\t"<<ld.b.imag()<<"\t"<<ld.w.real()<<"\t"<<ld.w.imag()
				<<"\n";
		}else{
			int ok = 0;
			fostr_max<<pid<<"\t"<<xyz.x()<<"\t"<<xyz.y()<<"\t"<<xyz.z()<<ok<<"\t"<<" --- Failed to find wave"<<"\n";

		}

		fostr_max.flush();


	}

	return;

};
//******************************************************************************
// temporal global search subroutines
//******************************************************************************

bool search_global_initial_ar_fixed_time(const task::TTaskParams& task_params, double a_ar, t_WCharsLoc& ret_wave){

	// solvers Contexts must be already set here!

	stab::t_LSBase* const stab_solver = g_pStabSolver;

	stab::t_GSBase* const gs_solver = g_pGSSolverTime;

	const t_StabScales& stab_scales = stab_solver->get_stab_scales();

	const int n_bt=g_taskParams.N_b;

	double bt_min = g_taskParams.b_ndim_min;
	double bt_max = g_taskParams.b_ndim_max;

	double dbt = (bt_max - bt_min)/double(n_bt);

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



void task::search_max_instab_fixed_point_time(const task::TTaskParams& task_params){

	double a_min = g_taskParams.a_ndim_min;
	double a_max = g_taskParams.a_ndim_max;

	int Na = g_taskParams.N_a;

	double da = (a_max - a_min)/double(Na);

	std::vector<t_WCharsLoc> max_waves_time;

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

		max_waves_time.resize(0);max_waves_time.clear();

		mf::t_GeomPoint xyz = g_pStabDB->get_pave_pt(pid).xyz;

		g_pStabSolver->setContext(xyz);

		g_pGSSolverTime->setContext(xyz);

		for (int j=0; j<Na; j++){

			int perc_done = double(j)/double(Na)*100;
			wxLogMessage(_T("\n=============SearchMax Loc : %d perc done =============\n"), perc_done);

			double cur_a = a_min + da*j;

			bool ok = search_global_initial_ar_fixed_time(task_params, cur_a, cur_max_wave);

			if (ok) {
				max_waves_time.push_back(cur_max_wave);
			} else
			{continue;}

			// debug
			const t_WCharsLoc& lw = cur_max_wave;
			fostr_all<<xyz.x()<<"\t"<<xyz.y()<<"\t"<<xyz.z()<<"\t"
				<<lw.a.real()<<"\t"<<lw.a.imag()<<"\t"<<lw.b.real()<<"\t"<<lw.b.imag()<<"\t"<<lw.w.real()<<"\t"<<lw.w.imag()
				<<"\n";
			fostr_all.flush();

		}

		if (max_waves_time.size()>0){
			cur_max_wave = t_WCharsLoc::find_max_instab_time(max_waves_time);

			t_WCharsLocDim ret_wave = cur_max_wave.make_dim();

			const t_WCharsLoc& lw = cur_max_wave;
			const t_WCharsLocDim& ld = ret_wave;

			// after xyz write flag - is point ok or not
			// to simplify reading for retrace

			int ok = 1;

			fostr_max<<pid<<"\t"<<xyz.x()<<"\t"<<xyz.y()<<"\t"<<xyz.z()<<"\t"<<ok<<"\t"
				<<lw.a.real()<<"\t"<<lw.a.imag()<<"\t"<<lw.b.real()<<"\t"<<lw.b.imag()<<"\t"<<lw.w.real()<<"\t"<<lw.w.imag()<<"\t"
				<<ld.a.real()<<"\t"<<ld.a.imag()<<"\t"<<ld.b.real()<<"\t"<<ld.b.imag()<<"\t"<<ld.w.real()<<"\t"<<ld.w.imag()
				<<"\n";
		}else{
			int ok = 0;
			fostr_max<<pid<<"\t"<<xyz.x()<<"\t"<<xyz.y()<<"\t"<<xyz.z()<<ok<<"\t"<<" --- Failed to find wave"<<"\n";

		}

		fostr_max.flush();


	}

	return;
};


// serializer for a gs record

class t_MsgGSRec : public t_MsgBase{

public:
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

// testing mpi send-receive, remove when done

void generate_msg_pack_gs(t_MsgPack & msg_pack){

	// alias
	const task::TTaskParams& task_params = g_taskParams;

	double a_min = g_taskParams.a_ndim_min;
	double a_max = g_taskParams.a_ndim_max;

	int Na = g_taskParams.N_a;

	double da = (a_max - a_min)/double(Na);

	std::vector<t_WCharsLoc> max_waves_time;

	t_WCharsLoc cur_max_wave;

	int pave_pnt_ind = g_taskParams.pave_point_id;
	if ( pave_pnt_ind>=g_pStabDB->get_npoints())
	{
		wxLogError(_T("StabDb: point_id is out of range: point_id=%d"), pave_pnt_ind);
		return;
	}

	// debug
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
	wxLogMessage(_T("task params pid on worker %d : %d"), mpi_rank, g_taskParams.pave_point_id);

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

		max_waves_time.resize(0);max_waves_time.clear();

		mf::t_GeomPoint xyz = g_pStabDB->get_pave_pt(pid).xyz;

		g_pStabSolver->setContext(xyz);

		g_pGSSolverTime->setContext(xyz);

		for (int j=0; j<Na; j++){

			int perc_done = double(j)/double(Na)*100;
			wxLogMessage(_T("\n=============SearchMax Loc : %d perc done =============\n"), perc_done);

			double cur_a = a_min + da*j;

			bool ok = search_global_initial_ar_fixed_time(task_params, cur_a, cur_max_wave);

			if (ok) {
				max_waves_time.push_back(cur_max_wave);
			} else
			{continue;}

		}

		t_MsgGSRec a_msg(msg_pack.get_msg(pid));

		if (max_waves_time.size()>0){
			cur_max_wave = t_WCharsLoc::find_max_instab_time(max_waves_time);

			const int ok = 1;

			t_WCharsLocDim ret_wave = cur_max_wave.make_dim();

			const t_WCharsLoc& lw = cur_max_wave;
			const t_WCharsLocDim& ld = ret_wave;

			// fill msg_flat, init msg and push to msg pack
			a_msg[0] = g_pStabDB->get_global_pid(pid);
			a_msg[1] = xyz.x();
			a_msg[2] = xyz.y();
			a_msg[3] = xyz.z();
			a_msg[4] = ok;
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
			const int ok = 0; double false_val = -1.0;
			a_msg[0] = g_pStabDB->get_global_pid(pid);
			a_msg[1] = xyz.x();
			a_msg[2] = xyz.y();
			a_msg[3] = xyz.z();
			a_msg[4] = ok;
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

void task::mpi_test(){

	/*
	double w_min = g_taskParams.w_ndim_min;
	double w_max = g_taskParams.w_ndim_max;

	int Nw = g_taskParams.N_w;

	double dw = (w_max - w_min)/double(Nw);

	std::vector<t_WCharsLoc> max_waves_spat;

	t_WCharsLoc cur_max_wave;

	int mpi_rank, mpi_size;

	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

	int pid_s_glob, pid_e_glob, pid_n_glob;
	int pid_s_loc, pid_e_loc;

	pid_s_glob = 0;
	pid_n_glob = g_pStabDB->get_npoints();
	pid_e_glob = pid_n_glob - 1;

	wxLogMessage(_T("rank=%d, size=%d"), mpi_rank, mpi_size);

	const int nAvg   = pid_n_glob / mpi_size;
	const int nResdl = pid_n_glob % mpi_size;

	wxLogMessage(_T("avg=%d, resid=%d"), nAvg, nResdl);

	for( int r = 0; r < mpi_size; ++r )
	{
		int bs = r * nAvg + ((r<nResdl) ?r :nResdl);
		int be = bs + nAvg + ((r<nResdl) ?1 :0) - 1;

		pid_s_loc = bs;
		pid_e_loc = be;

		wxLogMessage( _("* MPI rank %d owns range : %d-%d"), mpi_rank, bs, be );

	}
	//==========================

	char fname[33];

	for (int pid=pid_s_loc; pid<=pid_e_loc; pid++){

		// do some job, collect on master and output...


	}

	return;
	*/

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

		std::wofstream ofstr(_T("output.txt"));


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

	gs_solver->writeSpectrum(_T("output/spectrum.dat"));

}
