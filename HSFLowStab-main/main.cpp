#include "stdafx.h"

#include "gen_exception.h"

// support console
#ifndef _USE_OLD_OSTREAMS
using namespace std;
#endif

#include "cmd_parse.h"

#include "tests.h"

#include "tasks.h"

#include "solvers_glob.h"

#include "mpi.h"

bool load_Settings_n_Plugins();

//debug


int main(int argc, char* argv[]){
	int err = 0;

	int rc = MPI_Init(&argc,&argv);
	if (rc != MPI_SUCCESS) {
		printf ("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

#ifdef wxUSE_UNICODE
	wxChar **wxArgv = new wxChar*[argc + 1];
	{
		int n = 0;
		for(n=0; n<argc; n++)
		{
			wxMB2WXbuf warg = wxConvertMB2WX(argv[n]);
			wxArgv[n] = wxStrdup(warg);
		}
		wxArgv[n] = NULL;
	}
#else
#define wxArgv argv
#endif // wxUSE_UNICODE/!wxUSE_UNICODE

	processCmdLine(argc, wxArgv);
	wxLogMessage(
		_("\t-\n\t- Run %s\n\t-\n"),
		wxDateTime::Now().Format(_T("%Y-%m-%d %H:%M:%S")).c_str()
		);

	try{
		load_Settings_n_Plugins();
	}catch(t_GenException e){
		wxLogMessage(e.what_detailed());
		return false;
	}
	catch(...){
		wxLogMessage(_T("Unhandled Exception in load_Setting_n_Plugins. Aborting..."));
		return false;
	}


	//test::selfsim_M45_second_mode();

	//test::selfsim_M3_first_mode();

	//test::selfsim_M2_CF();

	//test::transhyb_base_wartman();

	//test::test_smat();

	//test::test_conj_grad_min2D();

	
	try{
		//test::king_m35();
		//test::king_m35_eN();
		//test::king_m35_generate_pave_points();goto finish;
		//test::king_m35_eN_time_envelope();

		//test::king_m35_eN_spat_fixedB();
		//return 0;
		
		task::init_glob_solvers();

		task::init_stab_dbs();

		stab::t_WPRetraceMode retrace_mode = g_pWPLine->get_retrace_mode();

		switch (g_taskParams.id)
		{
		case task::TTaskType::SearchInstabLoc:
			task::do_global_search();
			break;
		case task::TTaskType::SearchMaxInstabLoc:
			// using fixed debug value of pid
			task::do_global_search_find_max(0);
			break;
		case task::TTaskType::Retrace:
			switch (g_taskParams.spattime){

			case (task::TSpatTime::Spat):
				task::retrace_wplines_cond_spat(retrace_mode);
				break;
			case (task::TSpatTime::Time):
				wxLogError(_T("Error: retrace with TIME approach disabled"));
				break;
			}
			break;
		case task::RetraceMPI:
			task::retrace_MPI(retrace_mode);
			break;
		case task::RetraceStreamlines:
			task::retrace_streamlines();
			break;
		case task::GetWallGridLine:
			task::get_wall_gridline();
			break;
		case task::PostProcRetrace:
			task::postproc_retrace();
			break;
		case task::TTaskType::GetProfiles:
			task::get_profiles();
			break;
		case task::MPITest:
			task::do_global_search();
			break;
		case task::GetAmplitudeFuncs:
			task::get_amplitude_funcs();
			break;
		case task::GetMFChars:
			task::calc_Cp_etc();
			break;

		case task::CalcScalProd_LST_test:
			//task::calc_scal_prod_particle_test();
			break;

		case task::CalcScalProd_DNS_test:
			task::calc_scal_prod_particle_test(); 
			break;

		case task::GetBLSpectrum:
			task::get_bl_spectrum();
			break;
		case task::CalcMeanFlowRecDerivs:
			task::calc_mean_flow_rec_derivs();
			break;

		case task::CalcNeutralCurve:
			task::calc_neutral_curve();
			break;
		case task::Test:
			task::test();
			break;
		default:
			wxString errMsg(_T("Error: provided task not implemented"));
			wxLogError(errMsg);
			break;
		}
	}
	catch(t_GenException e){
		wxLogError(e.what());
	}

	catch (...) {
		wxLogMessage(_T("Unhandled Exception while doing task. Aborting..."));
		return false;
	}
	// todo - catch & dump other exceptions
	/*
	catch(...){
		wxLogError(_T("Something went wrong...see log"));
	}*/

finish:

	task::destroy_glob_solvers();

	MPI_Finalize();
}
/*
int WINAPI WinMain(	HINSTANCE	hInstance,			// Instance
					HINSTANCE	hPrevInstance,		// Previous Instance
					LPSTR		lpCmdLine,			// Command Line Parameters
					int			nCmdShow)			// Window Show State
{
	// all info redirected to debug console
	// THIS IS NECESSARY BECAUSE 
	// PAUSE IN FORTRAN LIBS WON'T WORK
	// 

	RedirectIOToConsole();
	load_Settings_n_Plugins();

	//test::king_al_2_new();

	//test::transhyb_base_08();

	//test::itam_hz();

	//test::profile_compar();

	test::selfsim_M45_second_mode();

	return 0;
}
*/