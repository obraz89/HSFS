#include "stdafx.h"
#include "log.h"

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
		//test::king_m35_generate_pave_points();getchar();return 0;
		//test::king_m35_eN_time_envelope();

		//test::king_m35_eN_spat_fixedB();
		//return 0;

		task::init_glob_solvers();

		// tmp. debug
		//mf::t_GeomPoint test_xyz(0.08859,0.022,0.0);
		//g_pStabSolver->setContext(test_xyz);return 0;
		//g_pStabSolver->dumpProfileStab(_T("profiles_stab_test.dat"));

		//double ff = g_pMFDomain->calc_bl_thick(test_xyz);
		//std::cout<<"BL Thick:"<<ff<<"\n";

		//std::cout<<"Dels="<<g_pStabSolver->get_stab_scales().Dels<<";"
		//	     <<"UeDim="<<g_pStabSolver->get_stab_scales().UeDim
		//		 <<"Me"<<g_pStabSolver->get_stab_scales().Me<<"\n";

		//g_pStabSolver->setContext(mf::t_GeomPoint(0.458244,0.054810,0.0));
		//return 0;

		// tmp, debugging
		//test::gs_lapack_vs_petsc();return 0;

		task::init_stab_dbs();

		switch (g_taskParams.id)
		{
		case task::TTaskType::SearchInstabLoc:
			switch (g_taskParams.spattime)
			{
			case task::TSpatTime::Spat:
				task::search_max_instab_fixed_point_spat(g_taskParams);
				break;
			case task::TSpatTime::Time:
				task::search_max_instab_fixed_point_time(g_taskParams);
				break;
			default :
				ssuGENTHROW(_T("Search Max Instab Local: Wrong Mode"));
				break;
			}
			break;
		case task::TTaskType::Retrace:
			task::retrace_wplines_wfixed();
			break;
		case task::TTaskType::GetProfiles:
			task::get_profiles();
			break;
		case task::MPITest:
			task::mpi_test();
			break;
		default:
			break;
		}
	}
	catch(t_GenException e){
		wxLogError(e.what());
	}
	// todo - catch & dump other exceptions
	/*
	catch(...){
		wxLogError(_T("Something went wrong...see log"));
	}*/

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