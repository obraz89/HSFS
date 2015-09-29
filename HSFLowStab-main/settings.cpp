///////////////////////////////////////////////////////////////////////////////
// Name:        main_settings.cpp
// Purpose:     Common parameters for HSFlow application
// Author:      Andrey V. Novikov
// Modified by: A. Obraz
///////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#pragma hdrstop

#include <wx/fileconf.h>

#include "common_data.h"
#include "PluginsManager.h"

#include "settings.h"
//-----------------------------------------------------------------------------


task::TTaskParams g_taskParams;

using namespace hsstab;

int get_task_id(const wxString& strName){

	for (int i=0; i<task::TaskNum; i++){
		if (strName==task::TaskNames[i])
		{
			return i;
		}
		
	}

	wxLogMessage(_T("Wrong task name, allowed options:"));
	for (int i=0; i<task::TaskNum; i++){
		wxLogMessage(task::TaskNames[i]);
	}
	ssuGENTHROW(_T("Wrong task name"));
	return -1;

}

int get_spattime_id(const wxString& strName){

	for (int i=0; i<task::SpatTimeNum; i++){
		if (strName==task::SpatTimeNames[i])
		{
			return i;
		}

	}
	ssuGENTHROW(_T("Wrong Stab Approach Name: Spat or Time allowed"));
	return -1;

}

bool load_Settings_n_Plugins()
{
	//
	// Load config
	//
	if( ! wxFileName::DirExists(CASE_SETTINGS_DIR) )
		if( ! wxFileName::Mkdir(CASE_SETTINGS_DIR, 0755) )
		{
			wxLogError(_("Can't create case settings dir '%s'"), hsstab::CASE_SETTINGS_DIR.c_str());
			exit(EXIT_FAILURE);
		}

		if( ! wxFileName::DirExists(OUTPUT_DIR) )
			if( ! wxFileName::Mkdir(OUTPUT_DIR, 0755) )
			{
				wxLogError(_("Can't create output dir '%s'"), hsstab::OUTPUT_DIR.c_str());
				exit(EXIT_FAILURE);
			}

	wxFileConfig* conf = new wxFileConfig(
		_T("HSStab Main"), _T("Obraz"),
		CASE_SETTINGS_DIR+_T("/main.ini"), wxEmptyString,
		wxCONFIG_USE_RELATIVE_PATH | wxCONFIG_USE_NO_ESCAPE_CHARACTERS
	);
	conf->SetRecordDefaults(); //write defaults to config file
	
	//---
	conf->SetPath(_T("/plugins") );

	wxString name[hsstab::plgNUM];
	name[hsstab::plgMF] = conf->Read(_T("mean_flow"), cmpnts::MF_HSFLOW3D_NAME);
	name[hsstab::plgLS] = conf->Read(_T("loc_search"), cmpnts::PF_LOCSRCH_NAME);
	name[hsstab::plgGS] = conf->Read(_T("glob_search"), cmpnts::PF_GLOBSRCH_NAME);
	name[hsstab::plgWPTrack] = conf->Read(_T("wptrack"), cmpnts::WPTRACK_DEFAULT_NAME);

	conf->Flush();


	//
	// Load selected plugins
	//

	bool ok = G_Plugins.load_plugin(hsstab::plgMF, name[hsstab::plgMF]);
	if( ! ok ) return false;

	conf->Flush();

	// 
	ok = G_Plugins.load_plugin(hsstab::plgLS, name[hsstab::plgLS]);
	if( ! ok ) return false;

	conf->Flush();
	

	ok = G_Plugins.load_plugin(hsstab::plgGS, name[hsstab::plgGS]);
	if( !ok )	return false;
	
	ok = G_Plugins.load_plugin(hsstab::plgWPTrack, name[hsstab::plgWPTrack]);
	if( !ok ) return false;
	
	conf->Flush();

	//---
	// configure stability task
	conf->SetPath(_T("/task") );

	task::TaskNames[task::SearchInstabLoc] = _T("SearchInstabLoc");
	task::TaskNames[task::Retrace] = _T("Retrace");
	task::TaskNames[task::MaxInstabLine] = _T("MaxInstabLine");
	task::TaskNames[task::AnalyzeWChars] = _T("AnalyzeWChars");
	task::TaskNames[task::GetProfiles] = _T("GetProfiles");
	task::TaskNames[task::MPITest] = _T("MPITest");
	task::TaskNames[task::GetAmplitudeFuncs] = _T("GetAmplitudeFuncs");
	task::TaskNames[task::GetMFChars] = _T("CalcCp");
	task::TaskNames[task::CalcScalProd] = _T("CalcScalProd");

	task::SpatTimeNames[task::Spat] = _T("Spat");
	task::SpatTimeNames[task::Time] = _T("Time");

	wxString strTaskType; 

	conf->Read(_T("task_type"), &strTaskType, task::TaskNames[task::SearchInstabLoc]);

	conf->Flush();
	
	g_taskParams.id = get_task_id(strTaskType);

	wxString strSpatTime;

	conf->Read(_T("SpatOrTimeApproach"), &strSpatTime, task::SpatTimeNames[task::Spat]);

	g_taskParams.spattime = get_spattime_id(strSpatTime);

	conf->Read(_T("pave_points_fname"), &g_taskParams.pave_grd_fname, _T("pave_points.dat"));

	const int zero=0;

	// TODO: split params into task-specific groups
	//if (g_taskParams.id==task::SearchInstabLoc || g_taskParams.id == task::MPITest){
	if (true){

		conf->Read(_T("a_ndim_min"), &g_taskParams.a_ndim_min, 1.0e-06);
		conf->Read(_T("a_ndim_max"), &g_taskParams.a_ndim_max ,1.0);
		g_taskParams.N_a = conf->Read(_T("N_a"), 10);


		conf->Read(_T("b_ndim_min"), &g_taskParams.b_ndim_min, 0.0e+00);
		conf->Read(_T("b_ndim_max"), &g_taskParams.b_ndim_max, 1.0e+00);
		g_taskParams.N_b = conf->Read(_T("N_b"), 10);

		conf->Read(_T("w_ndim_min"), &g_taskParams.w_ndim_min, 0.0e+00);
		conf->Read(_T("w_ndim_max"), &g_taskParams.w_ndim_max, 1.0e+00);
		g_taskParams.N_w = conf->Read(_T("N_w"), 10);

		conf->Read(_T("pave_point_id"), &g_taskParams.pave_point_id, zero);
	}

	//if (g_taskParams.id == task::GetAmplitudeFuncs){
	if (true){
			conf->Read(_T("pave_point_id"), &g_taskParams.pave_point_id, zero);
	}

	//if (g_taskParams.id==task::Retrace)
	if (true){
		conf->Read(_T("retrace_mode"), &g_taskParams.retrace_mode, zero);

		conf->Read(_T("pave_point_id"), &g_taskParams.pave_point_id, zero);

		conf->Flush();

	}

	delete conf;

	return true;
}
//-----------------------------------------------------------------------------


//--- Solver parameters setting
/*
static void SolverSettings()
{
	BLOCKD();  //default params

	//setInfo("NOERR/NOWARN/MEMORY/ALGEQ/LINEQ/PERMUT");
	setInfo("ERROR/WARNIN/MEMORY/ALGEQ/LINEQ/PERMUT");  // диагностика всех ошибок
	setInfo("NOSTR1/NOSTR2/NOSTR3");  // отключение сообщений о структуре матрицы Якоби
	
	// Дискретизация уравнений в частных производных
	// ПЕРЕНЕСЕНО В TDiscretPlugin

	//Итерационное решение нелинейной системы алгебраических уравнений методом Ньютона:
	setNwd(1, 1); //номер режима:0 - Ньютон-Рфасон с пересчётом матрицы Якоби
	//                                     1 - модиф. Ньютон-Рафсон с автоматическим пересч. Якоби
	//                                     2 - сначала стандартный Ньютон-Рафсон, затем модиф.
	setNwd(3,  g_genOpts.iterOut); // макс число итераций Ньютона
	setNwd(5,  g_genOpts.itsFalseJac); //допустимое число итераций без пересчёта матрицы Якоби
	setNwd(6,  g_genOpts.epsOut);  // предельная (желаемая) погрешность итерации
	setNwd(9,  g_genOpts.startNwtM);   //нач значение параметра регуляризации метода Ньютона
	setNwd(10, g_genOpts.startNwtM);   //мин значение параметра регуляризации метода Ньютона
	setNwd(12, 0.9);     //критерий обновления матрицы Якоби (знам до велич)

	//Итерационное решение системы линейных алгебраических уравнений методом GMRES:
	setItD("GMR01D", 1, 0);      //номер режима: 0 - конец по погрешности, 1 - конец по числу итераций
	setItD("GMR01D", 2, g_genOpts.iterIn);  //максимальное число итераций
	setItD("GMR01D", 3, g_genOpts.epsIn);   //предельная погрешность
	setItD("GMR01D", 4, 4);      //размерность подпространства Крылова
}
*/
//-----------------------------------------------------------------------------