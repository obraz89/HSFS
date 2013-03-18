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

using namespace hsstab;

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
	name[hsstab::plgWPTrack] = conf->Read(_T("wptrack"), _T("WPTrack"));

	delete conf;


	//
	// Load selected plugins
	//

	bool ok = G_Plugins.load_plugin(hsstab::plgMF, name[hsstab::plgMF]);
	if( ! ok ) return false;

	// 
	ok = G_Plugins.load_plugin(hsstab::plgLS, name[hsstab::plgLS]);
	if( ! ok ) return false;
	//---

	ok = G_Plugins.load_plugin(hsstab::plgGS, name[hsstab::plgGS]);
	if( !ok ) return false;

	ok = G_Plugins.load_plugin(hsstab::plgWPTrack, name[hsstab::plgWPTrack]);
	if( !ok ) return false;

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