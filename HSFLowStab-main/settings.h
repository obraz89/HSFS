///////////////////////////////////////////////////////////////////////////////
// Name:        settings.h
// Purpose:     Local settings for main module of HSFlow app
// Author:      Andrey V. Novikov
// Modified by:
///////////////////////////////////////////////////////////////////////////////

#pragma once

struct TgenericSettings
{
	//--- Init ---
	wxString strInitFieldFN; //файл начального поля, если пустой - то однородное
	wxString strPrevFieldFN; //файл поля на предыдущем шаге по времени, если пустой - то апроксимация
	bool initTimeFromFile; //для поля из файла: 0 - обнуление времени, 1 - взять из файла
	wxString strGridFN;  //имя файла сетки

};

namespace task{

	enum TTaskType { SearchInstabLoc=0, Retrace, MaxInstabLine, AnalyzeWChars, GetProfiles, TaskNum};

	extern wxString TaskNames[TaskNum];

	enum TSpatTime {Spat=0, Time, SpatTimeNum};

	extern wxString SpatTimeNames[SpatTimeNum];

	struct TTaskParams{

		int id;
		wxString pave_grd_fname;
		int spattime;

		double a_ndim_min, a_ndim_max,
			   b_ndim_min, b_ndim_max,
			   w_ndim_min, w_ndim_max;

		int pave_point_id;
		int N_a, N_b, N_w;

		int retrace_mode;

	};

}

extern task::TTaskParams g_taskParams;

//-----------------------------------------------------------------------------

extern TgenericSettings g_genOpts;

extern const wxChar strFN_HALF_STEP[];  //имя файла для 1/2 шага по времени
//-----------------------------------------------------------------------------
