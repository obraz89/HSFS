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

	enum TTaskType { 

		// Global search of instabilities
		SearchInstabLoc=0, 

		// Test: do global search and find max with local search (try Newton method)
		// ai = ai(b), w=fixed
		SearchMaxInstabLoc,

		// with initial wave from global search, try steepest descent to find max of ai
		// ai = ai(w,b)
		SearchMaxInstabLocGrad,

		// Search wchars (a_end, b_end, w_end) from (a_start, b_start, w_start)
		// keep br/ar fixed
		SearchWCharsLocWVecFix, 

		// Search wchars (a_end, b_end, w_end) from (a_start, b_start, w_start)
		// move directly from (b_start, w_start) to (b_end, w_end)
		SearchWCharsLocWBShift,

		// get amplification rates over wave paths with various closures
		Retrace, 

		// experimental - parallel retrace
		// when tested, replace Retrace 
		RetraceMPI,

		// retrace external (inviscid) streamlines only
		RetraceStreamlines,

		// "retrace" wall gridline
		GetWallGridLine,

		// posprocess data obtained in retraceMPI
		PostProcRetrace,

		// ?
		MaxInstabLine, 

		// ?
		AnalyzeWChars, 

		// Extract mean flow profiles
		GetProfiles, 

		// 
		MPITest, 

		// reconstruct amplitude functions for a given wave
		GetAmplitudeFuncs, 

		// get boundary layer spectrum from global search routine
		GetBLSpectrum, 

		// calculate Cp and other useful functions in pave points
		CalcMFChars, 

		// calculate scalar product of 2 modes, both are taken from LST analysis
		CalcScalProd_LST_test, 

		// calculate scalar product of 1 mode from LST and spectral component of 
		// disturbance froom DNS

		CalcScalProd_DNS_test,

		// test: calculate derivs of primitive variables at particular point
		CalcMeanFlowRecDerivs,

		// calculate neutral curve
		CalcNeutralCurve,

		// For current tests
		Test, 

		TaskNum};

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

		double a_dim_min, a_dim_max,
		       b_dim_min, b_dim_max,
			   w_dim_min, w_dim_max;

		int pave_point_id;
		int N_a, N_b, N_w;

	};

}

extern task::TTaskParams g_taskParams;

//-----------------------------------------------------------------------------

extern TgenericSettings g_genOpts;

extern const wxChar strFN_HALF_STEP[];  //имя файла для 1/2 шага по времени
//-----------------------------------------------------------------------------
