#include "TaskParameters.h"
#include "MF_Field.h"
#include "StabField.h"
#include "ODES_Stab.h"
#include "StabSolver.h"


#include "WavePackLine.h"
#include "StabField.h"

// debug
#include "EigenGS.h"

// for console io operations
#include <iostream>
#include <fstream>
#include <sstream>
#include <conio.h>
#include <stdio.h>

// support console
#ifndef _USE_OLD_OSTREAMS
using namespace std;
#endif

// enable console in GUI application
#include "console.h"

// GUI application
#include "windows.h"
//#include "AppManager.h"

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

	int nx=81, ny=161, nz=51;
	std::string NSFieldName;					// std::cin>> NSFieldName;
	NSFieldName = "input/new/04.61500.dat";		// now it is for al=2
	FILE* file = fopen("output/transitions.dat", "a+");
// read-process raw field
	MF_Field field(NSFieldName,nx,ny,nz);
	t_StabField stab_field(nx, nz);
	field.trans_to_cyl();
// set up ODES & StabSolver
	t_StabSolver stab_solver(field, stab_field);
// iterate over start positions for wave pack lines
	// DEBUG
	stab_solver.set3DContext(70, 50, 150);
	t_WaveChars w_init;
	w_init.w = t_Complex(3.85e-2, 2.55e-3);
	w_init.a = 6.26e-2;
	w_init.b = 0.1573;
	t_Complex base_resid = stab_solver.solve(w_init);
	std::cout<<"\nBase Resid:"<<base_resid<<std::endl;
	stab_solver.adjustLocal(w_init, 2);
	/*t_EigenGS gs_solver(field, 5);
	int gs_nnodes = 81;
	gs_solver.setContext(70, 50, w_init.a.real(), w_init.b.real(), gs_nnodes);
	gs_solver.getSpectrum();
	std::string f_name;
	std::ostringstream _str;
	_str<<"spectrum_N="<<gs_nnodes<<".dat";
	f_name = _str.str();
	gs_solver.writeSpectrum(&f_name[0]);
	*/
/*	
	for (int k_start = 50; k_start>2; k_start--){
		file = fopen("output/transitions.dat", "a+");
		std::cout<<"------------------------------------k_start="<<k_start<<"\n";
		int i_start = 5;
		double x_tr, t_tr;
//		WavePackLine wp_line(field, stab_field,70, 50, k_start);
//		wp_line.find_transition_location(x_tr, t_tr);
		//wp_line.print_line_to_file();
		//to_f_trans<<x_tr<<"\t"<<t_tr<<"\n";
//		fprintf(file,"%f\t%f\n", x_tr, t_tr);
//		fclose(file);
	}
*/	
	//to_f_trans.close();

	getchar();
	return 0;
}