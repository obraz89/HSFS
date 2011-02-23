#include "TaskParameters.h"
#include "MF_Field.h"
//#include "StreamLine.h"		// streamline concept dropped
#include "WavePackLine.h"
#include "StabField.h"



// for console io operations
#include <iostream>
#include <fstream>
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
	MF_Field field(NSFieldName,nx,ny,nz);
	StabField stab_field(nx, nz);
	field.trans_to_cyl();
	// initialize Streamline
	for (int k_start = 50; k_start>2; k_start--){
		file = fopen("output/transitions.dat", "a+");
		std::cout<<"------------------------------------k_start="<<k_start<<"\n";
		int i_start = 5;
		double x_tr, t_tr;
		/*StreamLine str_line(field, i_start,80, k_start);
		for (int i=0; i<1000; i++) str_line.add_node();
		str_line.find_transition_location(x_tr, t_tr);*/
		WavePackLine wp_line(field, stab_field,70, 50, k_start);
		wp_line.find_transition_location(x_tr, t_tr);
		wp_line.print_line_to_file();
		//to_f_trans<<x_tr<<"\t"<<t_tr<<"\n";
		fprintf(file,"%f\t%f\n", x_tr, t_tr);
		fclose(file);
		//str_line.print_line("v");
	}
	//to_f_trans.close();
	getchar();
	return 0;
}