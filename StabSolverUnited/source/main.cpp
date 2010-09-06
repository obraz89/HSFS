#include "TaskParameters.h"
#include "MF_Field.h"
#include "StreamLine.h"
#include "SolverCore.h"
#include "SmProfile.h"
#include "Smooth.h"
#include <iostream>

int main(){
	int nx=81, ny=161, nz=51;
	std::string NSFieldName;					// std::cin>> NSFieldName;
	NSFieldName = "input/new/07.85000.dat";		// now it is for al=1
	MF_Field field(NSFieldName,nx,ny,nz);
	field.trans_to_cyl();
	// initialize Streamline
	int i_start = 5, k_start = 25;
	StreamLine str_line(field, i_start, field.get_bound_index(i_start, k_start)+25, k_start);
	for (int i=0; i<10; i++) str_line.add_node();
	str_line.print_line("v");
	SmProfile cur_profile(field, nx-2, nz/2);
	cur_profile.smooth();
	cur_profile.setSolverParameters();
	SEARCH_MAX_INSTAB_TIME();
	getchar();
	return 0;
}