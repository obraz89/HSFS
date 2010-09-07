#include "TaskParameters.h"
#include "MF_Field.h"
#include "StreamLine.h"
#include <iostream>

int main(){
	int nx=81, ny=161, nz=51;
	std::string NSFieldName;					// std::cin>> NSFieldName;
	NSFieldName = "input/new/07.85000.dat";		// now it is for al=1
	MF_Field field(NSFieldName,nx,ny,nz);
	field.trans_to_cyl();
	// initialize Streamline
	int i_start = 5, k_start = 50;
	StreamLine str_line(field, i_start, field.get_bound_index(i_start, k_start)+5, k_start);
	for (int i=0; i<200; i++) str_line.add_node();
	Index pos_ind = str_line.find_transition_location();
	str_line.print_line("v");
	getchar();
	return 0;
}