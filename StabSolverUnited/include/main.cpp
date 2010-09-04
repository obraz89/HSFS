#include "TaskParameters.h"
#include "MF_Field.h"
#include <iostream>
int main(){
	int nx=81, ny=161, nz=51;
	std::string NSFieldName;
	// std::cin>> NSFieldName;
	NSFieldName = "input/new/07.85000.dat";		// now it is for al=1
	MF_Field field(NSFieldName,nx,ny,nz);
	field.trans_to_cyl();
	// all battle will be here
	getchar();
	return 0;
}