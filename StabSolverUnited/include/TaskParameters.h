#include "MF_Field.h"
const double MF_Field::Theta = //0.122173047; 	// 7deg
					           //0.13962634;	// 8deg
					             0.087266462;	// 5deg 
const double MF_Field::Re = 1.424E+07;
const double MF_Field::L_ref =0.381;
const double MF_Field::Gamma = 1.4;
const double MF_Field::Mach = 3.5;
const double MF_Field::Alpha = 1.0/57.29577951;
const double MF_Field::T_inf = 90.318;
const double MF_Field::T_wall = 300.0;
//612.0/(1 + 0.5*(gamma -1.)*Mach*Mach);             	// if from Stetson
const double MF_Field::T_mju = 110.4/T_inf;				//Sutherland law
const double MF_Field::Mju_pow = 0.7;					// power law
