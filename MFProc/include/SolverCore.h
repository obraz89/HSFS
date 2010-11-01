#ifndef __SOLVERCORE
#define __SOLVERCORE
#include "Elems.h"
extern "C"{
	const int SmProfSize = 541; // to match fortran common (TODO: how to automate ???)
	extern void SEARCH_MAX_INSTAB_TIME();
	extern void SEARCH_INITIAL_INSTAB_TIME();
	extern void NAVSTOK(int&, double&);
	extern struct {
	double	YNS[SmProfSize],UNS[SmProfSize],UNS1[SmProfSize],UNS2[SmProfSize],
							TNS[SmProfSize],TNS1[SmProfSize],TNS2[SmProfSize],
							WNS[SmProfSize],WNS1[SmProfSize],WNS2[SmProfSize],
							ENTNS[SmProfSize],FINS[SmProfSize],DFINS[SmProfSize],XST[SmProfSize];	
	} NS;
	extern struct{
		double RHO[SmProfSize];
	} RHONS; // density is in separate common for some reason
	extern struct{
		double	AMINF,REINF,TINF,XLL,REE,REE1,UEE;
	} DNS;
	extern struct{
		CompVal VA, VB, VR;
	} VGRC;
	extern struct{
		int  MF_I, MF_J;
		double	MF_X, MF_WE, MF_ANG;
	} CF_PROFILE_DATA;
	extern struct{
		int X_DIM,Y_DIM;
	} ADDITIONAL_NS;
	extern struct{
		CompVal SIGMA_SPAT, A_SPAT, B_SPAT, W_SPAT;
	} SOLVER_OUTPUT;
	extern struct{
		double K,SIG,G,M,XI,PB;
	} OUT;	// additional mean flow parameters
	extern struct {
		int IUPT;
	} IUPT;	// viscosity : 0-Sutherland, 1 - power
	extern struct {
		double A6[10], BVB, DL1; // BVB is for viscosity power coef ( eg 0.75)
	} BASA6;
	extern struct {
		int MAB, MAB1, MAB2;	// seems only MAB is used: see STAB.for (defines asymptotics)
	} ABOCT;
	extern struct {
		double Y1, D;
		int NS, L1;
	} HADY1;					// NS used;
	extern struct {
		double Y1O, DO, C, YC;
		int NO;					// NO used;
	} MACKY1;
	extern struct{
		CompVal A,B,W,R;
	} HADY2;
}
#endif // __SOLVERCORE