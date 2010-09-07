extern "C"{
	const int SmProfSize = 541; // to match fortran common (TODO: how to automate ???)
	void SEARCH_MAX_INSTAB_TIME();
	void TS_GLOBAL_TIME();
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
		int  MF_I, MF_J;
		double	MF_X, MF_WE, MF_ANG;
	} CF_PROFILE_DATA;
	extern struct{
		int X_DIM,Y_DIM;
	} ADDITIONAL_NS;
	extern struct{
		int REQ_TS_GLOB;
	} CONTROL;
	extern struct{
		double SIGMA_SPAT;
	} SOLVER_OUTPUT;
}