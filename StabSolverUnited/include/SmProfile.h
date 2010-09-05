extern "C"{
	const int SmProfSize = 541; // to match fortran common (TODO: automate ???)
	extern struct SmProfile{
	double	YNS[SmProfSize],UNS[SmProfSize],UNS1[SmProfSize],UNS2[SmProfSize],
							TNS[SmProfSize],TNS1[SmProfSize],TNS2[SmProfSize],
							WNS[SmProfSize],WNS1[SmProfSize],WNS2[SmProfSize],
							ENTNS[SmProfSize],FINS[SmProfSize],DFINS[SmProfSize],XST[SmProfSize];	
	} NS;
	extern struct{
		double RHO[SmProfSize];
	} RHONS; // density is in separate common for some reason
}