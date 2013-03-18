#include "stdafx.h"
#include "log.h"

// support console
#ifndef _USE_OLD_OSTREAMS
using namespace std;
#endif

#include "cmd_parse.h"

#include "tests.h"

bool load_Settings_n_Plugins();

int main(int argc, char* argv[]){
	int err = 0;

#ifdef wxUSE_UNICODE
	wxChar **wxArgv = new wxChar*[argc + 1];
	{
		int n = 0;
		for(n=0; n<argc; n++)
		{
			wxMB2WXbuf warg = wxConvertMB2WX(argv[n]);
			wxArgv[n] = wxStrdup(warg);
		}
		wxArgv[n] = NULL;
	}
#else
#define wxArgv argv
#endif // wxUSE_UNICODE/!wxUSE_UNICODE

	processCmdLine(argc, wxArgv);
	wxLogMessage(
		_("\t-\n\t- Run %s\n\t-\n"),
		wxDateTime::Now().Format(_T("%Y-%m-%d %H:%M:%S")).c_str()
		);

	load_Settings_n_Plugins();

	test::selfsim_M45_second_mode();
}
/*
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
	load_Settings_n_Plugins();

	//test::king_al_2_new();

	//test::transhyb_base_08();

	//test::itam_hz();

	//test::profile_compar();

	test::selfsim_M45_second_mode();

	return 0;
}
*/