#include "stdafx.h"

#include "cmd_parse.h"
#include <wx/cmdline.h>

#include <wx/file.h>
#include <wx/utils.h>
#include <wx/datetime.h>

#include "log.h"

#include "mpi.h"

static const wxString strTITLE = _("\n\tHSFlowStab: High Speed Flow Stability solver. (c) 2010-2013 TsAGI, NIO-8\n");

void processCmdLine(int argc, wxChar* wxArgv[])
{

	wxLog* logCout = new TLogStream(&std::_tcout);
	wxLog::SetActiveTarget(logCout);

	wxLogMessage(strTITLE);

	static const wxCmdLineEntryDesc cmdLineDesc[] =
	{
		{ wxCMD_LINE_SWITCH, _("h"), _("help"), _("show this help message"),	wxCMD_LINE_VAL_NONE, wxCMD_LINE_OPTION_HELP },
		{ wxCMD_LINE_SWITCH, _("v"), _("version"), _("print version and exit") },
		{ wxCMD_LINE_OPTION, _("l"), _("log"), _("log file") },
		{ wxCMD_LINE_PARAM, NULL, NULL, _("case_directory") },
		{ wxCMD_LINE_NONE }
	};

	wxCmdLineParser cmdline(cmdLineDesc, argc, wxArgv);
	cmdline.SetSwitchChars(_T("-"));
	wxMessageOutput::Set(new wxMessageOutputStderr);  // for CmdLineParser msgs output

	switch( cmdline.Parse() )
	{
	case 0:		break;
	case -1:	exit(EXIT_SUCCESS);  // help
	default:	exit(EXIT_FAILURE);  // syntax error
	}

	if( cmdline.Found(_T("v")) )
	{
		wxLogMessage(_("\tBuild date: %s"), __TDATE__);
		exit(EXIT_SUCCESS);
	}


	// Case dir
	wxFileName projDir(cmdline.GetParam(), wxEmptyString);   projDir.MakeAbsolute();
	if( ! projDir.DirExists() )
	{
		wxLogError(_("Case dir '%s' doesn't exist."), projDir.GetFullPath().c_str());
		exit(EXIT_FAILURE);
	}
	projDir.SetCwd();
	wxFileName::SetCwd(projDir.GetFullPath());


	// Logging to file
	wxString logFNBase;
	if( cmdline.Found(_T("l"), &logFNBase) )
	{
		int mpi_rank;
		wxString logFN;
		wxString LogDirStr = _T("Log");
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		if (!wxFileName::DirExists(LogDirStr))
			if (!wxFileName::Mkdir(LogDirStr, 0755))
			{
				wxLogError(_("Can't create Log dir "));
				exit(EXIT_FAILURE);
			}
		logFN.Printf(_T("%s/%s%d"), LogDirStr,logFNBase, mpi_rank);
		wxLog* logFile = new TLogFile(logFN);
		wxLog::SetActiveTarget(logFile);
		wxLogMessage(strTITLE);

		new wxLogChain(logCout);
	}
}