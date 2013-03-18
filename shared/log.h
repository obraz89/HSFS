#ifndef __LOGGER
#define __LOGGER

#include <iostream>
#include <fstream>

#include "wx/file.h"
#include "wx/log.h"
#include "wx/string.h"

#include "io_helpers.h"
namespace log_my{
	inline void wxLogMessageStd(const std::wstring& str){
		const wxChar* pStr = &(str[0]);
		wxLogMessage(pStr);
	};
}


//----------------------------from Nova-----------------------------------------


#ifdef _UNICODE
#define _tostream wostream
#define _tcout wcout
#define _tcerr wcerr
#else
#define _tostream ostream
#define _tcout cout
#define _tcerr cerr
#endif

//-----------------< For wxWidgets logging support >---------------------------

//Log to an "ostream", cerr by default
class IMPEXP_SHARED TLogStream : public wxLog
{
private:
	//using ptr here to avoid including <iostream.h> from this file
	std::_tostream* m_ostr;

public:
	// redirect log output to an ostream
	TLogStream(std::_tostream* ostr = NULL);

protected:
	void DoLogString(const wxChar* szString, time_t WXUNUSED(t));
};
//-----------------------------------------------------------------------------

//Log to file
class IMPEXP_SHARED TLogFile : public wxLog
{
	wxFile m_file;

public:
	TLogFile(const wxChar* fileName);

protected:
	void DoLogString(const wxChar* szString, time_t WXUNUSED(t));
};
//-----------------------------------------------------------------------------


#endif // __LOGGER