#ifndef __EXCEPTIONS_BASE__
#define __EXCEPTIONS_BASE__

#include "wx/string.h"
#include "dll_impexp_shared.h"

#define ssuGENTHROW(...)  \
	throw t_GenException( wxString::Format(__VA_ARGS__), __TFILE__, __LINE__ )

#define ssuTHROW(E, msg)  \
	throw E(msg, __TFILE__, __LINE__)
//------

class IMPEXP_SHARED t_GenException
{
protected:
	wxString    _what;
	wxString    _file;
	int         _line;

public:
	t_GenException(const wxString& what, const wxChar* szFile,  const int line);

	wxString what() const;
	wxString file() const;
	int line() const;
	wxString what_detailed() const;
	
	~t_GenException();

	IMPEXP_SHARED friend std::wostream& operator<<(std::wostream& ostr, const t_GenException& x);
};

#endif // __EXCEPTIONS_BASE__