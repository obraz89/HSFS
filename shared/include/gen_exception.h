#ifndef __EXCEPTIONS_BASE__
#define __EXCEPTIONS_BASE__
#include "wx/string.h"
#include <wx/wxprec.h>
#define ssuGENTHROW(...)  \
	throw t_GenException( wxString::Format(__VA_ARGS__), __TFILE__, __LINE__ )

#define ssuTHROW(E, msg)  \
	throw E(msg, __TFILE__, __LINE__)
//------

class t_GenException
{
protected:
	wxString    _what;
	wxString    _file;
	int         _line;

public:
	t_GenException(const wxString& what, const wxChar* szFile,  const int line)
		: _what(what), _file(szFile), _line(line) {   }

	wxString what() const  {  return _what;  }
	wxString file() const  {  return _file;  }
	int line() const       {  return _line;  }
	wxString what_detailed() const
	{
		return _what + _(". In file: ")+_file + _(", line: ")+wxString::Format(_T("%d"), _line);
	}
	virtual ~t_GenException(){};

	friend std::wostream& operator<<(std::wostream& ostr, const t_GenException& x){
		#ifdef _DEBUG
			return ostr<<_T("Exception")<<x.what_detailed().c_str();
		#else
			return ostr<<_T("Exception")<<x.what().c_str();
		#endif
	};
};

#endif // __EXCEPTIONS_BASE__