#include "stdafx.h"
#include "gen_exception.h"

t_GenException::t_GenException(const wxString& what, const wxChar* szFile,  const int line)
: _what(what), _file(szFile), _line(line) {};

wxString t_GenException::what() const{  return _what;  }

wxString t_GenException::file() const  {  return _file;  }

int t_GenException::line() const       {  return _line;  }

wxString t_GenException::what_detailed() const
{
	return _what + _(". In file: ")+_file + _(", line: ")+wxString::Format(_T("%d"), _line);
}

t_GenException::~t_GenException(){};

std::wostream& operator<<(std::wostream& ostr, const t_GenException& x){
	#ifdef _DEBUG
		return ostr<<_T("Exception")<<x.what_detailed().c_str();
	#else
		return ostr<<_T("Exception")<<x.what().c_str();
	#endif
};