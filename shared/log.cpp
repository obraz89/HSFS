#include "stdafx.h"
#include "log.h"

TLogStream::TLogStream(std::_tostream* ostr/* = NULL*/)
{
	if(ostr == NULL)  
		m_ostr = &std::_tcerr;
	else
		m_ostr = ostr;
};


void TLogStream::DoLogString(const wxChar* szString, time_t WXUNUSED(t))
{
#ifdef __WINDOWS__
	char* strBuf = new char[_tcslen(szString)+1];
	::CharToOem(szString, strBuf);
	std::cout << strBuf << std::endl;	
	delete[] strBuf;
#else
	(*m_ostr) << wxConvertWX2MB(szString) << std::endl;	
#endif
};

TLogFile::TLogFile(const wxChar* fileName)
{
	m_file.Open(fileName, wxFile::write_append);
}

void TLogFile::DoLogString(const wxChar* szString, time_t WXUNUSED(t))
{
	m_file.Write(wxString(szString)+_T("\n"), *wxConvCurrent);
}

