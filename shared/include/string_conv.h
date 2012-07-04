#include <iostream>
#include <stdlib.h>
#include <string>

#include "atlbase.h"
#include "atlstr.h"
#include "comutil.h"

#include "wx/string.h"

// Convert to a char*
inline std::string wx_to_stdstr(const wxString& wx_str){
	const wchar_t* orig = wx_str.c_str();
	size_t origsize = wcslen(orig) + 1;
	size_t newsize = 2*origsize;
	size_t convertedChars = 0;
	char* nstring = new char[newsize];
	wcstombs_s(&convertedChars, nstring, origsize, orig, _TRUNCATE);
	std::string ret(nstring);
	delete nstring;
	return ret;
};


