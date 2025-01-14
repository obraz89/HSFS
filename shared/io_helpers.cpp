#include "stdafx.h"

#include "io_helpers.h"

#include <stdlib.h>
#include <string>

#include "wx/tokenzr.h"

//#include "atlbase.h"
//#include "atlstr.h"
//#include "comutil.h"

#include "wx/string.h"

// fully templated)
const int std_manip::FIELD_WIDTH_DEFAULT=12;
const int std_manip::PRECISION_DEFAULT=6;

std::istream& io_hlp::eat_white(std::istream& istr){
	char ch;
	while (istr.get(ch)){
		if (isspace(ch)==0){
			istr.putback(ch);
			break;
		};
	};
	return istr;
}

std::wistream& io_hlp::eat_white(std::wistream& istr){
	wxChar ch;
	while (istr.get(ch)){
		if (isspace(ch)==0){
			istr.putback(ch);
			break;
		};
	};
	return istr;
}

std::string wx_to_stdstr(const wxString& wx_str){

	wxLogMessage(_T("wx_to_std_str removed! - use unicode io!"));
	return std::string();
	// unicode
	// win only
	/*
		const wxChar* orig = wx_str.c_str();
		size_t origsize = wcslen(orig) + 1;
		size_t newsize = 2*origsize;
		size_t convertedChars = 0;
		char* nstring = new char[newsize];
		wcstombs_s(&convertedChars, nstring, origsize, orig, _TRUNCATE);
		std::string ret(nstring);
		delete nstring;
		return ret;
	*/
};

void io_hlp::read_parse_str_array(
				const wxString& raw_str, wxChar Delim, std::vector<std::string>& dest){

	wxArrayString wxNames = wxStringTokenize(raw_str, Delim);

	dest.clear();

	for (int i=0; i<wxNames.Count(); i++) {

		wxString& rStr = wxNames[i];

		// trim from both left and right
		rStr.Trim(true);rStr.Trim(false);
		char strName[64];

		sprintf(strName, rStr.ToAscii());
		dest.push_back(std::string(strName));

	}


}




