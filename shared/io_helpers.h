#ifndef __IO_HELPERS
#define __IO_HELPERS

#include "dll_impexp_shared.h"

#include <iostream>

#include <vector>

class wxString;
IMPEXP_SHARED std::string wx_to_stdstr(const wxString& wx_str);

namespace io_hlp{
	IMPEXP_SHARED std::istream& eat_white(std::istream& istr);
	template<typename T> 
		inline std::istream& write_to_val(std::istream& istr, T& val){
			io_hlp::eat_white(istr);
			istr>>val;
			return istr;
		};

	IMPEXP_SHARED std::wistream& eat_white(std::wistream& istr);
	template<typename T> 
	inline std::wistream& write_to_val(std::wistream& istr, T& val){
		io_hlp::eat_white(istr);
		istr>>val;
		return istr;
	};

	IMPEXP_SHARED void read_parse_str_array(
		const wxString& raw_str, wxChar Delim, std::vector<std::string>& dest);


}

inline std::wostream& operator<<(std::wostream& ostr, const std::string& str){
	return ostr<<&(str[0]);
};


namespace std_manip{
	IMPEXP_SHARED extern const int FIELD_WIDTH_DEFAULT;
	IMPEXP_SHARED extern const int PRECISION_DEFAULT;
	template<typename T> inline std::wostream& _format_fixed(std::wostream& os, T val){
		os.width(FIELD_WIDTH_DEFAULT);
		os.precision(PRECISION_DEFAULT);
		std::ios_base::fmtflags old_flags = os.flags(std::ios::left|std::ios::fixed);
		os<<val;
		os.flags(old_flags);
		return os;
	};

	template<typename T> inline std::wostream& _format_sci(std::wostream& os, T val){
		os.width(FIELD_WIDTH_DEFAULT);
		os.precision(PRECISION_DEFAULT);
		std::ios_base::fmtflags old_flags = os.flags(std::ios::left|std::ios::scientific);
		os<<val;
		os.flags(old_flags);
		return os;
	};
	template<typename T>class t_Omanip{
		T _val;
		std::wostream& (*_formatter)(std::wostream &os, T val);
	public:
		t_Omanip(T val, std::wostream& (*formatter)(std::wostream&, T))
			:_val(val), _formatter(formatter){};
		friend std::wostream& operator<<(std::wostream& os, t_Omanip<T> m){return m._formatter(os, m._val);};
	};
	// TODO:think how to introduce options
	// and make this all useful)))
	template<typename T>inline t_Omanip<T> 
		std_format_fixed(T val){return t_Omanip<T>(val,_format_fixed<T>);}; 

	template<typename T>inline t_Omanip<T> 
		std_format_sci(T val){return t_Omanip<T>(val,_format_sci<T>);}; 
};
#endif // __IO_HELPERS