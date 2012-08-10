#ifndef __IO_HELPERS
#define __IO_HELPERS
#include <iostream>
#include "math_operands.h"
namespace std_manip{
	extern const int FIELD_WIDTH_DEFAULT;
	extern const int PRECISION_DEFAULT;
	template<typename T> inline std::ostream& _format_fixed(std::ostream& os, T val){
		os.width(FIELD_WIDTH_DEFAULT);
		os.precision(PRECISION_DEFAULT);
		int old_flags = os.flags(std::ios::left|std::ios::fixed);
		os<<val;
		os.flags(old_flags);
		return os;
	};
	template<typename T>class t_Omanip{
		T _val;
		std::ostream& (*_formatter)(std::ostream &os, T val);
	public:
		t_Omanip(T val, std::ostream& (*formatter)(std::ostream&, T))
			:_val(val), _formatter(formatter){};
		friend std::ostream& operator<<(std::ostream& os, t_Omanip<T> m){return m._formatter(os, m._val);};
	};
	// TODO:think how to introduce options
	// and make this all useful)))
	template<typename T>inline t_Omanip<T> std_format_fixed(T val){return t_Omanip<T>(val,_format_fixed<T>);}; 
};
#endif // __IO_HELPERS