#include "boost/lexical_cast.hpp"

#include <string>
//#include <iostream>
#include <sstream>

namespace parse{
	template<typename T>
	inline T get_val(std::string a_str){
		std::stringstream s_buf, s_out;
		s_buf<<a_str;
		char c;
	// eat leading spaces
		while (s_buf.get(c)){
			if (isspace(c)==0){
				s_buf.putback(c);
				break;
			}
		};
	// read until first comment character
	//	(until first generalized space)
		while (s_buf.get(c)){
			if (isspace(c)){
				break;
			}else{
				s_out<<c;
			};
		}
		try{
			return boost::lexical_cast<T, std::string>(s_out.str());
		}
		catch(boost::bad_lexical_cast){
			std::cerr<<"Error parsing input string:\n\t"<<s_out.str()<<std::endl;
			// abort ...
		}
	};
};