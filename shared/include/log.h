#ifndef __LOGGER
#define __LOGGER

#include <iostream>
#include <fstream>
#include "wx/string.h"

#include "string_conv.h"
class t_Log{
	static std::string _log_path;
	static bool _binded;
	static std::ofstream _ofstr;
public:
	static void bind(wxString log_path){
		if (!_binded){
			_log_path = wx_to_stdstr(log_path);
			_ofstr.open(&_log_path[0], std::ios::app);
			_binded = true;
		};
	};
	template<typename T> static void log(T val){
		_ofstr<<val;
		_ofstr.flush();
	};
	template<typename T> static void log_dup(T msg){
		_ofstr<<msg;
		_ofstr.flush();
		std::cout<<msg;
	};

	template<typename T> friend t_Log& operator<<(t_Log& log, T val){
		log_dup(val);
		return log;
	};
};

#endif // __LOGGER