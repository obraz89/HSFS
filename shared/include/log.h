#ifndef __LOGGER
#define __LOGGER

#include <iostream>
#include <fstream>
#include "wx/string.h"

#include "io_helpers.h"

class t_Log{
	static std::wstring _log_path;
	static bool _binded;
	static std::wofstream _ofstr;
public:
	enum t_Mode{Dup=0, Cout} mode;
	t_Log::t_Log():mode(t_Mode::Dup){};
	t_Log& set_mode(t_Mode a_mode){mode=a_mode; return *this;};
	static void bind(wxString log_path){
		if (!_binded){
			_log_path = log_path.c_str();
			_ofstr.open(&_log_path[0], std::ios::app);
			_binded = true;
		};
	};
	template<typename T> static void log(T val){
		_ofstr<<val;
		_ofstr.flush();
	};
	template<typename T> static void log_dup(T msg){
		log(msg);
		std::wcout<<msg;
	};

	//std::wstring

	template<typename T> friend t_Log& operator<<(t_Log& log, T val){
		log_dup(val);
		return log;
	};
};

static t_Log Log;

#endif // __LOGGER