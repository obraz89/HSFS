#include "log.h"

bool t_Log::_binded = false;
std::wstring t_Log::_log_path(_T(""));
std::wofstream t_Log::_ofstr;