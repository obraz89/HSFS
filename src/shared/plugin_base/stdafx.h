// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include <wx/wxprec.h>

#include <wx/intl.h>
#include <wx/log.h>

#include <wx/filename.h>

#include <iostream>
#include <string>
#include <map>

#include <math.h>

#ifdef _UNICODE
#define _tostream wostream
#define _tcout wcout
#define _tcerr wcerr
#else
#define _tostream ostream
#define _tcout cout
#define _tcerr cerr
#endif

// TODO: reference additional headers your program requires here
