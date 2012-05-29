#ifndef STABLOCAL_EXP
#define STABLOCAL_IMPEXP __declspec(dllimport)
#else
#define STABLOCAL_IMPEXP __declspec(dllexport)
#endif