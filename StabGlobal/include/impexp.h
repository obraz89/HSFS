#ifndef STABGLOBAL_EXP
#define STABGLOBAL_IMPEXP __declspec(dllimport)
#else
#define STABGLOBAL_IMPEXP __declspec(dllexport)
#endif