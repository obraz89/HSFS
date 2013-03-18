#ifndef MFPROC_EXP
#define MFPROC_IMPEXP __declspec(dllimport)
#else
#define MFPROC_IMPEXP __declspec(dllexport)
#endif