///////////////////////////////////////////////////////////////////////////////
// Name:        dll_import_export.h
// Purpose:     Preprocessor defines for import/export data from EXE
// Author:      Andrey V. Novikov
// Modified by:
///////////////////////////////////////////////////////////////////////////////

#pragma once

// SHARED_EXPORT should be defined only in HSFlow.exe

#if defined(SHARED_EXPORT)
#	define DLLIMPEXP WXEXPORT
#else
#	define DLLIMPEXP WXIMPORT
#endif
//-----------------------------------------------------------------------------
