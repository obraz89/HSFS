///////////////////////////////////////////////////////////////////////////////
// Name:        dll_import_export.h
// Purpose:     Preprocessor defines for import/export data from dll
// Author:      A.Obraz
///////////////////////////////////////////////////////////////////////////////

#include <wx/wxprec.h>

#pragma once

#if defined(SHARED_EXPORT)
#	define IMPEXP_SHARED WXEXPORT
#else
#	define IMPEXP_SHARED WXIMPORT
#endif
//-----------------------------------------------------------------------------
