///////////////////////////////////////////////////////////////////////////////
// Name:        dll_import_export.h
// Purpose:     Preprocessor defines for import/export data from dll
// Author:      A.Obraz
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <wx/wxprec.h>

#if defined(PHYSCOMMON_EXPORT)
#	define IMPEXP_PHYSCOMMON WXEXPORT
#else
#	define IMPEXP_PHYSCOMMON WXIMPORT
#endif
//-----------------------------------------------------------------------------
