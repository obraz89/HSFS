///////////////////////////////////////////////////////////////////////////////
// Name:        dll_import_export.h
// Purpose:     Preprocessor defines for import/export data from dll
// Author:      A.Obraz
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <wx/wxprec.h>

#if defined(PROFILE_EXPORT)
#	define IMPEXP_PROFILE WXEXPORT
#else
#	define IMPEXP_PROFILE WXIMPORT
#endif
//-----------------------------------------------------------------------------
