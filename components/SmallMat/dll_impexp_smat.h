///////////////////////////////////////////////////////////////////////////////
// Name:        dll_import_export.h
// Purpose:     Preprocessor defines for import/export data from dll
// Author:      A.Obraz
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <wx/wxprec.h>

#if defined(SMALLMAT_EXPORT)
#	define IMPEXP_SMALLMAT WXEXPORT
#else
#	define IMPEXP_SMALLMAT WXIMPORT
#endif
//-----------------------------------------------------------------------------