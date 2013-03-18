///////////////////////////////////////////////////////////////////////////////
// Name:        dll_import_export.h
// Purpose:     Preprocessor defines for import/export data from dll
// Author:      A.Obraz
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <wx/wxprec.h>

#if defined(MFHSFLOW_EXPORT)
#	define IMPEXP_MFHSFLOW WXEXPORT
#else
#	define IMPEXP_MFHSFLOW WXIMPORT
#endif
//-----------------------------------------------------------------------------
