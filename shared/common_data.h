#ifndef __COMMON_DATA
#define __COMMON_DATA

#include "dll_impexp_shared.h"

#include "wx/string.h"

namespace hsstab{

	IMPEXP_SHARED extern wxString CASE_SETTINGS_DIR;
	IMPEXP_SHARED extern wxString LOG_FILE;


	namespace cmpnts{
// mf plugins
		IMPEXP_SHARED extern const wxString MF_HSFLOW3D_NAME;
		IMPEXP_SHARED extern const wxString MF_HSFLOW2D_NAME;

// pf plugins
		IMPEXP_SHARED extern const wxString PF_LOCSRCH_NAME;	
		IMPEXP_SHARED extern const wxString PF_GLOBSRCH_NAME;

		IMPEXP_SHARED extern const wxString MF_CONF_DOMAIN;
		IMPEXP_SHARED extern const wxString EIGEN_CONF_DOMAIN;
		IMPEXP_SHARED extern const wxString STABSOLVER_CONF_DOMAIN;
	};
};

#endif	// __COMMON_DATA