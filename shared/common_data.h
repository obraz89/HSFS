#ifndef __COMMON_DATA
#define __COMMON_DATA

#include "dll_impexp_shared.h"

#include "wx/string.h"
#include <map>

namespace hsstab{

	IMPEXP_SHARED extern wxString CASE_SETTINGS_DIR;
	IMPEXP_SHARED extern wxString OUTPUT_DIR;
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

class IMPEXP_SHARED t_Enum{
	friend class t_CompParamInt;
protected:
	std::map<int, wxString> _mapVals;
	virtual void _init_map_vals()=0;
	int _curVal;
	int* _get_val_addr(){return &_curVal;};
public:
	virtual void set_value(int val){
		// enum)))
		if(_mapVals.find(val)==_mapVals.end()) return;
		_curVal=val;
	};
	virtual int get_value(){return _curVal;};
	virtual bool operator==(int val) const{return _curVal==val;};
	virtual void operator=(const int& val){set_value(val);};
};

//
// Problem solving state
//
struct TState
{
	int mpiRank, mpiNProcs;  // MPI rank, number of procs

	//int nTmStep;     // current time step number
	//int nwtIter;     // current Newton's iteration number

	//double residual; // L_inf norm of residual at previous non-linear iteration
};

IMPEXP_SHARED extern TState G_State;


#endif	// __COMMON_DATA