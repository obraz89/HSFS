#ifndef __COMMON_STRUCT
#define __COMMON_STRUCT

#include "dll_impexp_shared.h"

#include <map>

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

#endif		// __COMMON_STRUCT