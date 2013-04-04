#include "stdafx.h"
#include "LocSearchBase.h"

using namespace stab;


t_LSCond::t_LSCond(int cnd):_cond(cnd){};

t_LSCond::t_LSCond(int cnd, t_WaveChars a_wchars):_cond(cnd),wchars(a_wchars){};

int t_LSCond::get_mode() const{return _cond;};

t_LSBase::~t_LSBase(){};