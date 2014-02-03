#include "stdafx.h"
#include "LocSearchBase.h"

using namespace stab;

//---------------------------------------------------------------------t_LSCond
t_LSCond::t_LSCond(int cnd):_cond(cnd){};

t_LSCond::t_LSCond(int cnd, t_WaveChars a_wchars):_cond(cnd),wchars(a_wchars){};

int t_LSCond::get_mode() const{return _cond;};

//---------------------------------------------------------------------t_LSMode
t_LSMode::t_LSMode(int mode):_mode(mode){};

int t_LSMode::get_mode() const{return _mode;};

void t_LSMode::set_defaults(){_mode = ASYM_HOMOGEN;};

bool t_LSMode::is_flag_on(int flag) const{return _mode&flag;}
//---------------------------------------------------------------------t_LSBase

t_LSBase::~t_LSBase(){};

void t_LSBase::setLSMode(const t_LSMode& mode){_ls_mode = mode;}