#include "stdafx.h"
#include "LocSearchBase.h"

using namespace stab;

//---------------------------------------------------------------------t_LSCond
void t_LSCond::set(int cnd){_cond=cnd;}

t_LSCond::t_LSCond(){}
t_LSCond::t_LSCond(int cnd){set(cnd);};

void t_LSCond::set(int cnd, const t_WaveChars& a_wchars){
	_cond=cnd;
	wchars=a_wchars;
}
t_LSCond::t_LSCond(int cnd, const t_WaveChars& a_wchars){set(cnd, a_wchars);};

int t_LSCond::get_mode() const{return _cond;};

//---------------------------------------------------------------------t_LSMode
t_LSMode::t_LSMode(int mode):_mode(mode){};

int t_LSMode::get_mode() const{return _mode;};

// by default, direct problem with forced homogen asymptotics
void t_LSMode::set_defaults(){_mode = DIRECT|ASYM_HOMOGEN;};

bool t_LSMode::is_flag_on(int flag) const{ return _mode&flag;}
//---------------------------------------------------------------------t_LSBase

t_LSBase::~t_LSBase(){};

void t_LSBase::setLSMode(const t_LSMode& mode){_ls_mode = mode;}

void stab::dumpEigenFuncs(const std::string& fname, const int nnodes, const std::vector<double>& y_vec, const std::vector<t_VecCmplx>& amp_funcs) {

	const int STAB_MATRIX_DIM = 8;

	std::wofstream fstr(&fname[0]);
	fstr << _T("Y\tu_re\tu_im\tu'_re\tu'_im\tv_re\tv_im\tp_re\tp_im\tt_re\tt_im\tt'_re\tt'_im\tw_re\tw_im\tw'_re\tw'_im\n");

	// write out in reverse order (from wall to outer)
	for (int j = nnodes - 1; j >= 0; j--) {

		fstr << y_vec[j] << _T("\t");

		for (int k = 0; k<STAB_MATRIX_DIM; k++) {
			fstr << std_manip::std_format_sci<double>(amp_funcs[j][k].real())
				<< _T("\t")
				<< std_manip::std_format_sci<double>(amp_funcs[j][k].imag())
				<< _T("\t");
		}
		fstr << _T("\n");
	}

};