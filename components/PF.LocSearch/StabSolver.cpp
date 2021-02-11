#include "stdafx.h"

#include "StabSolver.h"
#include "common_data.h"
#include "log.h"

using namespace hsstab;
using namespace hsstab::cmpnts;

using namespace pf;

t_StabSolver::t_StabSolver(const mf::t_DomainBase& a_rFldNS):
_rFldNS(a_rFldNS), _profStab(), _math_solver(), 
_stab_matrix(STAB_MATRIX_DIM),
_scal_prod_matrix_H1(STAB_MATRIX_DIM), _scal_prod_matrix_HW(STAB_MATRIX_DIM),
_params(){

	_math_solver._pStab_solver = this;

	_ls_mode.set_defaults();

};

void t_StabSolver::init(const hsstab::TPlugin& g_plug){

	const hsstab::TPluginParamsGroup& g = g_plug.get_settings_grp_const("");

	_params.init(g);

	//_init();	// init solver

}

const t_StabScales& t_StabSolver::get_stab_scales() const{
	return _profStab.scales();
};

const std::vector<double>& t_StabSolver::get_y_distrib() const {
	return _math_solver.varRange;
}

// private, to be used in 3D context

/*const t_SqMatrix& t_StabSolver::_getStabMatrix3D() const{
	return _stab_matrix;
}*/

void t_StabSolver::_formRHS3D(const double& a_y, const t_VecCmplx& a_vars, t_VecCmplx& dest){

	_setStabMatrix3D(a_y);

	// direct problem
	// rhs = H*phi
	if (_ls_mode.is_flag_on(stab::t_LSMode::DIRECT)){

		matrix::base::mat_mul<t_Complex, t_Complex>(_stab_matrix, a_vars, dest);

		return;

	}

	// conjugate problem
	// dpsi/dy = H_conj*psi, 
	// H_conj = -1*(H_dir)_transposed_conjugate, H_dir - matrix of direct problem
	// we will solve dpsi_conjugate/dy = A*psi_conjugate
	// then A = -1*H_dir_transposed
	// the result is psi_conjugate : complex conjugate vetor of explicit conjugate vector psi 
	// rhs = (-1.0*H_tr)*psi

	if (_ls_mode.is_flag_on(stab::t_LSMode::CONJUGATE)){

		for (int i=0; i<STAB_MATRIX_DIM; i++){

			dest[i] = 0.0;

			for (int j=0; j<STAB_MATRIX_DIM; j++) {

				dest[i] = dest[i] - _stab_matrix[i][j]*a_vars[j];

			}

		}

		return;

	}

};

/************************************************************************/   
// Restore wave chars in global reference frame
// Use time approach
/************************************************************************/

t_WCharsGlob t_StabSolver::popGlobalWCharsTime(const mf::t_GeomPoint a_xyz){

	t_WCharsLoc restore_wave = _waveChars;
	t_WCharsLoc adjust_wave = _waveChars;

	if (_waveChars.get_treat()==stab::t_TaskTreat::SPAT){

		calcGroupVelocity(adjust_wave);

		adjust_wave.to_time();

		adjustLocal(adjust_wave, t_StabSolver::t_MODE::W_MODE);

	}

	t_WCharsGlob glob_wave(adjust_wave, _rFldNS.calc_jac_to_loc_rf(a_xyz), _profStab.scales());

	_waveChars = restore_wave;

	return glob_wave;

};


/************************************************************************/   
// Restore wave chars in global reference frame
// Use spat approach
/************************************************************************/

t_WCharsGlob t_StabSolver::popGlobalWCharsSpat(const mf::t_GeomPoint a_xyz){

	t_WCharsLoc restore_wave = _waveChars;
	t_WCharsLoc adjust_wave = _waveChars;

	if (_waveChars.get_treat()==stab::t_TaskTreat::TIME){

		calcGroupVelocity(adjust_wave);

		adjust_wave.to_spat();

		// TODO: here is another uncertainty:
		// A_MODE or B_MODE or mixed - how to choose
		adjustLocal(adjust_wave, t_StabSolver::t_MODE::A_MODE);

	}

	t_WCharsGlob glob_wave(adjust_wave, _rFldNS.calc_jac_to_loc_rf(a_xyz), _profStab.scales());

	_waveChars = restore_wave;

	return glob_wave;
};

/************************************************************************/   
// Set up the required context for a given point
// 
/************************************************************************/

void t_StabSolver::setContext(const mf::t_GeomPoint a_xyz){

	_cur_xyz = a_xyz;

	int nnodes_stab = _params.NNodes;

	_math_solver.setContext(nnodes_stab);

	t_ProfMFLoc profNS(_rFldNS);

	mf::t_ProfDataCfg prof_cfg;
	prof_cfg.ThickCoef = _params.ThickCoef;

	// NNodes can be used in INTERPOLATE types of initialization
	prof_cfg.NNodes = _params.NNodes;

	switch (_params.NSProfInit)
	{
	case (blp::NSINIT_EXTRACT):
		profNS.initialize(a_xyz, prof_cfg, blp::NSINIT_EXTRACT);
		break;
	case (blp::NSINIT_INTERPOLATE):
		profNS.initialize(a_xyz, prof_cfg, blp::NSINIT_INTERPOLATE);
		break;
	default:
		wxString msg(_T("PF.LocSearch: ProfNS Initialization type not supported"));
		wxLogError(msg); ssuGENTHROW(msg);
	}

	t_ProfStabCfg pstb_cfg;

	pstb_cfg.NNodes = nnodes_stab;
	pstb_cfg.NondimScaleType = _params.NondimScaleType;

	_profStab.initialize(profNS, pstb_cfg);

	if (_stab_matrix.nCols()!=STAB_MATRIX_DIM){
		ssuTHROW(t_GenException, 
		_T("PF.LocSearch Error:Wrong size of stab matrix while setting context"));
	}

	for (int j=0; j<nnodes_stab; j++){
		_math_solver.varRange[j] = _profStab.get_y(nnodes_stab-1-j);
	}

	_setCurvCoefs(a_xyz);

	return;
}

/************************************************************************/   
// Set up context using AVF profiles
// 
/************************************************************************/

void t_StabSolver::setContext(const t_ProfileStab* a_prof_stab){

	_profStab = *a_prof_stab;

	int nnodes = _profStab.get_nnodes();

	_math_solver.setContext(nnodes);

	if (_stab_matrix.nCols()!=STAB_MATRIX_DIM){
		ssuTHROW(t_GenException, 
			_T("PF.LocSearch Error:Wrong size of stab matrix while setting context"));
	}

	for (int j=0; j<nnodes; j++){
		_math_solver.varRange[j] = _profStab.get_y(nnodes-1-j);
	}

	return;
	
}

/************************************************************************/
// set the candidate wave
//
/************************************************************************/
void t_StabSolver::setWave(const t_WCharsLoc& wave) { _waveChars = wave; }

/************************************************************************/   
//COMPUTATION OF GROUP VELOCITY (VA,VB)=(DW/DA,DW/DB) 
//
/************************************************************************/

void t_StabSolver::calcGroupVelocity
(t_WCharsLoc &a_wave_chars){
	// empiric constant - tolerance
	const double dd = 0.001*a_wave_chars.a.real();//DELTA_SMALL;
	// ensure that we are in eigen 
	adjustLocal(a_wave_chars, W_MODE);

	t_WCharsLoc left_chars = a_wave_chars;
	t_WCharsLoc rght_chars = a_wave_chars;
	// dalpha
	left_chars.a-=dd;
	rght_chars.a+=dd;
    adjustLocal(left_chars, W_MODE);
	adjustLocal(rght_chars, W_MODE);
	if (left_chars.w == rght_chars.w){
		throw t_UnPhysWave();
	};
	a_wave_chars.vga = 0.5*(rght_chars.w - left_chars.w)/dd;

	// reset left, right
	left_chars = a_wave_chars;
	rght_chars = a_wave_chars;

	left_chars.b-=dd;
	rght_chars.b+=dd;
    adjustLocal(left_chars, W_MODE);
	adjustLocal(rght_chars, W_MODE);
	a_wave_chars.vgb = 0.5*(rght_chars.w - left_chars.w)/dd;
}


/************************************************************************/   
//Check wave chars - if wave is physical, group velo is good
//
/************************************************************************/
bool t_StabSolver::checkWCharsByGroupV(t_WCharsLoc& wchars){

	try{
		calcGroupVelocity(wchars);
	}catch(...){
		wxLogMessage(_T("WChars Check: Failed to calculate group velo - unphysical wave"));
		return false;
	}

	double vga_r = wchars.vga.real();
	double vgb_r = wchars.vgb.real();

	double vgr = sqrt(vga_r*vga_r + vgb_r*vgb_r);

	if (vgr>10.0){

		wxLogMessage(_T("WChars Check: group velo is big (vgr>10.0) - unphysical wave"));
		return false;

	} 

	if (vga_r<0.0){

		wxLogMessage(_T("WChars Check: vga<0 - upstream wave"));
		return false;

	}

	return true;

};
/************************************************************************/   
// Search maximum instability with keeping Re(w)=const
// Use time approach
/************************************************************************/

t_WCharsLoc t_StabSolver::_getStationaryMaxInstabTime
(const t_WCharsLoc& initial_guess){
	t_WCharsLoc nxt_wave = initial_guess;
	t_WCharsLoc prv_wave = initial_guess;
	t_Complex gv;
	double da, db, k, eps; 

	// TODO: configurable, empirics
	// eps max - max variation of wave vector 
	const double eps_max = 0.01;
	const int eps_it_max = 2;
	int ndebug=0;

	do{
		prv_wave = nxt_wave;

		//adjustLocal(prv_wave, W_MODE);   
		getEigenWFixed(initial_guess.w.real(), prv_wave, A_MODE);
		calcGroupVelocity(prv_wave);

		k = sqrt(pow(prv_wave.a.real(),2)+pow(prv_wave.b.real(),2));

		eps = eps_max;
		int eps_iter = 0;
		bool eps_success = false;

		do 
		{
			db = eps*k;
			double coef = prv_wave.vgb.real()/prv_wave.vga.real();
			da = -coef*db;

			nxt_wave.a = prv_wave.a + da;
			nxt_wave.b = prv_wave.b + db;
			nxt_wave.w = prv_wave.w;

			//adjustLocal(nxt_wave, W_MODE);
			eps_success = getEigenWFixed(initial_guess.w.real(), nxt_wave, A_MODE);

			if (eps_success && nxt_wave.w.imag()>prv_wave.w.imag()) break;

			db = -db;
			da = -coef*db;
			nxt_wave.a = prv_wave.a + da;
			nxt_wave.b = prv_wave.b + db;
			nxt_wave.w = prv_wave.w;
			//adjustLocal(nxt_wave, W_MODE);	
			eps_success = getEigenWFixed(initial_guess.w.real(), nxt_wave, A_MODE);

			if (eps_success && nxt_wave.w.imag()>prv_wave.w.imag()) break;

			eps/=2.;
		} while (eps_iter++<eps_it_max);

		if (eps_success==false) break;

		double wr_error = abs(1.0-nxt_wave.w.real()/initial_guess.w.real());
		std::wcout<<_T("getMaxWave Fixed Wr: iter err~")<<wr_error<<_T("\n Cur wave:")<<nxt_wave;
		
		if (ndebug++==100){
			ndebug=0;
			std::wcout<<_T("getStatMaxInstabTime 100 iters done, cur wave:\n")<<nxt_wave<<_T("\n");
			break;
		}
	} while(prv_wave.w.imag()<nxt_wave.w.imag());

	return prv_wave;
};

/************************************************************************/   
// Search maximum instability without any restrictions 
// Use time approach
/************************************************************************/

t_WCharsLoc t_StabSolver::_getMaxInstabTime_Stat
(const t_WCharsLoc &init_guess){
	// empiric half percent per iteration
	double dar = 0.005*init_guess.a.real();
	t_WCharsLoc base_wave = _getStationaryMaxInstabTime(init_guess);
	calcGroupVelocity(base_wave);
	t_WCharsLoc next_wave = base_wave;
	next_wave.a+=dar;
	next_wave.w+=base_wave.vga*dar;
	next_wave = _getStationaryMaxInstabTime(next_wave);
	// do we move in proper direction ?
	if (next_wave.w.imag()<base_wave.w.imag()){
		dar = -dar;
	};
	
	do{
		base_wave = next_wave;
		calcGroupVelocity(base_wave);
		next_wave.a+=dar;
		next_wave.w+=base_wave.vga*dar;
		next_wave = _getStationaryMaxInstabTime(next_wave);
	} while (base_wave.w.imag()>next_wave.w.imag());
	return base_wave;
};


/************************************************************************/   
// Search for exact eigen solution given initial approach a_wave_chars
// input: mode defines the value to be adjusted
// Use time approach
/************************************************************************/

bool t_StabSolver::adjustLocal
(t_WCharsLoc &a_wave_chars, t_StabSolver::t_MODE a_mode){
	// parameters
	const t_Complex d_arg = _params.AdjustStep;
	const int max_iters = _params.AdjustMaxIter;
	const double tol = _params.AdjustTol;
	// ensure that there is a base residual
	solve(a_wave_chars);
	t_WCharsLoc backup = a_wave_chars;
	t_Complex *choose_arg;
	if (a_mode==A_MODE){
		choose_arg = &(a_wave_chars.a);
	}else{
		if (a_mode==B_MODE){
			choose_arg = &(a_wave_chars.b);
		}else{
			choose_arg = &(a_wave_chars.w);
		};
	};
	t_Complex& arg = *choose_arg;
	t_Complex arg_base = arg;
	t_Complex res_base = a_wave_chars.resid;
	t_Complex deriv;
	for (int i=0; i<max_iters; i++){
		// check if we are converged
		if (smat::norm(a_wave_chars.resid)<tol){
			//std::wcout<<_T("Adjust Converged")<<res_base<<std::endl;
			return true;
		};
		// newton iteration
		arg+=d_arg;
		solve(a_wave_chars);
		deriv = d_arg/(a_wave_chars.resid - res_base);
		arg = arg_base - res_base*deriv;
		solve(a_wave_chars);
		arg_base = arg;
		res_base = a_wave_chars.resid;
	};
	wxLogMessage(_T("Local Search Error: no convergence\n"));
	// TODO: debug!
	ssuGENTHROW(_T("Test"));
	a_wave_chars = backup;
	return false;
};

std::vector<t_WCharsLoc> t_StabSolver::filterInitWaves
(const std::vector<t_WCharsLoc>& all_initials){
	const double close_tol = 0.1;
	std::vector<t_WCharsLoc> good_initials;
	std::vector<t_WCharsLoc>::const_iterator iter;
	t_WCharsLoc init_wave;
	for (iter = all_initials.begin(); iter<all_initials.end(); iter++){
		const t_WCharsLoc& init_wave = *iter;
		t_WCharsLoc adj_wave = init_wave;
		try{
			adjustLocal(adj_wave, t_MODE::W_MODE);
			if (abs(adj_wave.w.real()-init_wave.w.real())<
				close_tol*abs(adj_wave.w.real()+init_wave.w.real())){
					good_initials.push_back(adj_wave);
			};
		}catch(t_UnPhysWave){
			// oh shit!
		};
	};
	return good_initials;
};


/************************************************************************/   
//Search for nearest eigen solution
// with keeping Re(w)=const
/************************************************************************/
bool t_StabSolver::getEigenWFixed
(double wr_fixed, t_WCharsLoc& wave_chars, t_MODE mode){

	bool ok=true;

	t_WCharsLoc backup = wave_chars;
	t_Complex *pVg, *pArg;
	if (mode==A_MODE){
		pArg = &wave_chars.a;
		pVg  = &wave_chars.vga;
	}else{
		if (mode==B_MODE){
			pArg = &wave_chars.b;
			pVg  = &wave_chars.vgb;
		}else{
			throw t_WrongMode();
		}
	}
	double resid=1.0;
	double dwr, darg;
	t_Complex vg;
	// IMPORTANT TODO: do we think Vg change is small??
	// or is it a mistake - Vg is not recalculated!!!
	ok = ok && adjustLocal(wave_chars, W_MODE);
	calcGroupVelocity(wave_chars);
	int n_iter=0;
	do{
		dwr = wr_fixed - wave_chars.w.real();
		vg = *pVg;
		darg = (dwr/vg).real();
		*pArg+= darg;
		wave_chars.w+= dwr;
		ok = ok && adjustLocal(wave_chars, W_MODE);
		resid = abs((dwr)/wr_fixed);
		if (n_iter++>_params.AdjustMaxIter){
			wxLogMessage(_T("In GetEigenWFixed: no convergence"));
			return false;
		};
	} while (resid>=_params.AdjustTol);

	return ok;
};

/************************************************************************/   
// search for a wave chars with conditions provided
// TODO: this is to be developed futher 
/************************************************************************/

bool t_StabSolver::searchWave
(t_WCharsLoc& wchars, stab::t_LSCond cond, stab::t_TaskTreat task_mode)
{

	switch (cond.get_mode())
	{
	case stab::t_LSCond::W_FIXED:
		// TODO: maybe in some cases B_MODE, how to switch ?
		if (task_mode==stab::t_TaskTreat::SPAT){
			wxString msg(_T("PF.LS: W_FIXED implemented only for TIME approach"));
			wxLogError(msg); return false;
		}else{
			return getEigenWFixed(cond.wchars.w.real(), wchars, t_StabSolver::A_MODE);
		} ;
		break;

	case (stab::t_LSCond::A_FIXED|stab::t_LSCond::B_FIXED):
		return adjustLocal(wchars, W_MODE);
		//return adjustLocal_Grad(wchars, cond);

	case (stab::t_LSCond::A_FIXED|stab::t_LSCond::W_FIXED):
		return adjustLocal(wchars, B_MODE);
		//return adjustLocal_Grad(wchars, cond);

	case (stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED):
		return adjustLocal(wchars, A_MODE);
		//return adjustLocal_Grad(wchars, cond);

	default:
		throw t_NotImplemented();
	}

	return false;

};

void t_StabSolver::searchMaxWave(t_WCharsLoc& wchars, stab::t_LSCond cond, stab::t_TaskTreat task_mode){


	switch (cond.get_mode())
	{

	case stab::t_LSCond::FREE:
		if (task_mode==stab::t_TaskTreat::TIME){
			wchars = _getMaxInstabTime_Grad(wchars);
		}else{
			wchars = _getMaxInstabSpat_Grad(wchars);
		}
		return;

	case stab::t_LSCond::W_FIXED:
		if (task_mode==stab::t_TaskTreat::TIME){
			wchars.w.real(cond.wchars.w.real());
			wchars = _getStationaryMaxInstabTime(wchars);
		}else{
			//_getMaxInstabAlpha_Spat(wchars);
			throw t_NotImplemented();
		}
		return;

	case stab::t_LSCond::A_FIXED:
		throw t_NotImplemented();

	case stab::t_LSCond::B_FIXED:
		throw t_NotImplemented();

	default:
		throw t_NotImplemented();
	}
}

void t_StabSolver::calcAiDbDerivs(t_WCharsLoc& wave, double& dai_dbr, double& d2ai_dbr2, double a_darg) {

	stab::t_LSCond srch_cond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED);

	t_WCharsLoc base_wave = wave;

	searchWave(base_wave, srch_cond, stab::t_TaskTreat::SPAT);

	double c_val = base_wave.a.imag();

	t_WCharsLoc rgt_wave = base_wave;
	rgt_wave.b = base_wave.b.real() + a_darg; rgt_wave.resid = 1.0;
	searchWave(rgt_wave, srch_cond, stab::t_TaskTreat::SPAT);
	double r_val = rgt_wave.a.imag();

	t_WCharsLoc lft_wave = base_wave;
	lft_wave.b = base_wave.b.real() - a_darg; rgt_wave.resid = 1.0;
	searchWave(lft_wave, srch_cond, stab::t_TaskTreat::SPAT);
	double l_val = lft_wave.a.imag();

	dai_dbr = (r_val - l_val) / (2.0*a_darg);

	d2ai_dbr2 = (r_val - 2.0*c_val + l_val) / (a_darg*a_darg);

	wxLogMessage(_T("Fun=%f; Deriv=%f"), dai_dbr, d2ai_dbr2);

	wave = base_wave;

};

bool t_StabSolver::searchMaxAiSpat(const t_WCharsLoc& w_init, t_WCharsLoc& w_max) {

	bool ok = true;
	const int max_iters = 100;

	// TODO: emprics with d_rg, move to config when tested
	// not working with small d_arg (like 1.0e-07)
	const double d_arg = 1.0e-03;
	const double tol = 1.0e-06;

	w_max = w_init;

	double fun, fun_deriv;

	for (int i = 0; i<max_iters; i++) {

		calcAiDbDerivs(w_max, fun, fun_deriv, d_arg);

		// check if we are converged
		if (abs(fun)<tol) {
			std::wcout << _T("search dai_db=0 Converged\n") << w_max << std::endl;
			return true;
		};

		t_WCharsLoc base_wave = w_max;

		w_max.b = base_wave.b - fun / fun_deriv;
	};

	wxLogMessage(_T("Error: search max ai vs br - no convergence\n"));

	w_max = w_init;
	return false;


};

std::vector<t_WCharsLoc> t_StabSolver::filter_gs_waves_spat(const std::vector<t_WCharsLoc> wcands, stab::t_LSCond cond){

	std::vector<t_WCharsLoc> ret_waves;
	t_WCharsLoc cur_wave;

	for (int i=0; i<wcands.size(); i++){

		cur_wave = wcands[i];

		//std::wcout<<_T("GS Init:")<<cur_wave;
		bool good_init;
		try
		{

			wxLogMessage(_T("=============GS============="));

			wxLogMessage(_T("GS Estimate:%s"), &cur_wave.to_wstr()[0]);

			// first make raw estimate - good wave or not
			if (!( stab::check_wchars_c_phase(cur_wave) && 
				stab::check_wchars_increment(cur_wave)
				)) continue;

			good_init = searchWave(cur_wave, cond, stab::t_TaskTreat::SPAT);

			if (good_init && cur_wave.a.real() >= 0)
				wxLogMessage(_T("Discrete mode candidate:%s"), &cur_wave.to_wstr()[0]);

			if (good_init && stab::check_wchars_c_phase(cur_wave) && 
				stab::check_wchars_increment(cur_wave)){

				cur_wave.set_scales(get_stab_scales());

				t_WCharsLocDim dim_wave = cur_wave.make_dim();

				// TODO: nice checking that wave is physical
				if ( stab::check_wchars_c_phase(cur_wave) && 
					this->checkWCharsByGroupV(cur_wave) &&
					stab::check_wchars_increment(cur_wave)
					){

					wxLogMessage(_T("Checks for group velocity & increment: ok"));

					ret_waves.push_back(cur_wave);
				}else
					wxLogMessage(_T("Checks for group velocity & increment: failed (!)"));
			}
		}
		catch (...)
		{
			continue;

		}

	}

	return ret_waves;

}

std::vector<t_WCharsLoc> t_StabSolver::filter_gs_waves_time(const std::vector<t_WCharsLoc> wcands, stab::t_LSCond cond){

	std::vector<t_WCharsLoc> ret_waves;
	t_WCharsLoc cur_wave;

	for (int i=0; i<wcands.size(); i++){

		cur_wave = wcands[i];

		std::wcout<<_T("GS Init:")<<cur_wave;
		bool good_init;
		try
		{

			good_init = searchWave(cur_wave, cond, stab::t_TaskTreat::TIME);

			if (good_init && cur_wave.a.real()>=0)
				std::wcout<<_T("Discrete mode found:")<<cur_wave;

			if (good_init && cur_wave.a.real()>0 && cur_wave.w.imag()>0.0){

				cur_wave.set_scales(get_stab_scales());

				t_WCharsLocDim dim_wave = cur_wave.make_dim();

				// TODO: nice checking that wave is physical
				if ( stab::check_wchars_c_phase(cur_wave) && 
					this->checkWCharsByGroupV(cur_wave) &&
					stab::check_wchars_increment(cur_wave)
					){

						std::wcout<<_T("Instab found:")<<cur_wave;

						ret_waves.push_back(cur_wave);
				}
			}
		}
		catch (...)
		{
			continue;

		}

	}

	return ret_waves;

}

void t_StabSolver::getAmpFuncs(std::vector<t_VecCmplx>& amp_funcs){

	if (amp_funcs.size()!=getNNodes()) wxLogError(_T("GetAmpFuncs: wrong size of input vec"));

	std::vector<t_MatCmplx> solutions = _math_solver.reconstruct();

	int nvecs = getTaskDim();
	t_VecCmplx wall_coefs(4);

	_calcResidual(&wall_coefs[0]);

	for (int j=0; j<getNNodes(); j++){

		// amp_funcs = solutions[j]*wall_coefs;
		matrix::base::mat_mul(solutions[j], wall_coefs, amp_funcs[j]);

	}
}
// get index of amplitude function in amp funcs array by the name of function
int t_StabSolver::getFuncIndInAmpFuncs(char func_name) {

	if (func_name == 'u') return 0;
	if (func_name == 'v') return 2;
	if (func_name == 'p') return 3;
	if (func_name == 't') return 4;
	if (func_name == 'w') return 6;

	wxLogError(_T("Wrong func_name in t_StabSolver::getFuncIndInAmpFuncs: only u,v,w,p,t supported!"));
	return 0;
}

// compute <H1x, z>
// notation is according to AIAA-2002-2846 : 
// H1 = -i*dH0/da
// x - amp fun of direct problem (eigenval wchars_A_in)
// z - amp fun of conjugate problem (eigenval wchars_B_in)
t_Complex t_StabSolver::calcScalarProd_H1(
	const t_WCharsLoc& wchars_A_in, const t_WCharsLoc& wchars_B_in,
	std::vector<t_VecCmplx>* dns_vec_ptr){

	// important - first solve conjugate problem, then direct
	// because matrix of scalar prod should be computed for a direct task:
	// H1 = -i*dH0/da - for direct problem

	int nnodes = _math_solver.getNNodes();

	t_WCharsLoc wchars_A, wchars_B;

	wchars_A = wchars_A_in;
	wchars_B = wchars_B_in;

	if ((wchars_A.get_treat()!=stab::t_TaskTreat::SPAT)||
		(wchars_B.get_treat()!=stab::t_TaskTreat::SPAT)){
		wxLogError(_T("CalcScalProd Error: spat waves only, update needed for time waves!"));
		return HUGE_VAL;
	}

	std::vector<t_VecCmplx> sol_dir_A(nnodes, STAB_MATRIX_DIM);
	std::vector<t_VecCmplx> sol_con_B(nnodes, STAB_MATRIX_DIM);

	// conjugate problem

	setLSMode(stab::t_LSMode(stab::t_LSMode::CONJUGATE|stab::t_LSMode::ASYM_HOMOGEN));

	searchWave(wchars_B, stab::t_LSCond(stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED), 
		stab::t_TaskTreat::SPAT);


	wxLogMessage(_T("Wave chars B (conj):%s"), &(wchars_B.to_wstr()[0]));

	wxLogMessage(_T("Conjugate problem residual:%f"), smat::norm(wchars_B.resid));

	getAmpFuncs(sol_con_B);

	// direct problem

	setLSMode(stab::t_LSMode(stab::t_LSMode::DIRECT|stab::t_LSMode::ASYM_HOMOGEN));

	searchWave(wchars_A, stab::t_LSCond(stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED), 
		stab::t_TaskTreat::SPAT);

	wxLogMessage(_T("Wave chars A (direct):%s"), &(wchars_A.to_wstr()[0]));

	wxLogMessage(_T("Direct problem residual:%f"), smat::norm(wchars_A.resid));

	// if amplitude vector is provided from dns, use it as vector for wchars_A
	// otherwise use both vectors from LST
	if (dns_vec_ptr==NULL)
	{
		getAmpFuncs(sol_dir_A);

	} 
	else
	{
		if (dns_vec_ptr->size()!=getNNodes()){

			wxLogError(_T("Error: in calc scalar product - DNS vec size doesn't match LST vec size"));
			return 1.0e+09;

		}
			
		sol_dir_A = *dns_vec_ptr;
	}

	// amp funcs are calculated, compute scalar product
	// stab context is set for "direct" wave wchars_A

	std::vector<t_Complex> fun(nnodes, 0.0);
	std::vector<double>arg(nnodes, 0.0);

	t_VecCmplx v1(STAB_MATRIX_DIM), v2(STAB_MATRIX_DIM);

	// use order from wall to outer region !
	// solution order is from outer region to wall

	for (int i=0; i<nnodes; i++){

		int ind_r = nnodes - 1 - i;
		arg[i] = _math_solver.varRange[ind_r];

		const t_ProfRec& rec = _profStab.get_rec(i);
		_setScalProdMatrix_H1(rec);
		//matrix::base::mat_mul(_scal_prod_matrix, sol_dir[ind_r], v1);
		// using plain product because sol_conj is already a conjugated vector of conjugate task
		//fun[i] = vector::plain_prod<t_Complex, t_Complex>(v1, sol_conj[ind_r]);

//		DO 12 I=1,MF
//			DO 12 J=1,MF
//			12 GG(II)=GG(II)+H2(I,J)*A0(J,II)*DZ(I,II)

		fun[i]=0.0;
		for (int j=0; j<STAB_MATRIX_DIM; j++)
			for (int k=0; k<STAB_MATRIX_DIM; k++)
				fun[i] = fun[i] + _scal_prod_matrix_H1[k][j]*sol_dir_A[ind_r][k]*sol_con_B[ind_r][j];


	}

	//return smat::fun_integrate_simp4_uniform(arg, fun);

	return smat::fun_integrate(arg, fun);

}
// compute <H1*fun_direct, fun_conj>
// H1 = -i*dH/da
// NB:do not forget to set context of loc solver before
t_Complex t_StabSolver::calcScalarProd_H1(std::vector<t_VecCmplx>& fun_direct, std::vector<t_VecCmplx>& fun_conj) {

	int nnodes = _math_solver.getNNodes();

	if (nnodes != fun_direct.size() || nnodes != fun_conj.size())
		wxLogError(_T("t_StabSolver::calcScalarProd_H1: wrong size of amplitude functions"));

	std::vector<t_Complex> fun(nnodes, 0.0);
	std::vector<double>arg(nnodes, 0.0);

	for (int i = 0; i<nnodes; i++) {

		int ind_r = nnodes - 1 - i;
		arg[i] = _math_solver.varRange[ind_r];

		const t_ProfRec& rec = _profStab.get_rec(i);
		_setScalProdMatrix_H1(rec);

		fun[i] = 0.0;
		for (int j = 0; j<STAB_MATRIX_DIM; j++)
			for (int k = 0; k<STAB_MATRIX_DIM; k++)
				// using plain product because sol_conj is already a conjugated vector of conjugate task (!)
				fun[i] = fun[i] + _scal_prod_matrix_H1[k][j] * fun_direct[ind_r][k] * fun_conj[ind_r][j];


	}

	return smat::fun_integrate(arg, fun);

}

// compute <HW*fun_direct, fun_conj>
// HW = i*dH/dw
// NB:do not forget to set context of loc solver before
t_Complex t_StabSolver::calcScalarProd_HW(std::vector<t_VecCmplx>& fun_direct, std::vector<t_VecCmplx>& fun_conj) {

	int nnodes = _math_solver.getNNodes();

	if (nnodes != fun_direct.size() || nnodes != fun_conj.size())
		wxLogError(_T("t_StabSolver::calcScalarProd_H1: wrong size of amplitude functions"));

	std::vector<t_Complex> fun(nnodes, 0.0);
	std::vector<double>arg(nnodes, 0.0);

	for (int i = 0; i<nnodes; i++) {

		int ind_r = nnodes - 1 - i;
		arg[i] = _math_solver.varRange[ind_r];

		const t_ProfRec& rec = _profStab.get_rec(i);
		_setScalProdMatrix_H1(rec);

		fun[i] = 0.0;
		for (int j = 0; j<STAB_MATRIX_DIM; j++)
			for (int k = 0; k<STAB_MATRIX_DIM; k++)
				// using plain product because sol_conj is already a conjugated vector of conjugate task (!)
				fun[i] = fun[i] + _scal_prod_matrix_HW[k][j] * fun_direct[ind_r][k] * fun_conj[ind_r][j];


	}

	return smat::fun_integrate(arg, fun);

}
/************************************************************************/
// Calculate <H2*fun_direct, fun_conj>
// context must be set before
/************************************************************************/
t_Complex t_StabSolver::calcScalarProd_H2(std::vector<t_VecCmplx>& fun_direct, std::vector<t_VecCmplx>& fun_conj) {

	int nnodes = _math_solver.getNNodes();

	if (nnodes != fun_direct.size() || nnodes != fun_conj.size())
		wxLogError(_T("t_StabSolver::calcScalarProd_H1: wrong size of amplitude functions"));

	std::vector<t_Complex> fun(nnodes, 0.0);
	std::vector<double>arg(nnodes, 0.0);

	t_VecCmplx H2z(STAB_MATRIX_DIM);

	for (int i = 0; i<nnodes; i++) {

		int ind_r = nnodes - 1 - i;
		arg[i] = _math_solver.varRange[ind_r];

		_calc_H2z(i, fun_direct, H2z);

		fun[i] = 0.0;
		//for (int j = 0; j<STAB_MATRIX_DIM; j++)
		//	for (int k = 0; k<STAB_MATRIX_DIM; k++)
				// using plain product because sol_conj is already a conjugated vector of conjugate task (!)
		//		fun[i] = fun[i] + _scal_prod_matrix_H2[k][j] * fun_direct[ind_r][k] * fun_conj[ind_r][j];

		for (int j = 0; j < STAB_MATRIX_DIM; j++)
			fun[i] = fun[i] + H2z[j] * fun_conj[ind_r][j];

	}

	return smat::fun_integrate(arg, fun);

}

/************************************************************************/
// Search for lower and upper points of neutral curve at given xyz
// for now beta = fixed and seek for a = a(w), im(a)=0, lower and upper
// starting wchars must be unstable (usually most unstable wave from global search)
/************************************************************************/
void t_StabSolver::calcNeutPoints(const mf::t_GeomPoint& xyz, const t_WCharsLoc& wave_start, 
	t_WCharsLoc& wave_lower, t_WCharsLoc& wave_upper) {

	setContext(xyz);

	int n_max_steps = 100000;

	stab::t_LSCond srch_cond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED);

	const double eps = 0.0025;

	double dw = 0.0;
	double dw_calc = 0.0;

	t_WCharsLoc wave = wave_start;
	t_WCharsLoc wave_prv;

	// search upper point of neutral curve
	while (wave.a.imag() < 0.0) {

		dw_calc = eps*wave.w.real();
		// when w is nearly zero relative step is too small, use absolute step
		dw = (abs(dw_calc) > 1.0e-05) ? dw_calc : 1.0e-05;

		wave_prv = wave;

		wave.w += eps*wave.w.real();
		searchWave(wave, srch_cond, stab::t_TaskTreat::SPAT);

		// debug
		wxLogMessage(_T("w=%lf, -Im(a)=%lf"), wave.w.real(), -1.0*wave.a.imag());

	}
	double slope = (wave.a.imag() - wave_prv.a.imag()) / (wave.w.real() - wave_prv.w.real());
	double darg = -1.0/slope*wave_prv.a.imag();
	double val = wave_prv.a.imag() + slope*darg;

	wxLogMessage(_T("Left point:w=%lf, -Im(a)=%lf"), wave_prv.w.real(), -1.0*wave_prv.a.imag());
	wxLogMessage(_T("right point:w=%lf, -Im(a)=%lf"), wave.w.real(), -1.0*wave.a.imag());
	wxLogMessage(_T("Zero:w=%lf, -Im(a)=%lf"), wave_prv.w.real() + darg, val);

	wave_upper = wave_prv;
	// linear interpolation for frequency
	wave_upper.w = t_Complex(wave_prv.w.real() + darg, val);

	wave_upper.set_scales(get_stab_scales());

	wave = wave_start;

	// search upper point of neutral curve
	while (wave.a.imag() < 0.0) {

		dw_calc = eps*wave.w.real();
		// when w is nearly zero relative step is too small, use absolute step
		dw = (abs(dw_calc) > 1.0e-05) ? dw_calc : 1.0e-05;

		wave_prv = wave;

		wave.w -= eps*wave.w.real();
		searchWave(wave, srch_cond, stab::t_TaskTreat::SPAT);

		// debug
		wxLogMessage(_T("w=%lf, -Im(a)=%lf"), wave.w.real(), -1.0*wave.a.imag());

	}
	slope = (wave.a.imag() - wave_prv.a.imag()) / (wave.w.real() - wave_prv.w.real());
	darg = -1.0 / slope*wave_prv.a.imag();
	val = wave_prv.a.imag() + slope*darg;

	wxLogMessage(_T("Left point:w=%lf, -Im(a)=%lf"), wave_prv.w.real(), -1.0*wave_prv.a.imag());
	wxLogMessage(_T("right point:w=%lf, -Im(a)=%lf"), wave.w.real(), -1.0*wave.a.imag());
	wxLogMessage(_T("Zero:w=%lf, -Im(a)=%lf"), wave_prv.w.real() + darg, val);

	wave_lower = wave_prv;
	// linear interpolation for frequency
	wave_lower.w = t_Complex(wave_prv.w.real() + darg, val);

	wave_lower.set_scales(get_stab_scales());

};

// calculate dv_dy and dp_dy for an amplitude function
// i is index in numeration from outer to wall : i = 0 - outer rec, i = nnodes - 1 - wall rec
// second order approximation
void t_StabSolver::_calc_amp_fun_dv_dy_dp_dy(int i, std::vector<t_VecCmplx>& amp_fun, t_Complex& dv_dy, t_Complex& dp_dy) const {

	double dx_inv = 1.0/(_math_solver.varRange[1] - _math_solver.varRange[0]);
	t_Complex f0, f1, f2;
	t_Complex f_r[2];

	// indices of v and p are 2 and 3
	// (u, u', v, p ,t, t', w, w')
	if (i == 0) {
		for (int k = 0; k < 2; k++) {
			f0 = amp_fun[0][2+k];
			f1 = amp_fun[1][2+k];
			f2 = amp_fun[2][2+k];

			f_r[k] = 0.5*dx_inv*(-3.0 * f0 + 4.0*f1 - f2);
		}
		dv_dy = f_r[0];
		dp_dy = f_r[1];
		return;
	}
	int imax = _math_solver.getNNodes() - 1;
	if (i == imax) {
		for (int k = 0; k < 2; k++) {
			f0 = amp_fun[imax - 0][2+k];
			f1 = amp_fun[imax - 1][2+k];
			f2 = amp_fun[imax - 2][2+k];

			f_r[k] = 0.5*dx_inv*(3.0 * f0 - 4.0*f1 + f2);
		}
		dv_dy = f_r[0];
		dp_dy = f_r[1];
		return;
	}

	for (int k = 0; k < 2; k++) {
		f0 = amp_fun[i - 1][2+k];
		f1 = amp_fun[i + 1][2+k];

		f_r[k] = 0.5*dx_inv*(f1 - f0);
	}
	dv_dy = f_r[0];
	dp_dy = f_r[1];
	return;

};

void t_StabSolver::normalizeAmpFuncsByPressureAtWall(std::vector<t_VecCmplx>& amp_funcs){

	int nnodes = amp_funcs.size();

	t_Complex pw = amp_funcs[nnodes - 1][3];

	for (int i = 0; i < nnodes; i++)
		for (int j = 0; j < STAB_MATRIX_DIM; j++)
			amp_funcs[i][j] /= pw;

}
