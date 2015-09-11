#include "stdafx.h"

#include "StabSolver.h"
#include "common_data.h"
#include "log.h"

using namespace hsstab;
using namespace hsstab::cmpnts;

using namespace pf;

t_StabSolver::t_StabSolver(const mf::t_DomainBase& a_rFldNS):
_rFldNS(a_rFldNS), _profStab(), _math_solver(), 
_stab_matrix(STAB_MATRIX_DIM),_scal_prod_matrix(STAB_MATRIX_DIM), _params(){

	_math_solver._pStab_solver = this;

	_ls_mode.set_defaults();

};

void t_StabSolver::init(const hsstab::TPlugin& g_plug){

	const hsstab::TPluginParamsGroup& g = g_plug.get_settings_grp_const("");

	ls::_init_ortho_ls_base_params(_params, g);

	//_init();	// init solver

}

const t_StabScales& t_StabSolver::get_stab_scales() const{
	return _profStab.scales();
};

// private, to be used in 3D context

/*const t_SqMatrix& t_StabSolver::_getStabMatrix3D() const{
	return _stab_matrix;
}*/

// old - can be used only in local spectral procedures
/*void t_StabSolver::setMatSplitByW3D_88(const double& a_y, t_SqMatrix& mat_no_w, t_SqMatrix& mat_w) const{
	t_CompVal imagUnity(0.0, 1.0);

	const double& stabRe = _profStab.stabRe;
	const double& Me = _profStab.Me;
	const double gMaMa = MF_Field::Gamma*Me*Me;
	const double g_1MaMa = (MF_Field::Gamma-1.0)*Me*Me;

	const double u = _profStab.getValue(a_y, _profStab.u);
	const double u1 = _profStab.getValue(a_y, _profStab.u1);
	const double u2 = _profStab.getValue(a_y, _profStab.u2);

	const double w = _profStab.getValue(a_y, _profStab.w);
	const double w1 = _profStab.getValue(a_y, _profStab.w1);
	const double w2 = _profStab.getValue(a_y, _profStab.w2);

	const double t = _profStab.getValue(a_y, _profStab.t);
	const double t1 = _profStab.getValue(a_y, _profStab.t1);
	const double t2 = _profStab.getValue(a_y, _profStab.t2);

	const double mu = _profStab.getValue(a_y, _profStab.mu);
	const double mu1 = _profStab.getValue(a_y, _profStab.mu1);
	const double mu2 = _profStab.getValue(a_y, _profStab.mu2);

	const double inv_t = 1.0/t;
	const double inv_mu = 1.0/mu;

	const double vCoefL = 2.0/3.0*(0.0+2.0); //TODO: second Visc coef instead of 0.0
	const double vCoefS = 2.0/3.0*(0.0-1.0); // 
	const double vCoefM = 1.0+vCoefS;

	const t_CompVal& alpha = _waveChars.a;
	const t_CompVal& beta = _waveChars.b;

	const t_CompVal dEicon_no_W = alpha*u + beta*w ;
	const t_CompVal dEicon_W = -1.0;
	const double dMu2 = mu2*t1*u1 + mu1*u2;
	const t_CompVal xi = 1.0/(stabRe*inv_mu+imagUnity*vCoefL+gMaMa*dEicon_no_W);
// set no_W matrix
// of size 5rowsx8cols
	// first <--> "2"
	mat_no_w[0][0] = imagUnity*stabRe*inv_t*inv_mu*dEicon_no_W +
		pow(alpha,2) + pow(beta,2);

	mat_no_w[1][0] =  - inv_mu*mu1*t1;

	mat_no_w[2][0] = stabRe*inv_t*inv_mu*u1 - 
		imagUnity*alpha*(inv_mu*mu1*t1+vCoefM*inv_t*t1);

	mat_no_w[3][0] = alpha*(imagUnity*stabRe*inv_mu -
		vCoefM*gMaMa*dEicon_no_W);

	mat_no_w[4][0] = vCoefM*alpha*inv_t*dEicon_no_W - inv_mu*dMu2;

	mat_no_w[5][0] = -inv_mu*mu1*u1;
	// second <--> "3"
	mat_no_w[0][1] = -imagUnity*alpha;

	mat_no_w[2][1] = inv_t*t1;

	mat_no_w[3][1] = -imagUnity*gMaMa*dEicon_no_W;

	mat_no_w[4][1] = imagUnity*inv_t*dEicon_no_W;

	mat_no_w[6][1] = -imagUnity*beta;
	// third <--> "4"
	mat_no_w[0][2] = -imagUnity*xi*alpha*t1*
		(2.0*inv_mu*mu1 + vCoefL*inv_t);

	mat_no_w[1][2] = -imagUnity*xi*alpha;

	mat_no_w[2][2] = xi*(-alpha*alpha-beta*beta+
		vCoefL*inv_t*inv_mu*mu1*t1*t1+
		vCoefL*inv_t*t2-
		imagUnity*stabRe*inv_t*inv_mu*dEicon_no_W);

	mat_no_w[3][2] = -imagUnity*xi*vCoefL*gMaMa*
		(
		(inv_mu*mu1*t1+inv_t*t1)*dEicon_no_W+
		alpha*u1+beta*w1
		);

	mat_no_w[4][2] = imagUnity*xi*
		(
		(inv_mu*mu1+vCoefL*inv_t)*(alpha*u1+beta*w1)+
		vCoefL*inv_t*inv_mu*mu1*t1*dEicon_no_W
		);

	mat_no_w[5][2] = imagUnity*xi*vCoefL*inv_t*dEicon_no_W;

	mat_no_w[6][2] = -imagUnity*xi*beta*
		(2.0*inv_mu*mu1*t1 + vCoefL*inv_t*t1);

	mat_no_w[7][2] = -imagUnity*xi*beta;
	// fourth <--> "6"
	mat_no_w[1][3] = -2.0*MF_Field::Pr*g_1MaMa*u1;

	mat_no_w[2][3] = MF_Field::Pr*
		(
		stabRe*inv_t*inv_mu*t1 - 
		2.0*imagUnity*g_1MaMa*(alpha*u1+beta*w1)
		);

	mat_no_w[3][3] = -imagUnity*stabRe*MF_Field::Pr*inv_mu*g_1MaMa*dEicon_no_W;

	mat_no_w[4][3] = imagUnity*stabRe*MF_Field::Pr*inv_t*inv_mu*dEicon_no_W +
		alpha*alpha + beta*beta - 
		g_1MaMa*MF_Field::Pr*inv_mu*mu1*(u1*u1+w1*w1)-
		inv_mu*(mu2*t1*t1+mu1*t2);
	mat_no_w[5][3] = -2.0*inv_mu*mu1*t1;

	mat_no_w[7][3] = -2.0*MF_Field::Pr*g_1MaMa*w1;
	// fifth <--> "8"
	mat_no_w[2][4] = -imagUnity*beta*(inv_mu*mu1*t1+vCoefM*inv_t*t1)+
		stabRe*inv_mu*inv_t*w1;

	mat_no_w[3][4] = imagUnity*stabRe*beta*inv_mu-
		vCoefM*beta*gMaMa*dEicon_no_W;

	mat_no_w[4][4] = vCoefM*beta*inv_t*dEicon_no_W-
		inv_mu*(mu2*t1*w1+mu1*w2);

	mat_no_w[5][4] = -inv_mu*mu1*w1;

	mat_no_w[6][4] = imagUnity*stabRe*inv_t*inv_mu*dEicon_no_W +
		alpha*alpha + beta*beta;

	mat_no_w[7][4] = -inv_mu*mu1*t1;
// set W matrix
// of the same sizes
// 5rows x 8 cols
	// first <--> part of "2", 3 non-zeros
	mat_w[0][0] = imagUnity*stabRe*inv_t*inv_mu*dEicon_W;

	mat_w[3][0] = alpha*( -	vCoefM*gMaMa*dEicon_W);

	mat_w[4][0] = vCoefM*alpha*inv_t*dEicon_W;

	// second <--> part of "3", 2 non-zeros

	mat_w[3][1] = -imagUnity*gMaMa*dEicon_W;

	mat_w[4][1] = imagUnity*inv_t*dEicon_W;

	// third <--> "4", 4 non-zeros, a lot)

	mat_w[2][2] = xi*(-imagUnity*stabRe*inv_t*inv_mu*dEicon_W);

	mat_w[3][2] = -imagUnity*xi*vCoefL*gMaMa*
		(
		(inv_mu*mu1*t1+inv_t*t1)*dEicon_W
		);

	mat_w[4][2] = imagUnity*xi*
		(
		vCoefL*inv_t*inv_mu*mu1*t1*dEicon_W
		);

	mat_w[5][2] = imagUnity*xi*vCoefL*inv_t*dEicon_W;

	// fourth <--> part of "6", 2 non-zeros

	mat_w[3][3] = -imagUnity*stabRe*MF_Field::Pr*inv_mu*g_1MaMa*dEicon_W;

	mat_w[4][3] = imagUnity*stabRe*MF_Field::Pr*inv_t*inv_mu*dEicon_W;

	// fifth <--> "8", 3 non-zeros

	mat_w[3][4] = -vCoefM*beta*gMaMa*dEicon_W;

	mat_w[4][4] = vCoefM*beta*inv_t*dEicon_W;

	mat_w[6][4] = imagUnity*stabRe*inv_t*inv_mu*dEicon_W;

	return;
};*/

void t_StabSolver::_formRHS3D(const double& a_y, const t_VecCmplx& a_vars, t_VecCmplx& dest){
	_setStabMatrix3D(a_y);
	matrix::base::mat_mul<t_Complex, t_Complex>(_stab_matrix, a_vars, dest);
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

	t_ProfileNS profNS(_rFldNS);

	mf::t_ProfDataCfg prof_cfg;
	prof_cfg.ThickCoef = _params.ThickCoef;

	// NNodes can be used in INTERPOLATE types of initialization
	prof_cfg.NNodes = _params.NNodes;

	switch (_params.NSProfInit)
	{
	case (blp::t_NSInit::EXTRACT):
		profNS.initialize(a_xyz, prof_cfg, blp::t_NSInit::EXTRACT);
		break;
	case (blp::t_NSInit::INTERPOLATE):
		profNS.initialize(a_xyz, prof_cfg, blp::t_NSInit::INTERPOLATE);
		break;
	default:
		wxString msg(_T("PF.LocSearch: ProfNS Initialization type not supported"));
		wxLogError(msg); ssuGENTHROW(msg);
	}

	_profStab.initialize(profNS, nnodes_stab);

	if (_stab_matrix.nCols()!=STAB_MATRIX_DIM){
		ssuTHROW(t_GenException, 
		_T("PF.LocSearch Error:Wrong size of stab matrix while setting context"));
	}

	for (int j=0; j<nnodes_stab; j++){
		_math_solver.varRange[j] = _profStab.get_y(nnodes_stab-1-j);
	}

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

	if (vgr>1.0){

		wxLogMessage(_T("WChars Check: group velo is big - unphysical wave"));
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

void t_StabSolver::setInitWaves(const std::vector<t_WCharsLoc>& a_inits){
	this->_initWaves = a_inits;
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


t_WCharsLoc t_StabSolver::getMaxWave(const std::vector<t_WCharsLoc>& a_inits, int a_nnodesStab){
    // TODO:to implement somehow "default" parameter
	// TODO: nnodes_stab not used here?
	int nnodes_stab = (a_nnodesStab==0) ? _params.NNodes : a_nnodesStab;

	setInitWaves(a_inits);
	t_WCharsLoc res_wave;
	//ensure
	res_wave.w = 0.0;
	t_WCharsLoc cur_wave;
	for (int i=0; i<_initWaves.size(); i++){
		cur_wave = _getMaxInstabTime_Grad(_initWaves[i]);
		if (cur_wave.w.imag()>res_wave.w.imag()){
			res_wave = cur_wave;
		}
	}
	return res_wave;
}

std::vector<t_WCharsLoc> t_StabSolver::filter_gs_waves_spat(const std::vector<t_WCharsLoc> wcands, stab::t_LSCond cond){

	std::vector<t_WCharsLoc> ret_waves;
	t_WCharsLoc cur_wave;

	for (int i=0; i<wcands.size(); i++){

		cur_wave = wcands[i];

		std::wcout<<_T("GS Init:")<<cur_wave;
		bool good_init;
		try
		{

			// first make raw estimate - good wave or not
			if (!( stab::check_wchars_c_phase(cur_wave) && 
				stab::check_wchars_increment(cur_wave)
				)) continue;

			good_init = searchWave(cur_wave, cond, stab::t_TaskTreat::SPAT);

			if (good_init && cur_wave.a.real()>=0)
				std::wcout<<_T("Discrete mode found:")<<cur_wave;

			if (good_init && stab::check_wchars_c_phase(cur_wave) && 
				stab::check_wchars_increment(cur_wave)){

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
