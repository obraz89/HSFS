#include "stdafx.h"

#include "common_data.h"

#include "EigenGS_Spatial.h"

#include <iostream>
#include <fstream>

using namespace hsstab;
using namespace pf;

inline MKL_Complex16 getMKLCmplx(double ar, double ai){
	MKL_Complex16 r;
	r.real = ar; r.imag = ai;
	return r;}

inline MKL_Complex16 getMKLCmplx(t_CompVal val){
	MKL_Complex16 r;
	r.real = val.real();
	r.imag = val.imag();
	return r;
}

t_GlobSrchSpat::t_GlobSrchSpat(const mf::t_DomainBase& a_blk):_rBlk(a_blk),
_A(GS_NVARS_SPAT),
_B(GS_NVARS_SPAT),
_C(GS_NVARS_SPAT),
_A_AL(GS_NVARS_SPAT),
_B_AL(GS_NVARS_SPAT),
_C_AL(GS_NVARS_SPAT){
	_A.setToZero(); 
	_B.setToZero();
	_C.setToZero();
	_A_AL.setToZero();
	_B_AL.setToZero();
	_C_AL.setToZero();
};


void t_GlobSrchSpat::init(const hsstab::TPlugin& g_plug){

	const hsstab::TPluginParamsGroup& g = g_plug.get_settings_grp_const("");

	_params.init(g);

	const int& new_nnodes = _params.NNodes;

	_grid.resize(new_nnodes);

	_grid_y_stab.resize(new_nnodes);

	// container for insert vals
	// dim is max insert dim
	_insert_vals = new MKL_Complex16[GS_SO_BC_OUT_INS_SIZE];
	_insert_inds = new MKL_INT[GS_SO_BC_OUT_INS_SIZE];

	// global matrices & vectors

	_NDIM_G = _params.NNodes*(_params.NVars);

	int N_sq = _NDIM_G*_NDIM_G;

	_A_G = new MKL_Complex16[N_sq];
	_B_G = new MKL_Complex16[N_sq];

	// set matrices to zero
	for (int i=0; i<N_sq; i++){
		_A_G[i].real = 0.0;
		_A_G[i].imag = 0.0;

		_B_G[i].real = 0.0;
		_B_G[i].imag = 0.0;
	}

	_Alpha_G = new MKL_Complex16[_NDIM_G];
	_Beta_G  = new MKL_Complex16[_NDIM_G];

	// see MKL 11.1 user guide for "work" size reccomendations (lwork) 
	// zggev subroutine

	_LWork_G = 3*_NDIM_G;

	_Work_G = new MKL_Complex16[_LWork_G];
	_RWork_G = new double[3*_LWork_G];

	_spectrum.resize(_NDIM_G, 0.0);

}

t_GlobSrchSpat::~t_GlobSrchSpat(){

	delete[] _insert_vals;
	delete[] _insert_inds;

	delete[] _A_G, _B_G;

	delete[] _Alpha_G, _Beta_G, _Work_G, _RWork_G;

}

void t_GlobSrchSpat::setContext(const mf::t_GeomPoint& a_xyz){

	t_ProfMFLoc profNS(_rBlk);

	// TODO: control number of points in NS profile

	mf::t_ProfDataCfg prof_cfg;
	prof_cfg.ThickCoef = _params.ThickCoef;

	// NNodes can be used in INTERPOLATE types of initialization
	prof_cfg.NNodes = _params.NNodes;

	if (_params.NSProfInit == blp::NSINIT_SELFSIM_FROM_FILE) {
		prof_cfg.LoadFromAVFProfile = true;
	}

	prof_cfg.UseGlobalRFAsLocal = _rBlk.get_mf_params().BLUseGlobalRFAsLocal;

	switch (_params.NSProfInit)
	{
	case (blp::NSINIT_EXTRACT):
	case (blp::NSINIT_SELFSIM_FROM_FILE):
		profNS.initialize(a_xyz, prof_cfg, blp::NSINIT_EXTRACT);
		break;
	case (blp::NSINIT_INTERPOLATE):
		profNS.initialize(a_xyz, prof_cfg, blp::NSINIT_INTERPOLATE);
		break;
	default:
		wxString msg(_T("PF.GlobSearch: ProfNS Initialization type not supported"));
		wxLogError(msg); ssuGENTHROW(msg);
	}

	// TODO: play with coef
	// half nodes are placed into coef*bl_thick region
	const double Thick_HalfNodes_Coef = _params.ThickHalfNodesCoef;

	double y_max_mf = profNS.get_thick();
	double bl_thick_mf = y_max_mf/prof_cfg.ThickCoef;
	double y_i_mf = Thick_HalfNodes_Coef*bl_thick_mf;

	if (y_i_mf>=0.5*y_max_mf){
		wxString msg(_T("Spatial GS Error: computation domain too small, set larger Bl_Thick coefficient"));
		wxLogMessage(msg); ssuGENTHROW(msg);
	}

	double a_coef_mf = y_max_mf*y_i_mf/(y_max_mf - 2*y_i_mf);
	double b_coef_mf = 1.0 + a_coef_mf/y_max_mf;

	_deta = 1.0/(double)(_params.NNodes-1);	

	std::vector<double> grid_y_mf(_params.NNodes);

	for (int i=0; i<_params.NNodes; i++) {
		// computation grid
		_grid[i] = i*_deta;
		// physical grid
		grid_y_mf[i] = a_coef_mf*_grid[i]/(b_coef_mf - _grid[i]);
	}

	t_ProfStabCfg pstab_cfg;

	pstab_cfg.NNodes = _params.NNodes;
	pstab_cfg.NondimScaleType = _params.NondimScaleType;

	_profStab.initialize(profNS, grid_y_mf ,pstab_cfg);

	for (int i=0; i<_params.NNodes; i++) {
		_grid_y_stab[i] = _profStab.get_rec(i).y;
	}

	double y_max_stab = _profStab.get_thick();
	double bl_thick_stab = y_max_stab/prof_cfg.ThickCoef;
	double y_i_stab = Thick_HalfNodes_Coef*bl_thick_stab;

	//std::wcout<<"bl_thick_stab="<<y_max_stab/prof_cfg.ThickCoef;
	//_profStab.dump(_T("profStab_x0.2.dat"));

	_a_coef = y_max_stab*y_i_stab/(y_max_stab - 2*y_i_stab);
	_b_coef = 1.0 + _a_coef/y_max_stab;

	_set_curv_coefs_xyz(a_xyz);

}

void t_GlobSrchSpat::setContext(const t_ProfileStab* a_prof_stab){

		ssuGENTHROW(_T("GS Spatial: Check Init by profile_stab"));

}

/*void t_GlobSrchSpat::setContext(const t_ProfileStab* a_prof_stab, double bl_thick_scale){

	ssuGENTHROW(_T("GS Spatial: Check Init by profile_stab"));
	
	int nnodes = a_prof_stab->get_nnodes();

	_resize(nnodes);

	_profStab = *a_prof_stab;

	// TODO: configurable!
	const double Thick_HalfNodes_Coef = 2.0;

	double y_max = _profStab.get_thick();
	// Important: bl_thick_scale should be given in stab scales !
	double bl_thick = bl_thick_scale;
	double y_i = Thick_HalfNodes_Coef*bl_thick;

	if (y_i>=0.5*y_max){
		wxString msg(_T("Spatial GS Error: computation domain too small, set larger Bl_Thick coefficient"));
		wxLogMessage(msg); ssuGENTHROW(msg);
	}

	_a_coef = y_max*y_i/(y_max - 2*y_i);
	_b_coef = 1.0 + _a_coef/y_max;
	_deta = 1.0/(double)(_params.NNodes-1);	

	for (int i=0; i<nnodes; i++) {
		// computation grid
		_grid[i] = i*_deta;
		_grid_y_stab[i] = _profStab.get_rec(i).y;
	}
}*/

// semi-flag is true if it is k+1/2 point
void t_GlobSrchSpat::getMetricCoefs(const int a_nnode, double& f1, double& f2, double& f3, const bool semi_flag) const{

	double cur_eta = _grid[a_nnode];

	if (semi_flag) cur_eta+=0.5*(_grid[a_nnode+1]-_grid[a_nnode]);

	f3 = pow(_b_coef - cur_eta,2)/(_a_coef*_b_coef);

	f2 = -2.0*pow(_b_coef - cur_eta,3)/pow(_a_coef*_b_coef,2);

	f1 = f3*f3;
};

// semi-flag - is point like j+1/2 - first order equation
void t_GlobSrchSpat::setMatrices(const int a_nnode, const bool a_semi_flag){

	t_CompVal imagUnity(0.0, 1.0);

	const double& stabRe = _profStab.scales().ReStab;
	const double& Me = _profStab.scales().Me;
	const mf::t_FldParams& mf_params = _rBlk.get_mf_params();
	const double gMaMa = mf_params.Gamma*Me*Me;
	const double g_1MaMa = (mf_params.Gamma-1.0)*Me*Me;
	const double Pr = mf_params.Pr;

	double cur_y_stab = _grid_y_stab[a_nnode];
	if (a_semi_flag){
		cur_y_stab+=0.5*(_grid_y_stab[a_nnode+1]-_grid_y_stab[a_nnode]);
	};

	t_ProfRec rec = _profStab.get_rec(cur_y_stab);

	const double U = rec.u;
	const double U1 = rec.u1;
	const double U2 = rec.u2;

	const double W = rec.w;
	const double W1 = rec.w1;
	const double W2 = rec.w2;

	const double T = rec.t;
	const double T1 = rec.t1;
	const double T2 = rec.t2;

	const double mu = rec.mu;
	const double mu1 = rec.mu1;
	const double mu2 = rec.mu2;

	const double inv_t = 1.0/T;
	const double inv_mu = 1.0/mu;
	// k''/k:
	//const double mu_coef = mu1*t1*inv_mu;
	const double k2k = inv_mu*(mu2*pow(T1,2)+mu1*T2);

	const t_Complex dzeta_noA =  _beta*W - _w;
	const double dzeta_A = U;

// set _A
	_A.setToUnity();
	_A[2][2]=0.0;
// set _B
	const double sec_visc_ratio = _params.SecondViscRatio;

	const double L0 = 0.0 + sec_visc_ratio;
	const double L1 = 1.0 + sec_visc_ratio;
	const double L2 = 2.0 + sec_visc_ratio;

	const double inv_L2 = 1.0 / L2;

	const double m13 = _curv_coefs.m13;
	const double m31 = _curv_coefs.m31;
	const double m12 = _curv_coefs.m12;
	const double m32 = _curv_coefs.m32;

	// first row
	_B[0][0] = inv_mu*mu1*T1;
	_B[3][0] = inv_mu*mu1*U1;
	//second row
	_B[1][1] = _B[0][0];
	_B[2][1] = -stabRe/(L2*mu);
	_B[4][1] = imagUnity*_beta*L1/L2;
	// third
	_B[1][2] = 1.0;
	// fourth
	_B[0][3] = 2.0*g_1MaMa*Pr*U1;
	_B[3][3] = 2.0*inv_mu*mu1*T1;
	_B[4][3] = 2.0*g_1MaMa*Pr*W1;
	// last
	_B[1][4] = imagUnity*_beta*L1;
	_B[3][4] = inv_mu*mu1*W1;
	_B[4][4] = inv_mu*mu1*T1;

// set _C
	// first row, u-eq
	_C[0][0] = -imagUnity*stabRe*dzeta_noA*inv_mu*inv_t-_beta*_beta;

	_C[1][0] = -stabRe*U1*inv_mu*inv_t; 

	_C[2][0] = 0.0;

	_C[3][0] = inv_mu*(mu1*U2 + mu2*T1*U1);

	_C[4][0] = 0.0;

	// second row, v-eq
	_C[0][1] = 0.0;

	_C[1][1] = -imagUnity*dzeta_noA*stabRe/(L2*mu*T)-pow(_beta,2)/L2;

	_C[2][1] = 0.0;

	_C[3][1] = imagUnity*inv_mu*mu1*_beta*W1/L2;

	_C[4][1] = imagUnity*_beta*inv_mu*mu1*T1*L0/L2;

	// third row, p-eq

	_C[0][2] = 0.0;

	_C[1][2] = -T1*inv_t;

	_C[2][2] = imagUnity*gMaMa*dzeta_noA;

	_C[3][2] = -imagUnity*dzeta_noA*inv_t;

	_C[4][2] = imagUnity*_beta;

	//fourth row, T-eq
	_C[1][3] = Pr*(2.0*imagUnity*g_1MaMa*_beta*W1 - stabRe*T1*inv_mu*inv_t);

	_C[2][3] = imagUnity*dzeta_noA*g_1MaMa*Pr*stabRe*inv_mu;

	_C[3][3] = -imagUnity*dzeta_noA*stabRe*Pr*inv_mu*inv_t - pow(_beta,2)
			   +g_1MaMa*Pr*inv_mu*mu1*(pow(U1,2)+pow(W1,2)) + k2k;

	// fifth row, w-eq

	_C[0][4] = 0.0;

	_C[1][4] = imagUnity*_beta*inv_mu*mu1*T1 - stabRe*W1*inv_mu*inv_t;

	_C[2][4] = -imagUnity*_beta*stabRe*inv_mu;

	_C[3][4] = inv_mu*mu1*W2 + inv_mu*mu2*T1*W1;

	_C[4][4] = -imagUnity*dzeta_noA*stabRe*inv_mu*inv_t - L2*pow(_beta,2);

	const bool curv_on = _params.CurvTermsOn;

	// compressible additions
	/*
	if (curv_on) {

		// full additions after sousa's eqs
		// first row
		_C[0][0] += -stabRe*inv_mu*inv_t*W*m12;

		_C[1][0] += -stabRe*inv_mu*inv_t*U*m12;

		_C[2][0] += -stabRe*inv_mu*inv_t*gMaMa*(m12*U*W - m31*W*W);

		_C[3][0] += stabRe*inv_mu*inv_t*inv_t*(m12*U*W - m31*W*W);

		_C[4][0] += stabRe*inv_mu*inv_t*(2.0*W*m31 - U*m12);

		// second row

		_C[0][1] += -2.0*stabRe*inv_mu*inv_t*inv_L2*U*m12;

		_C[2][1] += stabRe*gMaMa*inv_mu*inv_t*inv_L2*(m12*U*U + m32*W*W);

		_C[3][1] += -stabRe*inv_mu*inv_L2*inv_t*inv_t*(m12*U*U + m32*W*W);

		_C[4][1] += 2.0*stabRe*inv_mu*inv_L2*inv_t*W*m32;

		// thrid row

		_C[0][2] += m31;

		_C[1][2] += m12 + m32;

		_C[2][2] += gMaMa*(m31*U + m12*W);

		_C[3][2] += -inv_t*(m31*U + m12*W);

		_C[4][2] += m12;

		// fourth row

		// fifth row
		_C[0][4] += stabRe*inv_mu*inv_t*(2.0*U*m12 - W*m31);

		_C[1][4] += -stabRe*inv_mu*inv_t*W*m32;

		_C[2][4] += -stabRe*inv_mu*inv_t*gMaMa*(m31*W*U - m12*U*U);

		_C[3][4] += stabRe*inv_mu*inv_t*inv_t*(m31*W*U - m12*U*U);

		_C[4][4] += -stabRe*inv_mu*inv_t*U*m31;

	} // compressible additions 
	*/

	// incompressible additions after Cebeci
	if (curv_on) {

		//wxLogMessage(_T("GS: Curv terms on!"));

		// first row
		_C[0][0] += -stabRe*W*m13;

		_C[1][0] += -stabRe*U*m12;

		_C[4][0] += stabRe*(2.0*W*m31 - U*m13);

		// second row

		_C[0][1] += 2.0*stabRe*U*m12;

		_C[4][1] += 2.0*stabRe*W*m32;

		// thrid row

		_C[0][2] += m31;

		_C[1][2] += m12 + m32;

		_C[4][2] += m13;

		// fourth row

		// fifth row
		_C[0][4] += stabRe*(2.0*U*m13 - W*m31);

		_C[1][4] += -stabRe*W*m32;

		_C[4][4] += -stabRe*U*m31;

	}
// Alpha matrices
	// _A_AL is always zero

	// set _B_AL
	// first row
	_B_AL[1][0] = imagUnity*L1;
	// second row
	_B_AL[0][1] = imagUnity*L1/L2;

	// set _C_AL
	// first row
	_C_AL[0][0] = -imagUnity*dzeta_A*stabRe*inv_mu*inv_t;

	_C_AL[1][0] = imagUnity*inv_mu*mu1*T1;

	_C_AL[2][0] = -imagUnity*stabRe*inv_mu;

	_C_AL[4][0] = -_beta*L1;

	//second row
	_C_AL[0][1] = imagUnity*inv_mu*mu1*T1*L0/L2;

	_C_AL[1][1] = -imagUnity*dzeta_A*stabRe/L2*inv_mu*inv_t;

	_C_AL[3][1] = imagUnity*inv_mu*mu1*U1/L2;

	// third row
	_C_AL[0][2] = imagUnity;

	_C_AL[2][2] = imagUnity*gMaMa*dzeta_A;

	_C_AL[3][2] = -imagUnity*dzeta_A*inv_t;

	// fourth row
	_C_AL[1][3] = 2.0*imagUnity*g_1MaMa*Pr*U1;

	_C_AL[2][3] = imagUnity*dzeta_A*g_1MaMa*Pr*stabRe*inv_mu;

	_C_AL[3][3] = -imagUnity*dzeta_A*stabRe*Pr*inv_mu*inv_t;

	// fifth row
	_C_AL[0][4] = -_beta*L1;

	_C_AL[4][4] = -imagUnity*dzeta_A*stabRe*inv_mu*inv_t;

};

// tmp, set curvature coefs by hand
// current values for sphere

void t_GlobSrchSpat::_set_curv_coefs_xyz(const mf::t_GeomPoint& a_xyz) {

	// sphere case, dimensional raduis
	const double R_dim = 0.089;

	_curv_coefs.set_coefs_sphere(_profStab.scales(), R_dim);


}
// a_eq_id={0,1,3,4} <--> {1,2,4,5} - SO
void t_GlobSrchSpat::fill_SO_row(const t_SqMatCmplx& a_LMat, const t_SqMatCmplx& a_MMat, 
								 const t_SqMatCmplx& a_RMat, const int a_nnode, const int a_eq_id, int& a_ins_nvals){

		
	// large and small index bases
	int l_base = (a_nnode-1)*GS_NVARS_SPAT;
	int s_base=0;
// uniform grid
	const double step = _deta;
	const double inv_step = 1.0/step;
	const double inv_step2 = inv_step*inv_step;
	double f1, f2, f3;
	getMetricCoefs(a_nnode, f1, f2, f3, false);

	a_ins_nvals = GS_SO_INS_SIZE;

	// first five elements in inserts: f~[j-1] and f^[j-1/2]
	// f~<=>{u,v,o,w,t} ; f^<=>{0,0,p,0,0}
	for (int i=0; i<GS_NVARS_SPAT; i++){
		if (i==2){
			_insert_vals[s_base+i] = getMKLCmplx(-f3*inv_step*a_MMat[i][a_eq_id]
									 +0.5*a_RMat[i][a_eq_id]);
		}else{
			_insert_vals[s_base+i] = getMKLCmplx(
				(f1*inv_step2-0.5*f2*inv_step)*a_LMat[i][a_eq_id]
				-0.5*inv_step*f3*a_MMat[i][a_eq_id]);
		};
		_insert_inds[s_base+i] = l_base + i;
	};
	l_base+=_params.NVars;
	s_base+=_params.NVars;
	
	// fill central 5 points:  f~[j] and f^[j+1/2]
	for (int i=0; i<_params.NVars; i++){
		if (i==2){
			_insert_vals[s_base+i] = getMKLCmplx(f3*inv_step*a_MMat[i][a_eq_id]
									 +0.5*a_RMat[i][a_eq_id]);
		}else{
			_insert_vals[s_base+i] = getMKLCmplx(-2.0*f1*inv_step2*a_LMat[i][a_eq_id]
									 +a_RMat[i][a_eq_id]);
		};	
		_insert_inds[s_base+i] = l_base + i;
	};
	l_base+=_params.NVars;
	s_base+=_params.NVars;
	// fill last five points : f~[j+1] + f^[j+3/2]
	for (int i=0; i<_params.NVars; i++){
		if (i==2){
			_insert_vals[s_base+i] = getMKLCmplx(0.0,0.0);
		}else{
			_insert_vals[s_base+i] = getMKLCmplx(
				(f1*inv_step2+0.5*f2*inv_step)*a_LMat[i][a_eq_id]
				+0.5*inv_step*f3*a_MMat[i][a_eq_id]);
		};	
		_insert_inds[s_base+i] = l_base + i;
	};
};

t_CompVal t_GlobSrchSpat::_calc_p_out_eigval() {

	// NSOL_VECS_ODES and STAB_MATRIX_DIM
	t_MatCmplx vecs(4,8);
	t_CompVal lambdas[4];

	t_WCharsLoc wave;
	wave.b = this->_beta;
	wave.w = this->_w;
	// tricky part 
	// a = 0.0 worked for M=3 parab wing test
	wave.a = 0.0;//t_CompVal(0.023531, -0.001602);// wave.w;
	_p_loc_solver->setWave(wave);

	_p_loc_solver->setAsymptotics(vecs, &lambdas[0]);

	const t_StabScales& ls_scales = _p_loc_solver->get_stab_scales();

	const t_StabScales& gs_scales = this->_profStab.scales();

	wxLogMessage(_T("Debug: ls Dels=%lf, gs Dels=%lf"), ls_scales.Dels, gs_scales.Dels);

	// 
	t_CompVal p_eig_val = lambdas[2];

	t_CompVal ret = p_eig_val * gs_scales.Dels / ls_scales.Dels;
	wxLogMessage(_T("Debug: computed p eigval: (%lf, %lf)"), ret.real(), ret.imag());

	return ret;
}

void t_GlobSrchSpat::fill_SO_row_BC_out(const t_SqMatCmplx& a_LMat, const t_SqMatCmplx& a_MMat,
	const t_SqMatCmplx& a_RMat, const int a_nnode, const int a_eq_id, int& a_ins_nvals) {

	// large and small index bases
	int l_base = (a_nnode - 3)*GS_NVARS_SPAT;
	int s_base = 0;
	// uniform grid
	const double step = _deta;
	const double inv_step = 1.0 / step;
	const double inv_step2 = inv_step*inv_step;
	double f1, f2, f3;
	getMetricCoefs(a_nnode, f1, f2, f3, false);

	a_ins_nvals = GS_SO_BC_OUT_INS_SIZE;

	// first five elements in inserts: f~[j-3] and f^[j-5/2]
	// f~<=>{u,v,o,w,t} ; f^<=>{0,0,p,0,0}
	for (int i = 0; i<GS_NVARS_SPAT; i++) {
		if (i == 2) {
			_insert_vals[s_base + i] = getMKLCmplx(0.0);
		}
		else {
			_insert_vals[s_base + i] = getMKLCmplx(
				(-f1*inv_step2/3.0)*a_LMat[i][a_eq_id]);
		};
		_insert_inds[s_base + i] = l_base + i;
	};
	l_base += _params.NVars;
	s_base += _params.NVars;

	// next five elements in inserts: f~[j-2] and f^[j-3/2]
	// f~<=>{u,v,o,w,t} ; f^<=>{0,0,p,0,0}
	for (int i = 0; i<GS_NVARS_SPAT; i++) {
		if (i == 2) {
			_insert_vals[s_base + i] = getMKLCmplx(0.0);
		}
		else {
			_insert_vals[s_base + i] = getMKLCmplx(
				(2.0*f1*inv_step2+0.5*f2*inv_step)*a_LMat[i][a_eq_id]
			+0.5*f3*inv_step*a_MMat[i][a_eq_id]);
		};
		_insert_inds[s_base + i] = l_base + i;
	};
	l_base += _params.NVars;
	s_base += _params.NVars;

	// next five elements in inserts: f~[j-1] and f^[j-1/2]
	// f~<=>{u,v,o,w,t} ; f^<=>{0,0,p,0,0}

	t_CompVal lambda = _calc_p_out_eigval();
	t_CompVal eps = (2.0*f3 + lambda*step) / (2.0*f3 - lambda*step);

	for (int i = 0; i<GS_NVARS_SPAT; i++) {
		if (i == 2) {
			_insert_vals[s_base + i] = getMKLCmplx(
				(eps-1.0)*f3*inv_step*a_MMat[i][a_eq_id] + 0.5*(eps+1.0)*a_RMat[i][a_eq_id]);
		}
		else {
			_insert_vals[s_base + i] = getMKLCmplx(
				(-3.0*f1*inv_step2 - 2.0*f2*inv_step)*a_LMat[i][a_eq_id]
				- 2.0*f3*inv_step*a_MMat[i][a_eq_id]);
		};
		_insert_inds[s_base + i] = l_base + i;
	};
	l_base += _params.NVars;
	s_base += _params.NVars;

	// last five elements in inserts: f~[j] and f^[j+1/2]
	// f~<=>{u,v,o,w,t} ; f^<=>{0,0,p,0,0}
	for (int i = 0; i<GS_NVARS_SPAT; i++) {
		if (i == 2) {
			_insert_vals[s_base + i] = getMKLCmplx(0.0);
		}
		else {
			_insert_vals[s_base + i] = getMKLCmplx(
				(4.0/3.0*f1*inv_step2 + 1.5*f2*inv_step)*a_LMat[i][a_eq_id]
				+ 1.5*f3*inv_step*a_MMat[i][a_eq_id] + a_RMat[i][a_eq_id]);
		};
		_insert_inds[s_base + i] = l_base + i;
	};
}

// the only first order continuity equation 
void t_GlobSrchSpat::fill_FO_row(const t_SqMatCmplx& a_MMat, const t_SqMatCmplx& a_RMat, 
								 const int a_nnode, int& a_ins_nvals){

	const int a_eq_id = 2;
	// large and small index bases
	int l_base = a_nnode*_params.NVars;
	int s_base=0;
// uniform grid
	const double step = _deta;
	const double inv_step = 1.0/step;
	const double inv_step2 = inv_step*inv_step;
	double f1, f2, f3;
	getMetricCoefs(a_nnode, f1, f2, f3, true);

	a_ins_nvals = GS_FO_INS_SIZE;

	// first block f~[j] and f^[j+1/2]:
	for (int i=0; i<_params.NVars; i++){
		if (i==2){
			_insert_vals[s_base+i] = getMKLCmplx(a_RMat[i][a_eq_id]);
		}else{
			_insert_vals[s_base+i] = getMKLCmplx(-f3*inv_step*a_MMat[i][a_eq_id]
									 +0.5*a_RMat[i][a_eq_id]);
		};
		_insert_inds[s_base+i] = l_base + i;
	};
	l_base+=_params.NVars;
	s_base+=_params.NVars;

	// second block f~[j+1] + f^[j+3/2]:
	for (int i=0; i<_params.NVars; i++){
		if (i==2){
			_insert_vals[s_base+i] = getMKLCmplx(0.0, 0.0);
		}else{
			_insert_vals[s_base+i] = getMKLCmplx(f3*inv_step*a_MMat[i][a_eq_id]
										+0.5*a_RMat[i][a_eq_id]);
		};
		_insert_inds[s_base+i] = l_base + i;
	};
};

inline int calcPlainIndByIJ(int row, int col, int ndim){
	return col+row*ndim;
}

inline void insertValsToMat(MKL_Complex16* mat, int mat_size, int row_num, MKL_Complex16* vals, MKL_INT* inds, int n_inserts){

	int g_ind;

	for (int j=0; j<n_inserts; j++){

		g_ind = calcPlainIndByIJ(row_num, inds[j], mat_size);
		mat[g_ind] = vals[j];

	}

}


int t_GlobSrchSpat::_solve(){

	static char help[]="Global Search\n";

	wxLogMessage(_("\nGlobal Eigensearch started: N=%d\n"), _params.NNodes);

	int n_inserts = 0;

	// IMPORTANT: reset global matrices -
	// they are corrupted after zggev call...
	for (int i=0; i<_NDIM_G*_NDIM_G; i++){

		_A_G[i] = getMKLCmplx(0.0,0.0);
		_B_G[i] = getMKLCmplx(0.0,0.0);
	}

	for (int i=0; i<_LWork_G; i++){
		_Work_G[i] = getMKLCmplx(0.0,0.0);;
	}

	for (int i=0; i<3*_LWork_G; i++){
		_RWork_G[i] = 0.0;
	}

	for (int i=0; i<_NDIM_G; i++){
		_Alpha_G[i] = getMKLCmplx(0.0,0.0);;
		_Beta_G[i] = getMKLCmplx(0.0,0.0);;
	}

	// fill A~ and B~ by rows
	// first block - bc,bc,fo,bc,bc
	int row_num=0;
	{
		for (int j=0; j<GS_NVARS_SPAT; j++){

			row_num=j;

			if (j!=2){

				// j=0,1,3,4<=> u~,v~,t~,w~ - apply homogeneos bc
				// insert_val is an arbitrary value : phi=0 <=> phi + alpha*eps*phi = 0
				// TODO: check that the result doesn't depend on eps
				// val_b=eps=1=> (1+alpha)*{u,v,w,T}=0 => {u,v,w,T}=0
				// TODO: try val_b=eps=-1
				// TODO: make different BC options, e.g. dt/dy=0 !

				int g_ind = calcPlainIndByIJ(row_num, row_num, _NDIM_G);

				MKL_Complex16 one; one.real = 1.0; one.imag = 0.0;

				_A_G[g_ind] = one;
				_B_G[g_ind] = one;

			}else{

				// ordinary FO row for p^[0+1/2]
				setMatrices(0, true);
				fill_FO_row(_B, _C, 0, n_inserts);

				insertValsToMat(_A_G, _NDIM_G, row_num, _insert_vals, _insert_inds, n_inserts);				

				fill_FO_row(_B_AL, _C_AL, 0, n_inserts);

				insertValsToMat(_B_G, _NDIM_G, row_num, _insert_vals, _insert_inds, n_inserts);

			}
		}
	}
	// main body of the matrices A, B - {SO, SO, FO, SO, SO} blocks
	for (int i=1; i<_params.NNodes-1; i++){
		for (int j=0; j<GS_NVARS_SPAT; j++){

			row_num = GS_NVARS_SPAT*i + j;

			if (j==2){
				setMatrices(i, true);

				fill_FO_row(_B, _C, i, n_inserts);
				insertValsToMat(_A_G, _NDIM_G, row_num, _insert_vals, _insert_inds, n_inserts);

				fill_FO_row(_B_AL, _C_AL, i, n_inserts);
				insertValsToMat(_B_G, _NDIM_G, row_num, _insert_vals, _insert_inds, n_inserts);

			}else{
				setMatrices(i, false);

				fill_SO_row(_A, _B, _C,i, j, n_inserts);
				insertValsToMat(_A_G, _NDIM_G, row_num, _insert_vals, _insert_inds, n_inserts);

				fill_SO_row(_A_AL, _B_AL, _C_AL,i, j, n_inserts);
				insertValsToMat(_B_G, _NDIM_G, row_num, _insert_vals, _insert_inds, n_inserts);
			};
		}
	};
	// last block 5x5
	// just set blocks to unity matrices
	{
		if (_params.BCOutKind==t_EigenGSParams::t_BCOutKind::BC_OUT_HOMOGEN){

			wxLogMessage(_T("homogen"));

			for (int j = 0; j<GS_NVARS_SPAT; j++) {

				row_num = (_params.NNodes - 1)*GS_NVARS_SPAT + j;

				int g_ind = calcPlainIndByIJ(row_num, row_num, _NDIM_G);

				MKL_Complex16 one; one.real = 1.0; one.imag = 0.0;

				_A_G[g_ind] = one;
				_B_G[g_ind] = one;

			}
		}
		
		if (_params.BCOutKind == t_EigenGSParams::t_BCOutKind::BC_OUT_P_ASYM) {
			wxLogMessage(_T("p_asym"));

			int i = _params.NNodes - 1;

			for (int j = 0; j<GS_NVARS_SPAT; j++) {

				row_num = GS_NVARS_SPAT*i + j;

				if (j == 2) {

					int g_ind = calcPlainIndByIJ(row_num, row_num, _NDIM_G);

					MKL_Complex16 one; one.real = 1.0; one.imag = 0.0;

					_A_G[g_ind] = one;
					_B_G[g_ind] = one;

				}else {
					setMatrices(i, false);

					fill_SO_row_BC_out(_A, _B, _C, i, j, n_inserts);
					insertValsToMat(_A_G, _NDIM_G, row_num, _insert_vals, _insert_inds, n_inserts);

					fill_SO_row_BC_out(_A_AL, _B_AL, _C_AL, i, j, n_inserts);
					insertValsToMat(_B_G, _NDIM_G, row_num, _insert_vals, _insert_inds, n_inserts);
				};
			}

		}

	}

	char jobvl = 'N'; 
	char jobvr = 'N';

	MKL_Complex16 vl[1]; 
	MKL_Complex16 vr[1];

	MKL_INT ldvl = 1; 
	MKL_INT ldvr = 1;


	MKL_INT info;

	zggev( &jobvl, &jobvr, &_NDIM_G, _A_G, &_NDIM_G, _B_G, &_NDIM_G,
		_Alpha_G, _Beta_G, vl, &ldvl, vr, &ldvr, _Work_G, &_LWork_G, _RWork_G, &info );



	for (int i=0; i<_NDIM_G; i++){

		t_CompVal a_s(_Alpha_G[i].real, _Alpha_G[i].imag);
		t_CompVal b_s(_Beta_G[i].real, _Beta_G[i].imag);

		// searching solution of Ah + wBh=0
		// solving Ah = lBh => w = -l !!!
		_spectrum[i] = -a_s/b_s;

	}

	return 0;
}

int t_GlobSrchSpat::getSpectrum(const t_WCharsLoc& init_wave){

  _beta  = init_wave.b;
  _w = init_wave.w;

  int err_code = _solve();
  return err_code;

};

void t_GlobSrchSpat::writeSpectrum(const std::string &a_filename) const{
	std::wofstream to_f(&a_filename[0]);
	for (int i=0; i<_spectrum.size(); i++){
		to_f<<_spectrum[i].real()<<_T("\t")<<_spectrum[i].imag()<<std::endl;
	};
}

void t_GlobSrchSpat::writeSpectrumPhase(const std::string &a_filename) const{

	wxLogMessage(_T("writeSpectrumPhase Warning: ci,cr comps are approximate!\n"));
	
	//double k = sqrt(_alpha*_alpha + _beta*_beta);
	t_Complex c;
	std::wofstream to_f(&a_filename[0]);
	for (int i=0; i<_spectrum.size(); i++){

		c = _w/_spectrum[i];
		
		to_f<<c.real()<<_T("\t")<<c.imag()<<_T("\n");
	};
	
}

std::vector<t_WCharsLoc> t_GlobSrchSpat::getInstabModes(const t_WCharsLoc& init_wave){

		std::vector<t_WCharsLoc> inits;
		std::vector<t_Complex>::const_iterator it;
		getSpectrum(init_wave);
		// TODO: empirics!!!
		for (it=_spectrum.begin(); it<_spectrum.end(); it++){
			// IMPORTANT TODO: fix parasitic solutions, ask AVF
			// now neglect parasitic branch manually
			if (it->imag()<_params.Arg_Threshold && it->real()>0.0){
				t_WCharsLoc wave;
				wave.a = *it;
				wave.b = _beta;
				wave.w = _w;
				wave.set_treat(stab::t_TaskTreat::SPAT);
				wave.set_scales(get_stab_scales());
				inits.push_back(wave);
			}
		}
		return inits;
};

t_WCharsLoc t_GlobSrchSpat::searchMaxInstab(const t_WCharsLoc& init_wave){
	const std::vector<t_WCharsLoc>& all_initials = getInstabModes(init_wave);
	return t_WCharsLoc::find_max_instab_spat(all_initials);
};

const t_StabScales& t_GlobSrchSpat::get_stab_scales() const{
	return _profStab.scales();
};