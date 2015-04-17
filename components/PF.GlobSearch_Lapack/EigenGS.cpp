#include "stdafx.h"

#include "common_data.h"

#include "EigenGS.h"

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

t_EigenGS::t_EigenGS(const mf::t_DomainBase& a_blk):_rBlk(a_blk),
_A(GS_NVARS_TIME),
_B(GS_NVARS_TIME),
_C(GS_NVARS_TIME),
_CW(GS_NVARS_TIME){
	_A.setToZero(); 
	_B.setToZero();
	_C.setToZero();
	_CW.setToZero();
};

void t_EigenGS::init(const hsstab::TPlugin& g_plug){

	const hsstab::TPluginParamsGroup& g = g_plug.get_settings_grp_const("");

	gs::_init_eigen_gs_base_params(_params, g);

	const int& new_nnodes = _params.NNodes;

	_grid.resize(new_nnodes);

	_grid_y_stab.resize(new_nnodes);

	// container for insert vals
	// dim is max insert dim
	_insert_vals = new MKL_Complex16[3*_params.NVars];
	_insert_inds = new MKL_INT[3*_params.NVars];


	// global matrices & vectors

	_NDIM_G = (_params.NNodes-1)*(_params.NVars);

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

t_EigenGS::~t_EigenGS(){

	delete[] _insert_vals;
	delete[] _insert_inds;

	delete[] _A_G, _B_G;

	delete[] _Alpha_G, _Beta_G, _Work_G, _RWork_G;

}

void t_EigenGS::setContext(const mf::t_GeomPoint& a_xyz){

	t_ProfileNS profNS(_rBlk);

	// TODO: control number of points in NS profile

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

	_profStab.initialize(profNS, grid_y_mf ,_params.NNodes);

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

}

// semi-flag is true if it is k-1/2 point
void t_EigenGS::getMetricCoefs(const int a_nnode, double& f1, double& f2, double& f3, const bool semi_flag) const{
	double cur_eta = _grid[a_nnode];
	if (semi_flag){
		cur_eta-=0.5*(_grid[a_nnode]-_grid[a_nnode-1]);
	};
	f3 = pow(_b_coef - cur_eta,2)/(_a_coef*_b_coef);
	f2 = -2.0*pow(_b_coef - cur_eta,3)/
		     pow(_a_coef*_b_coef,2);
	f1 = f3*f3;
};

// semi-flag - is point like j+1/2 - first order equation
void t_EigenGS::setMatrices(const int a_nnode, const bool a_semi_flag){
	t_CompVal imagUnity(0.0, 1.0);

	const double& stabRe = _profStab.scales().ReStab;
	const double& Me = _profStab.scales().Me;
	const mf::t_FldParams& mf_params = _rBlk.get_mf_params();
	const double gMaMa = mf_params.Gamma*Me*Me;
	const double g_1MaMa = (mf_params.Gamma-1.0)*Me*Me;
	const double Pr = mf_params.Pr;

	// rename a_y is not input but lazy to rewrite
	// get physical y

	double cur_y_stab = _grid_y_stab[a_nnode];
	if (a_semi_flag){
		cur_y_stab-=0.5*(_grid_y_stab[a_nnode]-_grid_y_stab[a_nnode-1]);
	};

	const double a_y = cur_y_stab;
	t_ProfRec rec = _profStab.get_rec(a_y);

	const double inv_t = 1.0/rec.t;
	const double inv_mu = 1.0/rec.mu;
	// k''/k:
	//const double mu_coef = mu1*t1*inv_mu;
	const double k2k = inv_mu*(rec.mu2*pow(rec.t1,2)+rec.mu1*rec.t2);	//2.0*pow(mu_coef,2)-mu_coef;

	const double dzeta_noW = _alpha*rec.u + _beta*rec.w;
	const double dzeta_W = 1.0;

// set _A
	_A.setToUnity();
	_A[2][2]=0.0;
// set _B
	const double lz = 0.0; // 0 + lambda/mu
	const double lf = 1.0; // 1 + lambda/mu
	const double ls = 2.0; // 2 + lambda/mu

	// first row
	_B[0][0] = inv_mu*rec.mu1*rec.t1;
	_B[1][0] = imagUnity*_alpha*lf;
	_B[3][0] = inv_mu*rec.mu1*rec.u1;
	//second row
	_B[0][1] = imagUnity*_alpha*lf/ls;
	_B[1][1] = _B[0][0];
	_B[2][1] = -stabRe/(ls*rec.mu);
	_B[4][1] = imagUnity*_beta*lf/ls;
	// third
	_B[1][2] = 1.0;
	// fourth
	_B[0][3] = 2.0*g_1MaMa*Pr*rec.u1;
	// k'/k ??
	_B[3][3] = 2.0*rec.mu1*rec.t1*inv_mu; //TODO: was -2.0*... - bug ?
	_B[4][3] = 2.0*g_1MaMa*Pr*rec.w1;
	// last
	_B[1][4] = imagUnity*_beta*lf;
	_B[3][4] = inv_mu*rec.mu1*rec.w1;
	_B[4][4] = inv_mu*rec.mu1*rec.t1;
// clear C - no freq terms
	// first row
	_C[0][0] = -imagUnity*dzeta_noW*stabRe*inv_mu*inv_t
		       -(ls*pow(_alpha,2)+pow(_beta,2));

	_C[1][0] = -stabRe*rec.u1*inv_mu*inv_t
		       +imagUnity*_alpha*inv_mu*rec.mu1*rec.t1;

	_C[2][0] = -imagUnity*_alpha*stabRe*inv_mu;

	_C[3][0] = inv_mu*(rec.mu1*rec.u2 + rec.mu2*rec.t1*rec.u1);
	_C[4][0] = -_alpha*_beta*lf;
	// second

	_C[0][1] = imagUnity*_alpha*inv_mu*rec.mu1*rec.t1*lz/ls;

	_C[1][1] = -imagUnity*dzeta_noW*stabRe/(ls*rec.mu*rec.t)
		       -(pow(_alpha,2)+pow(_beta,2))/ls;

	_C[3][1] = imagUnity*inv_mu*rec.mu1*(_alpha*rec.u1+_beta*rec.w1)/ls;
	_C[4][1] = imagUnity*_beta*inv_mu*rec.mu1*rec.t1*lz/ls;
	// third
	_C[0][2] = imagUnity*_alpha;
	_C[1][2] = -rec.t1*inv_t;
	_C[2][2] = imagUnity*gMaMa*dzeta_noW;
	_C[3][2] = -imagUnity*dzeta_noW*inv_t;
	_C[4][2] = imagUnity*_beta;
	//fourh
	_C[1][3] = Pr*
			   (
			    2.0*imagUnity*g_1MaMa*(_alpha*rec.u1 + _beta*rec.w1)
			    -stabRe*rec.t1*inv_mu*inv_t);

	_C[2][3] = imagUnity*dzeta_noW*g_1MaMa*Pr
		       *stabRe*inv_mu;

	_C[3][3] = -imagUnity*dzeta_noW*stabRe*Pr*inv_mu*inv_t
		       -(pow(_alpha,2)+pow(_beta,2))
			   +g_1MaMa*Pr*inv_mu*rec.mu1*(pow(rec.u1,2)+pow(rec.w1,2))
			   +k2k;
	// last at least
	_C[0][4] = -_alpha*_beta*lf;
	_C[1][4] = imagUnity*_beta*inv_mu*rec.mu1*rec.t1
		       -stabRe*rec.w1*inv_mu*inv_t;

	_C[2][4] = -imagUnity*_beta*stabRe*inv_mu;
	_C[3][4] = inv_mu*rec.mu1*rec.w2 + inv_mu*rec.mu2*rec.t1*rec.w1;

	_C[4][4] = -imagUnity*dzeta_noW*stabRe*inv_mu*inv_t
		       -(pow(_alpha,2)+ls*pow(_beta,2));
	// W matrix
	_CW[0][0] = -imagUnity*dzeta_W*stabRe*inv_mu*inv_t;
	_CW[1][1] = -imagUnity*dzeta_W*stabRe/ls*inv_mu*inv_t;
	_CW[2][2] = imagUnity*gMaMa*dzeta_W;
	_CW[3][2] = -imagUnity*dzeta_W*inv_t;
	_CW[2][3] = imagUnity*dzeta_W*g_1MaMa*Pr*stabRe*inv_mu;
	_CW[3][3] = -imagUnity*dzeta_W*stabRe*Pr*inv_mu*inv_t;
	_CW[4][4] = -imagUnity*dzeta_W*stabRe*inv_mu*inv_t;

};
// a_eq_id={0,1,3,4} <--> {1,2,4,5} - SO
void t_EigenGS::fill_SO_template(const t_SqMatCmplx& a_LMat, const t_SqMatCmplx& a_MMat, const t_SqMatCmplx& a_RMat, 
								 const int a_nnode, const int a_eq_id, int& a_ins_nvals){
	if (a_nnode==0){
		// TODO: to log
		std::wcerr<<_T("GS: SO template on bottom boundary!\n");
	};
	bool first_block=false;
	if (a_nnode==1){
		first_block=true;
		// small template 2 point
		a_ins_nvals = 2*_params.NVars;
	}else{
		if (a_nnode==_params.NNodes-1){
			// "fake" template
			// for diag matrix
			a_ins_nvals = 1;
			_insert_vals[0] = getMKLCmplx(1.0, 0.0);
			_insert_inds[0] = (a_nnode-1)*_params.NVars+a_eq_id;
			return;
		}else{
			// full 3-point template
			a_ins_nvals = 3*_params.NVars;
		};
		
	};
	// large and small index bases
	int l_base = first_block ? 0 : (a_nnode-2)*_params.NVars;
	int s_base=0;
// uniform grid
	const double step = _deta;
	const double inv_step = 1.0/step;
	const double inv_step2 = inv_step*inv_step;
	double f1, f2, f3;
	getMetricCoefs(a_nnode, f1, f2, f3, false);
	if (!first_block){
	// first five elements in inserts
	// f~[j-1] : u,v,0,w,t
		for (int i=0; i<_params.NVars; i++){
			if (i==2){
				_insert_vals[s_base+i] = getMKLCmplx(0.0, 0.0);
			// f~(k-1):
			}else{
				_insert_vals[s_base+i] = getMKLCmplx(
					(f1*inv_step2-0.5*f2*inv_step)*a_LMat[i][a_eq_id]
					-0.5*inv_step*f3*a_MMat[i][a_eq_id]);
			};
			_insert_inds[s_base+i] = l_base + i;
		};
		l_base+=_params.NVars;
		s_base+=_params.NVars;
	};
	
	// fill central 5 points: 
	// f~[j] + f^[j-1/2]
	for (int i=0; i<_params.NVars; i++){
		if (i==2){
			_insert_vals[s_base+i] = getMKLCmplx(
				-f3*inv_step*a_MMat[i][a_eq_id]
				+0.5*a_RMat[i][a_eq_id]);
		}else{
			_insert_vals[s_base+i] = getMKLCmplx(
				-2.0*f1*inv_step2*a_LMat[i][a_eq_id]
				+a_RMat[i][a_eq_id]);
		};	
		_insert_inds[s_base+i] = l_base + i;
	};
	l_base+=_params.NVars;
	s_base+=_params.NVars;
	// fill last five points
	// f~[j+1] + f^[j+1/2]
	for (int i=0; i<_params.NVars; i++){
		if (i==2){
			_insert_vals[s_base+i] = getMKLCmplx(
				f3*inv_step*a_MMat[i][a_eq_id]
				+0.5*a_RMat[i][a_eq_id]);
		}else{
			_insert_vals[s_base+i] = getMKLCmplx(
				(f1*inv_step2+0.5*f2*inv_step)*a_LMat[i][a_eq_id]
				+0.5*inv_step*f3*a_MMat[i][a_eq_id]);
		};	
		_insert_inds[s_base+i] = l_base + i;
	};
};
// we have only one first order continuity equation 
// a_eq_id = 2 
// this template is in staggered point a_nnode - 1/2
void t_EigenGS::fill_FO_template(const t_SqMatCmplx& a_MMat, const t_SqMatCmplx& a_RMat, 
								 const int a_nnode, const int a_eq_id, int &a_ins_nvals){
	bool first_block=false;
	if (a_nnode==1){
		first_block=true;
		// small template 
		a_ins_nvals = _params.NVars;
	}else{
		// full 2-point template
		a_ins_nvals = 2*_params.NVars;
	};
	// large and small index bases
	int l_base = first_block ? 0 : (a_nnode-2)*_params.NVars;
	int s_base=0;
// uniform grid
	const double step = _deta;
	const double inv_step = 1.0/step;
	const double inv_step2 = inv_step*inv_step;
	double f1, f2, f3;
	getMetricCoefs(a_nnode, f1, f2, f3, true);
	if (!first_block){
		for (int i=0; i<_params.NVars; i++){
		// staggered mesh)
		// f~(k-1):
			if (i==2){
				_insert_vals[s_base+i] = 
					getMKLCmplx(0.0, 0.0);
			}else{
				_insert_vals[s_base+i] = getMKLCmplx(
					-f3*inv_step*a_MMat[i][a_eq_id]
					+0.5*a_RMat[i][a_eq_id]);
			};
			_insert_inds[s_base+i] = l_base + i;
		};
		l_base+=_params.NVars;
		s_base+=_params.NVars;
	};
	// f~(k) + f^(k-1/2):
	for (int i=0; i<_params.NVars; i++){
		if (i==2){
			_insert_vals[s_base+i] = getMKLCmplx(a_RMat[i][a_eq_id]);
		}else{
			_insert_vals[s_base+i] = getMKLCmplx(
				f3*inv_step*a_MMat[i][a_eq_id]
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



int t_EigenGS::_solve(){

	static char help[]="Global Search\n";

	wxLogMessage(_("\nGlobal Eigensearch started: N=%d\n"), _params.NNodes);

	int n_inserts = 0;

	// IMPORTANT: reset global matrices -
	// they are corrupted after zggev call...
	for (int i=0; i<_NDIM_G*_NDIM_G; i++){

		_A_G[i] = getMKLCmplx(0.0,0.0);
		_B_G[i] = getMKLCmplx(0.0,0.0);
	}

	// fill A~ by rows
	for (int i=1; i<_params.NNodes; i++){
		for (int j=0; j<_params.NVars; j++){
			if (j==2){
				setMatrices(i, true);
				fill_FO_template(_B, _C, i,j, n_inserts);
			}else{
				setMatrices(i, false);
				fill_SO_template(_A, _B, _C, i,j, n_inserts);
			};
			int row_num = _params.NVars*(i-1)+j;
			// set row values

			insertValsToMat(_A_G, _NDIM_G, row_num, _insert_vals, _insert_inds, n_inserts);

		}

	};
	// fill B~ by rows, order must be the same as in A

	t_SqMatCmplx zero_A(_params.NVars);
	t_SqMatCmplx zero_B(_params.NVars);

	for (int i=1; i<_params.NNodes; i++){
		for (int j=0; j<_params.NVars; j++){
			if (j==2){
				setMatrices(i, true);
				fill_FO_template(zero_B, _CW, i,j, n_inserts);
			}else{
				// optimize setMatrices
				setMatrices(i, false);
				fill_SO_template(zero_A, zero_B, _CW, i,j, n_inserts);
			};
			int row_num = _params.NVars*(i-1)+j;

			insertValsToMat(_B_G, _NDIM_G, row_num, _insert_vals, _insert_inds, n_inserts);

		}
	}

	// matrices A,B filled, ready to compute eigenvalues
	// do not compute left and right eigenvectors
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

		_spectrum[i] = a_s/b_s;

	}



	return 0;
}

int t_EigenGS::getSpectrum(const t_WCharsLoc& init_wave){

	if (init_wave.a.imag()!=0 || init_wave.b.imag()!=0){
		wxLogMessage(_T("Complex wave number provided for Temporal Global Search, ignoring..."));
	}
	
	_alpha = init_wave.a.real();
	_beta  = init_wave.b.real();

	int err_code = _solve();
	return err_code;

};

void t_EigenGS::writeSpectrum(const std::string &a_filename) const{
	std::wofstream to_f(&a_filename[0]);
	for (int i=0; i<_spectrum.size(); i++){
		to_f<<_spectrum[i].real()<<_T("\t")<<_spectrum[i].imag()<<std::endl;
	};
}

void t_EigenGS::writeSpectrumPhase(const std::string &a_filename) const{

	double k = sqrt(_alpha*_alpha + _beta*_beta);
	t_Complex ca, ck;
	std::wofstream to_f(&a_filename[0]);
	for (int i=0; i<_spectrum.size(); i++){
		
		ca = _spectrum[i]/_alpha;
		ck = _spectrum[i]/k;
		
		to_f<<ca.real()<<_T("\t")<<ca.imag()<<_T("\t")
			<<ck.real()<<_T("\t")<<ck.imag()<<_T("\t")<<std::endl;
	};
}

std::vector<t_WCharsLoc> t_EigenGS::getInstabModes(const t_WCharsLoc& init_wave){

		std::vector<t_WCharsLoc> inits;
		std::vector<t_Complex>::const_iterator it;
		getSpectrum(init_wave);
		// TODO: empirics!!!
		for (it=_spectrum.begin(); it<_spectrum.end(); it++){
			if (it->imag()>_params.W_Threshold){
				t_WCharsLoc init_wave;
				init_wave.a = _alpha;
				init_wave.b = _beta;
				init_wave.w = *it;
				init_wave.set_treat(stab::t_TaskTreat::TIME);
				inits.push_back(init_wave);
			}
		}
		return inits;
};

t_WCharsLoc t_EigenGS::searchMaxInstab(const t_WCharsLoc& init_wave){
   const std::vector<t_WCharsLoc>& all_initials = getInstabModes(init_wave);
   return t_WCharsLoc::find_max_instab_time(all_initials);
};


t_WCharsLoc t_EigenGS::searchMaxInstabGlob(){
	//Important TODO: This is the most interesting question : ask AVF
	double a_min = 0.01;
	double a_max = 1.5;
	double b_min = 0.01;
	double b_max = 1.5;
	int n_a = 30;
	int n_b = 30;
	std::vector<t_WCharsLoc> all_initials;
	for (int i=0; i<n_a; i++){
		for(int j=0; j<n_b; j++){
			double a = a_min + (a_max-a_min)/double(n_a)*i;
			double b = b_min + (b_max-b_min)/double(n_b)*j;
			t_WCharsLoc init_wave;
			init_wave.a = a;
			init_wave.b = b;
			std::vector<t_WCharsLoc> inits = getInstabModes(init_wave);
			all_initials.push_back(t_WCharsLoc::find_max_instab_time(inits));
		}
	}
	return t_WCharsLoc::find_max_instab_time(all_initials);
};

std::vector<t_WCharsLoc> t_EigenGS::searchInstabFixed(t_Mode mode, double fixed_val){

	double a,b;
	double *pArg;
	std::vector<t_WCharsLoc> all_initials;
	if (mode==A_MODE){
		b = fixed_val;
		pArg = &a;
	}else{
		a = fixed_val;
		pArg = &b;
	};
	double arg_min = 0.01;
	double arg_max = 1.01;
	int n = 100;
	for (int i=0; i<n; i++){
		std::cout<<"GS fixed: "<<i<<"% done\n";
		*pArg = arg_min + (arg_max-arg_min)/double(n)*i;
		t_WCharsLoc init_wave;
		init_wave.a = a;
		init_wave.b = b;
		init_wave.set_treat(stab::t_TaskTreat::TIME);
		std::vector<t_WCharsLoc> inits = getInstabModes(init_wave);
		all_initials.push_back(t_WCharsLoc::find_max_instab_time(inits));
	};
	return all_initials;
};

t_WCharsLoc t_EigenGS::searchMaxInstabFixed(t_Mode mode, double fixed_val){

	const std::vector<t_WCharsLoc>& all_initials = 
		searchInstabFixed(mode, fixed_val);

	return t_WCharsLoc::find_max_instab_time(all_initials);

};

