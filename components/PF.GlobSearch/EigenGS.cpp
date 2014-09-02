#include "stdafx.h"

#include "common_data.h"

#include "EigenGS.h"
#include "slepceps.h"

#include <iostream>
#include <fstream>

using namespace hsstab;
using namespace pf;

t_EigenGS::t_EigenGS(const mf::t_DomainBase& a_blk):_rBlk(a_blk),
_A(GS_NVARS_TIME),
_B(GS_NVARS_TIME),
_C(GS_NVARS_TIME),
_CW(GS_NVARS_TIME){};

void t_EigenGS::init(const hsstab::TPlugin& g_plug){

	const hsstab::TPluginParamsGroup& g = g_plug.get_settings_grp_const("");

	gs::_init_eigen_gs_base_params(_params, g);

	_init();	// init solver
}


void t_EigenGS::_init(){

	//_A.resize(_params.NVars);
	//_B.resize(_params.NVars);
	//_C.resize(_params.NVars);
	//_CW.resize(_params.NVars);

	SlepcInitialize((int*)0,(char***)0,(char*)0,"EigenGS slepc context");

};

/*
void t_EigenGS::setContext(const mf::t_BlkInd a_ind, 
					  const double a_alpha, const double a_beta,
					  const int a_nnodes){
	if (a_nnodes>0){
		_params.NNodes = a_nnodes;
	}
	_alpha = a_alpha;
	_beta = a_beta;
	t_ProfileNS profNS(_rBlk);

	profNS.initialize(a_ind, _params.ThickCoef);
	_grid.resize(_params.NNodes);
	_profStab.initialize(profNS, _params.NNodes);
	
	double y_max = _profStab.get_thick() - _profStab.get_y(0);
	_a_coef = 1.0*y_max; // play with coef
	_b_coef = 1.0 + _a_coef/y_max;
	double del = 1.0/(double)(_params.NNodes-1);	
	for (int i=0; i<_params.NNodes; i++){
		_grid[i] = (double)(i)*del;
	};
};
*/

void t_EigenGS::setContext(const mf::t_GeomPoint& a_xyz){

	t_ProfileNS profNS(_rBlk);

	// TODO: control number of points in NS profile

	mf::t_ProfDataCfg prof_cfg;
	prof_cfg.ThickCoef = _params.ThickCoef;
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

	_grid.resize(_params.NNodes);
	_profStab.initialize(profNS, _params.NNodes);

	double y_max = _profStab.get_thick() - _profStab.get_y(0);
	_a_coef = 1.0*y_max; // TODO: play with coef
	_b_coef = 1.0 + _a_coef/y_max;
	double del = 1.0/(double)(_params.NNodes-1);	

	for (int i=0; i<_params.NNodes; i++){
		_grid[i] = (double)(i)*del;
	};

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
	double cur_eta = _grid[a_nnode];
	if (a_semi_flag){
		cur_eta-=0.5*(_grid[a_nnode]-_grid[a_nnode-1]);
	};
	const double a_y = _a_coef*cur_eta/(_b_coef-cur_eta);
	t_ProfRec rec = _profStab.get_rec(a_y);
	/*
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
	*/
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
								 const int a_nnode, const int a_eq_id){
	if (a_nnode==0){
		// TODO: to log
		std::wcerr<<_T("GS: SO template on bottom boundary!\n");
	};
	bool first_block=false;
	if (a_nnode==1){
		first_block=true;
		// small template 2 point
		_insert_vals.resize(2*_params.NVars,0.0);
		_insert_inds.resize(2*_params.NVars,0);
	}else{
		if (a_nnode==_params.NNodes-1){
			// "fake" template
			// for diag matrix
			_insert_vals.resize(1,0.0);
			_insert_inds.resize(1,0);
			_insert_vals[0] = 1.0;
			_insert_inds[0] = (a_nnode-1)*_params.NVars+a_eq_id;
			return;
		}else{
			// full 3-point template
			_insert_vals.resize(3*_params.NVars,0.0);
			_insert_inds.resize(3*_params.NVars,0);
		};
		
	};
	// large and small index bases
	int l_base = first_block ? 0 : (a_nnode-2)*_params.NVars;
	int s_base=0;
// uniform grid
	const double step = _grid[a_nnode] - _grid[a_nnode-1];
	const double inv_step = 1.0/step;
	const double inv_step2 = inv_step*inv_step;
	double f1, f2, f3;
	getMetricCoefs(a_nnode, f1, f2, f3, false);
	if (!first_block){
	// first five elements in inserts
	// f~[j-1] : u,v,0,w,t
		for (int i=0; i<_params.NVars; i++){
			if (i==2){
				_insert_vals[s_base+i] = 0.0;
			// f~(k-1):
			}else{
				_insert_vals[s_base+i] = 
					(f1*inv_step2-0.5*f2*inv_step)*a_LMat[i][a_eq_id]
					-0.5*inv_step*f3*a_MMat[i][a_eq_id];
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
			_insert_vals[s_base+i] = 
				-f3*inv_step*a_MMat[i][a_eq_id]
				+0.5*a_RMat[i][a_eq_id];
		}else{
			_insert_vals[s_base+i] = 
				-2.0*f1*inv_step2*a_LMat[i][a_eq_id]
				+a_RMat[i][a_eq_id];
		};	
		_insert_inds[s_base+i] = l_base + i;
	};
	l_base+=_params.NVars;
	s_base+=_params.NVars;
	// fill last five points
	// f~[j+1] + f^[j+1/2]
	for (int i=0; i<_params.NVars; i++){
		if (i==2){
			_insert_vals[s_base+i] = 
				f3*inv_step*a_MMat[i][a_eq_id]
				+0.5*a_RMat[i][a_eq_id];
		}else{
			_insert_vals[s_base+i] = 
				(f1*inv_step2+0.5*f2*inv_step)*a_LMat[i][a_eq_id]
				+0.5*inv_step*f3*a_MMat[i][a_eq_id];
		};	
		_insert_inds[s_base+i] = l_base + i;
	};
	// DEBUG
	double del = 1.0/(double)(_params.NNodes-1);
	for (int i=0; i<_insert_vals.size(); i++){
	//	_insert_vals[i]*=pow(del,2);
	};
	// DEBUG		
	// _insert_vals and _insert_inds filled
};
// we have only one first order continuity equation 
// a_eq_id = 2 
// this template is in staggered point a_nnode - 1/2
void t_EigenGS::fill_FO_template(const t_SqMatCmplx& a_MMat, const t_SqMatCmplx& a_RMat, 
								 const int a_nnode, const int a_eq_id){
	bool first_block=false;
	if (a_nnode==1){
		first_block=true;
		// small template 
		_insert_vals.resize(_params.NVars, 0.0);
		_insert_inds.resize(_params.NVars, 0);
	}else{
		// full 2-point template
		_insert_vals.resize(2*_params.NVars, 0.0);
		_insert_inds.resize(2*_params.NVars, 0);
	};
	// large and small index bases
	int l_base = first_block ? 0 : (a_nnode-2)*_params.NVars;
	int s_base=0;
// uniform grid
	const double step = _grid[a_nnode] - _grid[a_nnode-1];
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
					0.0;
			}else{
				_insert_vals[s_base+i] = 
					-f3*inv_step*a_MMat[i][a_eq_id]
					+0.5*a_RMat[i][a_eq_id];
			};
			_insert_inds[s_base+i] = l_base + i;
		};
		l_base+=_params.NVars;
		s_base+=_params.NVars;
	};
	// f~(k) + f^(k-1/2):
	for (int i=0; i<_params.NVars; i++){
		if (i==2){
			_insert_vals[s_base+i] = a_RMat[i][a_eq_id];
		}else{
			_insert_vals[s_base+i] = 
				f3*inv_step*a_MMat[i][a_eq_id]
				+0.5*a_RMat[i][a_eq_id];
		};
		_insert_inds[s_base+i] = l_base + i;
	};
	// DEBUG
	double del = 1.0/(double)(_params.NNodes-1);
	for (int i=0; i<_insert_vals.size(); i++){
		//_insert_vals[i]*=del;
	}
	// DEBUG
	// filled FO template
};


int t_EigenGS::_solve(){

	static char help[]="Global Search\n";

	Mat         	 A,B;		  
	EPS         	 eps;		  
	const EPSType  type;
	PetscReal   	 error, tol, re, im;
	PetscScalar 	 kr, ki;
	PetscErrorCode ierr;
	PetscInt    	 nev, maxit, Istart, Iend, col[3], its, lits, nconv;
	PetscBool     FirstBlock=PETSC_FALSE, LastBlock=PETSC_FALSE;
	PetscViewer 	 viewer;
	PetscBool  	 flg;
	// to fill matrices
	int mpi_rank, comm_size;

	//SlepcInitialize((int*)0,(char***)0,(char*)0,"SSU - high speed flow stability solver");

	ierr = PetscPrintf(PETSC_COMM_WORLD,"\nGlobal Eigensearch started: N=%d\n\n",_params.NNodes);CHKERRQ(ierr);

	MPI_Comm_size(PETSC_COMM_WORLD, &comm_size);
	MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank);

	if (comm_size!=1){
		std::cerr<<"This is first step: uniprocessor usage only!\n";
		return -1;
	};

	int large_matrix_size = (_params.NNodes-1)*(_params.NVars);

	ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,large_matrix_size, large_matrix_size);CHKERRQ(ierr);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
	ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,large_matrix_size, large_matrix_size);CHKERRQ(ierr);
	ierr = MatSetFromOptions(B);CHKERRQ(ierr);

	// main work goes here
	// these arrays are used to push non-zero elems in
	// matrix rows

	ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);

	//if (Istart<=1) FirstBlock=PETSC_TRUE;
	//if (Iend>=large_matrix_size-n_stab_vars+1) LastBlock=PETSC_TRUE;

	// TODO: decide how to distribute rows to processors
	// for now just fill matrices
	if ((Istart!=0)||(Iend!=large_matrix_size)){
		std::cerr<<"Wrong ownership range-one process; aborting\n";
		return -1;
	}else{
		// fill A~ by rows
		for (int i=1; i<_params.NNodes; i++){
			for (int j=0; j<_params.NVars; j++){
				if (j==2){
					setMatrices(i, true);
					fill_FO_template(_B, _C, i,j);
				}else{
					// optimize setMatrices
					setMatrices(i, false);
					fill_SO_template(_A, _B, _C, i,j);
				};
				int n_inserts = _insert_inds.size();
				int row_num = _params.NVars*(i-1)+j;
				ierr = MatSetValues(A,1,&row_num,n_inserts,&_insert_inds[0],&_insert_vals[0],INSERT_VALUES);CHKERRQ(ierr);
			}

		};
		// fill B~ by rows, order must be the same as in A

		t_SqMatCmplx zero_A(5);
		t_SqMatCmplx zero_B(5);
		for (int i=1; i<_params.NNodes; i++){
			for (int j=0; j<_params.NVars; j++){
				if (j==2){
					setMatrices(i, true);
					fill_FO_template(zero_B, _CW, i,j);
				}else{
					// optimize setMatrices
					setMatrices(i, false);
					fill_SO_template(zero_A, zero_B, _CW, i,j);
				};
				int n_inserts = _insert_inds.size();
				int row_num = _params.NVars*(i-1)+j;
				ierr = MatSetValues(B,1,&row_num,n_inserts,&_insert_inds[0],&_insert_vals[0],INSERT_VALUES);CHKERRQ(ierr);
			}
		}
	}
	// matrices A,B filled, assemble then
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	// view matrices if viewable)
	if (large_matrix_size<200){
		PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
		PetscViewerSetType(viewer, PETSCVIEWERASCII);
		PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_COMMON);
		PetscViewerASCIIOpen(PETSC_COMM_WORLD, "A_VIEW.txt", &viewer);
		ierr = MatView(A, viewer);CHKERRQ(ierr);
		PetscViewerASCIIOpen(PETSC_COMM_WORLD, "B_VIEW.txt", &viewer);
		ierr = MatView(B, viewer);CHKERRQ(ierr);
		PetscViewerDestroy(&viewer);
	}

	//  Create eigensolver context
	ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
	ierr = EPSSetType(eps, EPSARNOLDI);
	ierr = EPSSetDimensions(eps, large_matrix_size, PETSC_DECIDE, PETSC_DECIDE);
	ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD, "matrices assembly: OK");
	//   Set solver parameters at runtime
	//ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);


	ierr = EPSSolve(eps);CHKERRQ(ierr);

	ierr = EPSGetIterationNumber(eps, &its);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %d\n",its);CHKERRQ(ierr);
	ierr = EPSGetOperationCounters(eps,PETSC_NULL,PETSC_NULL,&lits);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of linear iterations of the method: %d\n",lits);CHKERRQ(ierr);
	ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
	ierr = EPSGetDimensions(eps,&nev,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %d\n",nev);CHKERRQ(ierr);
	ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%d\n",tol,maxit);CHKERRQ(ierr);


	ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged approximate eigenpairs: %d\n\n",nconv);CHKERRQ(ierr);
	_spectrum.resize(nconv, 0.0);
	if (nconv>0) {
		//ierr = PetscPrintf(PETSC_COMM_WORLD,
		// "           k             ||Ax-kBx||/||kx||\n"
		// "  --------------------- ------------------\n" );CHKERRQ(ierr);

		for(PetscInt i=0; i<nconv; i++ ) {
			ierr = EPSGetEigenpair(eps,i,&kr,&ki,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
			ierr = EPSComputeRelativeError(eps,i,&error);CHKERRQ(ierr);

#if defined(PETSC_USE_COMPLEX)
			re = PetscRealPart(kr);
			im = PetscImaginaryPart(kr);
			_spectrum[i] = std::complex<double>(re, im);
#else
			re = kr;
			im = ki;
#endif
			//ierr = PetscPrintf(PETSC_COMM_WORLD," % 12g\n",error);CHKERRQ(ierr);

		}
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n" );CHKERRQ(ierr);
	}


	ierr = EPSDestroy(&eps);CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = MatDestroy(&B);CHKERRQ(ierr);
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

void t_EigenGS::writeSpectrum(const std::wstring &a_filename) const{
	std::wofstream to_f(&a_filename[0]);
	for (int i=0; i<_spectrum.size(); i++){
		to_f<<_spectrum[i].real()<<_T("\t")<<_spectrum[i].imag()<<std::endl;
	};
}

void t_EigenGS::writeSpectrumPhase(const std::wstring &a_filename) const{

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

