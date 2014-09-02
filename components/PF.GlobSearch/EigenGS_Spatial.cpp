#include "stdafx.h"

#include "common_data.h"

#include "EigenGS_Spatial.h"
#include "slepceps.h"

#include <iostream>
#include <fstream>

using namespace hsstab;
using namespace pf;

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

	gs::_init_eigen_gs_base_params(_params, g);

	_init();	// init solver
}


void t_GlobSrchSpat::_init(){

	SlepcInitialize((int*)0,(char***)0,(char*)0,"EigenGS slepc context");

};

void t_GlobSrchSpat::setContext(const mf::t_GeomPoint& a_xyz){

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
	// TODO: play with coef a
	// a=1*y_max : half of nodes inside 1/3 of comp domain
	// a=infinity*y_max : even physical distribution of nodes
	_a_coef = 1.0*y_max;
	_b_coef = 1.0 + _a_coef/y_max;
	_deta = 1.0/(double)(_params.NNodes-1);	

	for (int i=0; i<_params.NNodes; i++) _grid[i] = i*_deta;

}

void t_GlobSrchSpat::setContext(const t_ProfileStab* a_prof_stab){
	
	_grid.resize(_params.NNodes);

	_profStab = *a_prof_stab;

	double y_max = _profStab.get_thick() - _profStab.get_y(0);

	// this is what Malik advises
	double y_semi = _profStab.get_y_by_velo(0.5);

	_a_coef = 2.0*y_semi; // play with coef
	_b_coef = 1.0 + _a_coef/y_max;
	double del = 1.0/(double)(_params.NNodes-1);	
	for (int i=0; i<_params.NNodes; i++){
		_grid[i] = (double)(i)*del;
	};

}

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

	// rename a_y is not input but lazy to rewrite
	// get physical y
	double cur_eta = _grid[a_nnode];
	if (a_semi_flag){
		cur_eta+=0.5*(_grid[a_nnode+1]-_grid[a_nnode]);
	};
	const double a_y = _a_coef*cur_eta/(_b_coef-cur_eta);
	t_ProfRec rec = _profStab.get_rec(a_y);

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
	// first row
	_C[0][0] = -imagUnity*stabRe*dzeta_noA*inv_mu*inv_t-_beta*_beta;

	_C[1][0] = -stabRe*U1*inv_mu*inv_t; 

	_C[3][0] = inv_mu*(mu1*U2 + mu2*T1*U1);

	// second row
	_C[1][1] = -imagUnity*dzeta_noA*stabRe/(L2*mu*T)-pow(_beta,2)/L2;

	_C[3][1] = imagUnity*inv_mu*mu1*_beta*W1/L2;

	_C[4][1] = imagUnity*_beta*inv_mu*mu1*T1*L0/L2;

	// third

	_C[1][2] = -T1*inv_t;

	_C[2][2] = imagUnity*gMaMa*dzeta_noA;

	_C[3][2] = -imagUnity*dzeta_noA*inv_t;

	_C[4][2] = imagUnity*_beta;

	//fourth
	_C[1][3] = Pr*(2.0*imagUnity*g_1MaMa*_beta*W1 - stabRe*T1*inv_mu*inv_t);

	_C[2][3] = imagUnity*dzeta_noA*g_1MaMa*Pr*stabRe*inv_mu;

	_C[3][3] = -imagUnity*dzeta_noA*stabRe*Pr*inv_mu*inv_t - pow(_beta,2)
			   +g_1MaMa*Pr*inv_mu*mu1*(pow(U1,2)+pow(W1,2)) + k2k;

	// last at least
	_C[1][4] = imagUnity*_beta*inv_mu*mu1*T1 - stabRe*W1*inv_mu*inv_t;

	_C[2][4] = -imagUnity*_beta*stabRe*inv_mu;

	_C[3][4] = inv_mu*mu1*W2 + inv_mu*mu2*T1*W1;

	_C[4][4] = -imagUnity*dzeta_noA*stabRe*inv_mu*inv_t - L2*pow(_beta,2);
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
// a_eq_id={0,1,3,4} <--> {1,2,4,5} - SO
void t_GlobSrchSpat::fill_SO_row(const t_SqMatCmplx& a_LMat, const t_SqMatCmplx& a_MMat, 
								 const t_SqMatCmplx& a_RMat, const int a_nnode, const int a_eq_id){

		
	// large and small index bases
	int l_base = (a_nnode-1)*GS_NVARS_SPAT;
	int s_base=0;
// uniform grid
	const double step = _deta;
	const double inv_step = 1.0/step;
	const double inv_step2 = inv_step*inv_step;
	double f1, f2, f3;
	getMetricCoefs(a_nnode, f1, f2, f3, false);

	// first five elements in inserts: f~[j-1] and f^[j-1/2]
	// f~<=>{u,v,o,w,t} ; f^<=>{0,0,p,0,0}
	for (int i=0; i<GS_NVARS_SPAT; i++){
		if (i==2){
			_insert_vals_SO[s_base+i] = -f3*inv_step*a_MMat[i][a_eq_id]
									 +0.5*a_RMat[i][a_eq_id];
		}else{
			_insert_vals_SO[s_base+i] = 
				(f1*inv_step2-0.5*f2*inv_step)*a_LMat[i][a_eq_id]
				-0.5*inv_step*f3*a_MMat[i][a_eq_id];
		};
		_insert_inds_SO[s_base+i] = l_base + i;
	};
	l_base+=_params.NVars;
	s_base+=_params.NVars;
	
	// fill central 5 points:  f~[j] and f^[j+1/2]
	for (int i=0; i<_params.NVars; i++){
		if (i==2){
			_insert_vals_SO[s_base+i] = f3*inv_step*a_MMat[i][a_eq_id]
									 +0.5*a_RMat[i][a_eq_id];
		}else{
			_insert_vals_SO[s_base+i] = -2.0*f1*inv_step2*a_LMat[i][a_eq_id]
									 +a_RMat[i][a_eq_id];
		};	
		_insert_inds_SO[s_base+i] = l_base + i;
	};
	l_base+=_params.NVars;
	s_base+=_params.NVars;
	// fill last five points : f~[j+1] + f^[j+3/2]
	for (int i=0; i<_params.NVars; i++){
		if (i==2){
			_insert_vals_SO[s_base+i] = 0.0;
		}else{
			_insert_vals_SO[s_base+i] = 
				(f1*inv_step2+0.5*f2*inv_step)*a_LMat[i][a_eq_id]
				+0.5*inv_step*f3*a_MMat[i][a_eq_id];
		};	
		_insert_inds_SO[s_base+i] = l_base + i;
	};
};
// the only first order continuity equation 
void t_GlobSrchSpat::fill_FO_row(const t_SqMatCmplx& a_MMat, const t_SqMatCmplx& a_RMat, 
								 const int a_nnode){

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

	// first block f~[j] and f^[j+1/2]:
	for (int i=0; i<_params.NVars; i++){
		if (i==2){
			_insert_vals_FO[s_base+i] = a_RMat[i][a_eq_id];
		}else{
			_insert_vals_FO[s_base+i] = -f3*inv_step*a_MMat[i][a_eq_id]
									 +0.5*a_RMat[i][a_eq_id];
		};
		_insert_inds_FO[s_base+i] = l_base + i;
	};
	l_base+=_params.NVars;
	s_base+=_params.NVars;

	// second block f~[j+1] + f^[j+3/2]:
	for (int i=0; i<_params.NVars; i++){
		if (i==2){
			_insert_vals_FO[s_base+i] = 0.0;
		}else{
			_insert_vals_FO[s_base+i] = f3*inv_step*a_MMat[i][a_eq_id]
										+0.5*a_RMat[i][a_eq_id];
		};
		_insert_inds_FO[s_base+i] = l_base + i;
	};
};


int t_GlobSrchSpat::_solve(){

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

	// one value - p^ at N+1/2 is added as arbitrary value
	// to keep data structure NNodes*{u,v,p,t,w}
	int large_matrix_size = _params.NNodes*GS_NVARS_SPAT;

	ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,large_matrix_size, large_matrix_size);CHKERRQ(ierr);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
	ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,large_matrix_size, large_matrix_size);CHKERRQ(ierr);
	ierr = MatSetFromOptions(B);CHKERRQ(ierr);

	ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);

	// TODO: decide how to distribute rows to processors
	// for now just fill matrices
	if ((Istart!=0)||(Iend!=large_matrix_size)){
		std::cerr<<"Wrong ownership range-one process; aborting\n";
		return -1;
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
				t_Complex insert_val_A[1];insert_val_A[0]=1.0; 
				t_Complex insert_val_B[1];insert_val_B[0]=1.0;
				int insert_ind[1]; insert_ind[0]=row_num;
				MatSetValues(A,1,&row_num,1,&insert_ind[0], &insert_val_A[0],INSERT_VALUES);CHKERRQ(ierr);
				MatSetValues(B,1,&row_num,1,&insert_ind[0], &insert_val_B[0],INSERT_VALUES);CHKERRQ(ierr);

			}else{

				// ordinary FO row for p^[0+1/2]
				setMatrices(0, true);
				fill_FO_row(_B, _C, 0);
				MatSetValues(A,1,&row_num,GS_FO_INS_SIZE,&_insert_inds_FO[0], &_insert_vals_FO[0],INSERT_VALUES);CHKERRQ(ierr);

				fill_FO_row(_B_AL, _C_AL, 0);
				MatSetValues(B,1,&row_num,GS_FO_INS_SIZE,&_insert_inds_FO[0], &_insert_vals_FO[0],INSERT_VALUES);CHKERRQ(ierr);

			}
		}
	}
	// main body of the matrices A, B - {SO, SO, FO, SO, SO} blocks
	for (int i=1; i<_params.NNodes-1; i++){
		for (int j=0; j<GS_NVARS_SPAT; j++){
			row_num = GS_NVARS_SPAT*i + j;
			if (j==2){
				setMatrices(i, true);

				fill_FO_row(_B, _C, i);
				MatSetValues(A,1,&row_num,GS_FO_INS_SIZE,&_insert_inds_FO[0], &_insert_vals_FO[0],INSERT_VALUES);CHKERRQ(ierr);

				fill_FO_row(_B_AL, _C_AL, i);
				MatSetValues(B,1,&row_num,GS_FO_INS_SIZE,&_insert_inds_FO[0], &_insert_vals_FO[0],INSERT_VALUES);CHKERRQ(ierr);

			}else{
				setMatrices(i, false);

				fill_SO_row(_A, _B, _C,i, j);
				MatSetValues(A,1,&row_num,GS_SO_INS_SIZE,&_insert_inds_SO[0], &_insert_vals_SO[0],INSERT_VALUES);CHKERRQ(ierr);

				fill_SO_row(_A_AL, _B_AL, _C_AL,i, j);
				MatSetValues(B,1,&row_num,GS_SO_INS_SIZE,&_insert_inds_SO[0], &_insert_vals_SO[0],INSERT_VALUES);CHKERRQ(ierr);
			};
		}
	};
	// last block 5x5
	// just set blocks to unity matrices
	{
		for (int j=0; j<GS_NVARS_SPAT; j++){

			row_num=(_params.NNodes-1)*GS_NVARS_SPAT + j;
			t_Complex insert_val_A[1];insert_val_A[0]=1.0; 
			t_Complex insert_val_B[1];insert_val_B[0]=1.0;
			int insert_ind[1]; insert_ind[0]=row_num;

			MatSetValues(A,1,&row_num,1,&insert_ind[0], &insert_val_A[0],INSERT_VALUES);CHKERRQ(ierr);
			MatSetValues(B,1,&row_num,1,&insert_ind[0], &insert_val_B[0],INSERT_VALUES);CHKERRQ(ierr);

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
#else
			re = kr;
			im = ki;
#endif
			// system solved in for Ah+lambda*Bh=0
			// so alpha = -lambda
			_spectrum[i] = std::complex<double>(-1.0*re, -1.0*im);
			//ierr = PetscPrintf(PETSC_COMM_WORLD," % 12g\n",error);CHKERRQ(ierr);

		}
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n" );CHKERRQ(ierr);
	}


	ierr = EPSDestroy(&eps);CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = MatDestroy(&B);CHKERRQ(ierr);
	return 0;
}

int t_GlobSrchSpat::getSpectrum(const t_WCharsLoc& init_wave){

  _beta  = init_wave.b;
  _w = init_wave.w;

  int err_code = _solve();
  return err_code;

};

void t_GlobSrchSpat::writeSpectrum(const std::wstring &a_filename) const{
	std::wofstream to_f(&a_filename[0]);
	for (int i=0; i<_spectrum.size(); i++){
		to_f<<_spectrum[i].real()<<_T("\t")<<_spectrum[i].imag()<<std::endl;
	};
}

void t_GlobSrchSpat::writeSpectrumPhase(const std::wstring &a_filename) const{

	wxLogMessage(_T("writeSpectrumPhase not implemented for Spatial GS!\n"));
	/*
	double k = sqrt(_alpha*_alpha + _beta*_beta);
	t_Complex ca, ck;
	std::wofstream to_f(&a_filename[0]);
	for (int i=0; i<_spectrum.size(); i++){
		
		ca = _spectrum[i]/_alpha;
		ck = _spectrum[i]/k;
		
		to_f<<ca.real()<<_T("\t")<<ca.imag()<<_T("\t")
			<<ck.real()<<_T("\t")<<ck.imag()<<_T("\t")<<std::endl;
	};
	*/
}

std::vector<t_WCharsLoc> t_GlobSrchSpat::getInstabModes(const t_WCharsLoc& init_wave){

		std::vector<t_WCharsLoc> inits;
		std::vector<t_Complex>::const_iterator it;
		getSpectrum(init_wave);
		// TODO: empirics!!!
		for (it=_spectrum.begin(); it<_spectrum.end(); it++){
			// IMPORTANT TODO: fix parasitic solutions, ask AVF
			// now neglect parasitic branch manually
			if (it->imag()<-1.0*_params.W_Threshold && it->real()>0.0){
				t_WCharsLoc wave;
				wave.a = *it;
				wave.b = _beta;
				wave.w = _w;
				wave.set_treat(stab::t_TaskTreat::SPAT);
				inits.push_back(wave);
			}
		}
		return inits;
};

t_WCharsLoc t_GlobSrchSpat::searchMaxInstab(const t_WCharsLoc& init_wave){
	const std::vector<t_WCharsLoc>& all_initials = getInstabModes(init_wave);
	return t_WCharsLoc::find_max_instab_spat(all_initials);
};
