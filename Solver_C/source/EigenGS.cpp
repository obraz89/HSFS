#include "EigenGS.h"
#include "slepceps.h"
// task dim =5 3D
// 4 for 2D
t_EigenGS::t_EigenGS(const MF_Field& a_rFld, const int a_task_dim):
_rFldNS(a_rFld), _n_vars(a_task_dim), _profStab(0),
_A(a_task_dim), _B(a_task_dim), _C(a_task_dim), _CW(a_task_dim){
	// keep max size : SO template 15 point for 3D
};

void t_EigenGS::setContext(const int a_i, const int a_k, 
					  const double& a_alpha, const double& a_beta,
					  const int a_nnodes){
	_alpha = a_alpha;
	_beta = a_beta;
	t_ProfileNS profNS(_rFldNS);
	profNS.setProfiles(a_i, a_k);
	_nnodes = a_nnodes;
	_profStab.resize(a_nnodes);
	_grid.resize(a_nnodes);

	_profStab.setProfiles(profNS);
	
	double y_max = _profStab.y.back() - _profStab.y[0];
	_a_coef = 1.0*y_max; // play with coef
	_b_coef = 1.0 + _a_coef/y_max;
	double del = 1.0/(double)(_nnodes-1);	
	for (int i=0; i<_nnodes; i++){
		_grid[i] = (double)(i)*del;
	};
};
// semi-flag is true if it is k-1/2 point
void t_EigenGS::getMetricCoefs(const int& a_nnode, double& f1, double& f2, double& f3, const bool semi_flag) const{
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

	const double& stabRe = _profStab.stabRe;
	const double& Me = _profStab.Me;
	const double gMaMa = MF_Field::Gamma*Me*Me;
	const double g_1MaMa = (MF_Field::Gamma-1.0)*Me*Me;

	// rename a_y is not input but lazy to rewrite
	// get physical y
	double cur_eta = _grid[a_nnode];
	if (a_semi_flag){
		cur_eta-=0.5*(_grid[a_nnode]-_grid[a_nnode-1]);
	};
	const double a_y = _a_coef*cur_eta/(_b_coef-cur_eta);

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
	// k''/k:
	const double mu_coef = mu1*t1*inv_mu;
	const double k2k = 2.0*pow(mu_coef,2)-mu_coef;

	const double dzeta_noW = _alpha*u + _beta*w;
	const double dzeta_W = -1.0;

// set _A
	_A.setToUnity();
	_A[2][2]=0.0;
// set _B
	const double lz = 0.0; // 0 + lambda/mu
	const double lf = 1.0; // 1 + lambda/mu
	const double ls = 2.0; // 2 + lambda/mu

	// first row
	_B[0][0] = inv_mu*mu1*t1;
	_B[1][0] = imagUnity*_alpha*lf;
	_B[3][0] = inv_mu*mu1*u1;
	//second row
	_B[0][1] = imagUnity*_alpha*lf/ls;
	_B[1][1] = _B[0][0];
	_B[2][1] = -stabRe/(ls*mu);
	_B[4][1] = imagUnity*_beta*lf/ls;
	// third
	_B[1][2] = 1.0;
	// fourth
	_B[0][3] = 2.0*g_1MaMa*MF_Field::Pr*u1;
	// k'/k ??
	_B[3][3] = -2.0*mu1*t1*inv_mu;
	_B[4][3] = 2.0*g_1MaMa*MF_Field::Pr*w1;
	// last
	_B[1][4] = imagUnity*_beta*lf;
	_B[3][4] = inv_mu*mu1*w1;
	_B[4][4] = inv_mu*mu1*t1;
// clear C - no freq terms
	// first row
	_C[0][0] = -imagUnity*dzeta_noW*stabRe*inv_mu*inv_t
		       -(ls*pow(_alpha,2)+pow(_beta,2));

	_C[1][0] = -stabRe*u1*inv_mu*inv_t
		       +imagUnity*_alpha*inv_mu*mu1*t1;

	_C[2][0] = -imagUnity*_alpha*stabRe*inv_mu;

	_C[3][0] = inv_mu*(mu1*u2 + mu2*t1*u1);
	_C[4][0] = -_alpha*_beta*lf;
	// second

	_C[0][1] = imagUnity*_alpha*inv_mu*mu1*t1*lz/ls;

	_C[1][1] = -imagUnity*dzeta_noW*stabRe/(ls*mu*t)
		       -(pow(_alpha,2)+pow(_beta,2))/ls;

	_C[3][1] = imagUnity*inv_mu*mu1*(_alpha*u1+_beta*w1)/ls;
	_C[4][1] = imagUnity*_beta*inv_mu*mu1*t1*lz/ls;
	// third
	_C[0][2] = imagUnity*_alpha;
	_C[1][2] = -t1*inv_t;
	_C[2][2] = imagUnity*gMaMa*dzeta_noW;
	_C[3][2] = -imagUnity*dzeta_noW*inv_t;
	_C[4][2] = imagUnity*_beta;
	//fourh
	_C[1][3] = MF_Field::Pr*
			   (
			    2.0*imagUnity*g_1MaMa*(_alpha*u1 + _beta*w1)
			    -stabRe*t1*inv_mu*inv_t);

	_C[2][3] = imagUnity*dzeta_noW*g_1MaMa*MF_Field::Pr
		       *stabRe*inv_mu;

	_C[3][3] = -imagUnity*dzeta_noW*stabRe*MF_Field::Pr*inv_mu*inv_t
		       -(pow(_alpha,2)+pow(_beta,2))
			   +g_1MaMa*MF_Field::Pr*inv_mu*mu1*(pow(u1,2)+pow(w1,2))
			   +k2k;
	// last at least
	_C[0][4] = -_alpha*_beta*lf;
	_C[1][4] = imagUnity*_beta*inv_mu*mu1*t1
		       -stabRe*w1*inv_mu*inv_t;

	_C[2][4] = -imagUnity*_beta*stabRe*inv_mu;
	_C[3][4] = inv_mu*mu1*w2 + inv_mu*mu2*t1*w1;

	_C[4][4] = -imagUnity*dzeta_noW*stabRe*inv_mu*inv_t
		       -(pow(_alpha,2)+ls*pow(_beta,2));
	// W matrix
	_CW[0][0] = -imagUnity*dzeta_W*stabRe*inv_mu*inv_t;
	_CW[1][1] = -imagUnity*dzeta_W*stabRe/ls*inv_mu*inv_t;
	_CW[2][2] = imagUnity*gMaMa*dzeta_W;
	_CW[3][2] = -imagUnity*dzeta_W*inv_t;
	_CW[2][3] = imagUnity*dzeta_W*g_1MaMa*MF_Field::Pr*stabRe*inv_mu;
	_CW[3][3] = -imagUnity*dzeta_W*stabRe*MF_Field::Pr*inv_mu*inv_t;
	_CW[4][4] = -imagUnity*dzeta_W*stabRe*inv_mu*inv_t;

};
// a_eq_id={0,1,3,4} <--> {1,2,4,5} - SO
void t_EigenGS::fill_SO_template(const t_SqMatrix& a_LMat, const t_SqMatrix& a_MMat, const t_SqMatrix& a_RMat, 
								 const int a_nnode, const int a_eq_id){
	if (a_nnode==0){
		std::cerr<<"GS: SO template on bottom boundary!\n";
	};
	bool first_block=false;
	if (a_nnode==1){
		first_block=true;
		// small template 2 point
		_insert_vals.resize(2*_n_vars,0.0);
		_insert_inds.resize(2*_n_vars,0);
	}else{
		if (a_nnode==_nnodes-1){
			// "fake" template
			// for diag matrix
			_insert_vals.resize(1,0.0);
			_insert_inds.resize(1,0);
			_insert_vals[0] = 1.0;
			_insert_inds[0] = (a_nnode-1)*_n_vars+a_eq_id;
			return;
		}else{
			// full 3-point template
			_insert_vals.resize(3*_n_vars,0.0);
			_insert_inds.resize(3*_n_vars,0);
		};
		
	};
	// large and small index bases
	int l_base = first_block ? 0 : (a_nnode-2)*_n_vars;
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
		for (int i=0; i<_n_vars; i++){
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
		l_base+=_n_vars;
		s_base+=_n_vars;
	};
	
	// fill central 5 points: 
	// f~[j] + f^[j-1/2]
	for (int i=0; i<_n_vars; i++){
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
	l_base+=_n_vars;
	s_base+=_n_vars;
	// fill last five points
	// f~[j+1] + f^[j+1/2]
	for (int i=0; i<_n_vars; i++){
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
	double del = 1.0/(double)(_nnodes-1);
	for (int i=0; i<_insert_vals.size(); i++){
	//	_insert_vals[i]*=pow(del,2);
	};
	// DEBUG		
	// _insert_vals and _insert_inds filled
};
// we have only one first order continuity equation 
// a_eq_id = 2 
// this template is in staggered point a_nnode - 1/2
void t_EigenGS::fill_FO_template(const t_SqMatrix& a_MMat, const t_SqMatrix& a_RMat, 
								 const int a_nnode, const int a_eq_id){
	bool first_block=false;
	if (a_nnode==1){
		first_block=true;
		// small template 
		_insert_vals.resize(_n_vars, 0.0);
		_insert_inds.resize(_n_vars, 0);
	}else{
		// full 2-point template
		_insert_vals.resize(2*_n_vars, 0.0);
		_insert_inds.resize(2*_n_vars, 0);
	};
	// large and small index bases
	int l_base = first_block ? 0 : (a_nnode-2)*_n_vars;
	int s_base=0;
// uniform grid
	const double step = _grid[a_nnode] - _grid[a_nnode-1];
	const double inv_step = 1.0/step;
	const double inv_step2 = inv_step*inv_step;
	double f1, f2, f3;
	getMetricCoefs(a_nnode, f1, f2, f3, true);
	if (!first_block){
		for (int i=0; i<_n_vars; i++){
		// staggered mesh)
		// f~(k-1):
			if (i==2){
				_insert_vals[s_base+i] = 
					0.0;
			}else{
				_insert_vals[s_base+i] = 
					-f3*inv_step*a_MMat[i][a_eq_id];
			};
			_insert_inds[s_base+i] = l_base + i;
		};
		l_base+=_n_vars;
		s_base+=_n_vars;
	};
	// f~(k) + f^(k-1/2):
	for (int i=0; i<_n_vars; i++){
		if (i==2){
			_insert_vals[s_base+i] = a_RMat[i][a_eq_id];
		}else{
			_insert_vals[s_base+i] = 
				f3*inv_step*a_MMat[i][a_eq_id];
		};
		_insert_inds[s_base+i] = l_base + i;
	};
	// DEBUG
	double del = 1.0/(double)(_nnodes-1);
	for (int i=0; i<_insert_vals.size(); i++){
		//_insert_vals[i]*=del;
	}
	// DEBUG
	// filled FO template
};

int t_EigenGS::search(){
  static char help[]="Global Search\n";
  // slepc locals
  Mat         	 A,B;		  
  EPS         	 eps;		  
  const EPSType  type;
  PetscReal   	 error, tol, re, im;
  PetscScalar 	 kr, ki;
  PetscErrorCode ierr;
  PetscInt    	 nev, maxit, Istart, Iend, col[3], its, lits, nconv;
  PetscTruth     FirstBlock=PETSC_FALSE, LastBlock=PETSC_FALSE;
  PetscViewer 	 viewer;
  PetscTruth  	 flg;
  // to fill matrices
  int mpi_rank, comm_size;

  // pass command line arguments
  SlepcInitialize((int*)0,(char***)0,(char*)0,help);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nGlobal Eigensearch started: N=%d\n\n",_nnodes);CHKERRQ(ierr);

  MPI_Comm_size(PETSC_COMM_WORLD, &comm_size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank);

  if (comm_size!=1){
	  std::cerr<<"This is first step: uniprocessor usage only!\n";
	  return -1;
  };

  int large_matrix_size = (_nnodes-1)*(_n_vars);

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
	    for (int i=1; i<_nnodes; i++){
		  for (int j=0; j<_n_vars; j++){
			  if (j==2){
				  setMatrices(i, true);
				fill_FO_template(_B, _C, i,j);
			  }else{
				  // optimize setMatrices
				  setMatrices(i, false);
				fill_SO_template(_A, _B, _C, i,j);
			  };
			int n_inserts = _insert_inds.size();
			int row_num = _n_vars*(i-1)+j;
			ierr = MatSetValues(A,1,&row_num,n_inserts,&_insert_inds[0],&_insert_vals[0],INSERT_VALUES);CHKERRQ(ierr);
		}
		
	    };
	  // fill B~ by rows, order must be the same as in A

		t_SqMatrix zero_A(5);
		t_SqMatrix zero_B(5);
	  for (int i=1; i<_nnodes; i++){
		  for (int j=0; j<_n_vars; j++){
			  if (j==2){
				  setMatrices(i, true);
				fill_FO_template(zero_B, _CW, i,j);
			  }else{
				  // optimize setMatrices
				  setMatrices(i, false);
				fill_SO_template(zero_A, zero_B, _CW, i,j);
			  };
			int n_inserts = _insert_inds.size();
			int row_num = _n_vars*(i-1)+j;
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
		PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);
		PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_COMMON);
		PetscViewerASCIIOpen(PETSC_COMM_WORLD, "A_VIEW.txt", &viewer);
		ierr = MatView(A, viewer);CHKERRQ(ierr);
		PetscViewerASCIIOpen(PETSC_COMM_WORLD, "B_VIEW.txt", &viewer);
		ierr = MatView(B, viewer);CHKERRQ(ierr);
		PetscViewerDestroy(viewer);
	}

	//  Create eigensolver context
	ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
	//ierr = EPSSetType(eps, EPSARNOLDI);
	ierr = EPSSetDimensions(eps, large_matrix_size, PETSC_DECIDE, PETSC_DECIDE);
	ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr);

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

	if (nconv>0) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,
         "           k             ||Ax-kBx||/||kx||\n"
         "  --------------------- ------------------\n" );CHKERRQ(ierr);
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
			if( im != 0.0 ) {
				ierr = PetscPrintf(PETSC_COMM_WORLD," % 6f %+6f i",re,im);CHKERRQ(ierr);
			} else {
				ierr = PetscPrintf(PETSC_COMM_WORLD,"       % 6f      ",re); CHKERRQ(ierr);
			}
			ierr = PetscPrintf(PETSC_COMM_WORLD," % 12g\n",error);CHKERRQ(ierr);
		}
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n" );CHKERRQ(ierr);
	}
  
  ierr = EPSDestroy(eps);CHKERRQ(ierr);
  ierr = MatDestroy(A);CHKERRQ(ierr);
  ierr = MatDestroy(B);CHKERRQ(ierr);
  ierr = SlepcFinalize();CHKERRQ(ierr);
  return 0;

};
