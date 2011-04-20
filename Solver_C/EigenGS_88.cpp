#include "EigenGS.h"
#include "slepceps.h"
t_EigenGS::t_EigenGS(const t_StabSolver &a_stab_slv):
_rStab_slvr(a_stab_slv), _y_range(_rStab_slvr._profStab.y), 
_step(_rStab_slvr._profStab.y[1] - _rStab_slvr._profStab.y[0]){
  if (_rStab_slvr._math_solver.getTaskDim()==4){
	  _n_vars = 5;
  }else{
	  if (_rStab_slvr._math_solver.getTaskDim()==3){
		  _n_vars = 4;
	  }else{
		  std::cerr<<"GSearch ERROR: Unknown task dimension\n";
	  };
  };
  _nnodes = _rStab_slvr._profStab.size();
  _vec_pack_order.resize(_nnodes);
  _vec_pack_order[0]=1;
  _vec_pack_order.back()=1;
  for (int i=1; i<_nnodes-1; i++){
		_vec_pack_order[i] = _n_indep_vars;
  };
}

int t_EigenGS::getVectorIndex(const int a_j, const int a_k){
	int ind=0;
	for (int i=0;i<a_j;i++){
		ind+=_vec_pack_order[i];
	}
	ind+=a_k;
	return ind;
}
	/*
	This subroutine implements modified Spectral Search
	procedure first described by Malik, 1980 [ref.doc]
	The goal is to achieve initia guesses in time approach
	for a complex wavenumber w with a given pair of
	wavenumbers (alpha,beta)
	the problem has a form
	Ax=wBx
	x-eigenvectors, w-set of eigenvalues
	packing
	EigenVector:
	from stability computations the vector of solution
	at a point is {u, u', v, p, w, w', t, t'}
	to form eigenproblem the uniform grid of a_nnodes nodes is
	utilized to represent the full discrete vector of solution
	with reduced set of vars s={u,v,p,w,t}
	{s0,s1,s2,s[a_nnodes-1]}
	for a global search boundary conditions
	{u,v,w,t}wall={u,v,w,t}outer={0,0,0,0}
	thus for computations vector x is
	{p[0], s[1], ..., s[n-2], p[n-1]}
	this order is used to fill non-zero elements
	in A and B which are of size
	5(n-2)+2 for 3D and 4(n-2)+2 for 2D
	fill order: as in stab comps: 1,3,4,5,7 <--> 1,2,3,4,5

    */

int t_EigenGS::search(){
  static char help[]="Global Search\n";
  int _task_dim = _rStab_slvr.getTaskDim();
  // split matrices are packed 
  // and have sizes 5rowsX8cols for 3D
  // 2D: 4rowsX6cols
  t_Matrix stab_mat_no_w(2*_task_dim,_task_dim+1);
  t_Matrix stab_mat_w(2*_task_dim,_task_dim+1);
  // slepc locals
  Mat         	 A,B;		  /* matrices */
  EPS         	 eps;		  /* eigenproblem solver context */
  const EPSType  type;
  PetscReal   	 error, tol, re, im;
  PetscScalar 	 kr, ki;
  PetscErrorCode ierr;
  PetscInt    	 nev, maxit, i, Istart, Iend, col[3], its, lits, nconv;
  PetscTruth     FirstBlock=PETSC_FALSE, LastBlock=PETSC_FALSE;
  PetscViewer 	 viewer;
  PetscTruth  	 flg;
  // to fill matrices
  int mpi_rank, comm_size;

  // pass command line arguments
  SlepcInitialize((int*)0,(char***)0,(char*)0,help);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nGlobal Eigensearch started%d\n\n",_nnodes);CHKERRQ(ierr);

  MPI_Comm_size(PETSC_COMM_WORLD, &comm_size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank);

  if (comm_size!=1){
	  std::cerr<<"This is first step: uniprocessor usage only!\n";
	  return -1;
  }

  int large_matrix_size = _n_vars*(_nnodes-2)+2;

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,large_matrix_size, large_matrix_size);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
  ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,large_matrix_size, large_matrix_size);CHKERRQ(ierr);

  // main work goes here
  // these arrays are used to push non-zero elems in
  // matrix rows
  double insert_vals[16];
  int insert_ind[16];
  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);

  //if (Istart<=1) FirstBlock=PETSC_TRUE;
  //if (Iend>=large_matrix_size-n_stab_vars+1) LastBlock=PETSC_TRUE;

  // TODO: decide how to distribute rows to processors
  // for now just fill matrices
  if ((Istart!=0)||(Iend!=large_matrix_size)){
	  std::cerr<<"Wrong ownership range-one process; aborting\n";
	  return -1;
  }else{
	  double cur_y = _y_range[0];
	  double inv_step = 1.0/step;
	  _rStab_slvr.setMatSplitByW3D(cur_y, stab_mat_no_w, stab_mat_w);
	  // fill A
	  // first string boundary-wall k="3":
	  // no freq: 5 non-zeros
	  // h[k][1]
	  // это тупо но ты уже придумал как сделать умно)))
	  insert_vals[0]=stab_mat_no_w[0][1];
	  insert_vals[1]=stab_mat_no_w[2][1]+inv_step;
	  insert_vals[2]=stab_mat_no_w[3][1]
	  for (int i=Istart;i<Iend; i++){

	  }
  }
  // matrices A,B filled

  //  Create eigensolver context
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr);

  //   Set solver parameters at runtime
  //ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  ierr = EPSSolve(eps);CHKERRQ(ierr);

  /*
     Optional: Get some information from the solver and display it
  */
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

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* 
     Get number of converged eigenpairs
  */
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged approximate eigenpairs: %d\n\n",nconv);CHKERRQ(ierr);

  if (nconv>0) {
    /*
       Display eigenvalues and relative errors
    */
    ierr = PetscPrintf(PETSC_COMM_WORLD,
         "           k             ||Ax-kBx||/||kx||\n"
         "  --------------------- ------------------\n" );CHKERRQ(ierr);
    for( i=0; i<nconv; i++ ) {
      /* 
         Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
         ki (imaginary part)
      */
      ierr = EPSGetEigenpair(eps,i,&kr,&ki,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);

      /*
         Compute the relative error associated to each eigenpair
      */
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
  
  /* 
     Free work space
  */
  ierr = EPSDestroy(eps);CHKERRQ(ierr);
  ierr = MatDestroy(A);CHKERRQ(ierr);
  ierr = MatDestroy(B);CHKERRQ(ierr);
  ierr = SlepcFinalize();CHKERRQ(ierr);
  return 0;

}