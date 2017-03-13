#include "stdafx.h"

#include "fun_interpolate.h"
#include "mkl.h"

//#include "errcheck.inc"
//#include "generatedata.inc"
//#include "rescheck.inc"

// in: argument x array of size nx {increasing}
// inout: function y array {when smoothed, data is modified}
// in: nx - sizeof x,y,y1,y2
// out: 1-st and 2-nd derivatives arrays y1 and y2

#define NDORDER        3 // size of array describing derivative orders
                         // to compute
#define NDER           2 // number of derivatives to compute

// use natural cubic interpolation
void smat::interpolate_profile_sm_deriv_cubic(double* a_x, double *a_y, const int a_nx, double* a_y1, double* a_y2){

	DFTaskPtr task;                     // Data Fitting task descriptor
    MKL_INT sorder;                     // spline order
    MKL_INT stype;                      // spline type

    MKL_INT xhint;                      // additional info about break points
    const MKL_INT ny = 1;                         // number of functions
    MKL_INT yhint;                      // additional info about function
    MKL_INT nscoeff = (a_nx-1)*DF_PP_CUBIC;				   // number of spline coefficients
    MKL_INT scoeffhint;                 // additional info about spline
                                        // coefficients
    MKL_INT bc_type;                    // boundary conditions type
    MKL_INT nbc;                        // number of boundary conditions
    MKL_INT ic_type;                    // internal conditions type
    MKL_INT nsite_bl;                   // number of interpolaton sites
                                        // in one block
    MKL_INT nblocks = 1;                    // number of blocks

    MKL_INT sitehint;                   // additional info about interpolation
                                        // sites
    MKL_INT ndorder;                    // size of array describing derivative
                                        // orders
    MKL_INT dorder[] = { 0, 1, 1 };     // 1st and 2nd derivatives
                                        // will be computed
    MKL_INT rhint;                      // interpolation results storage format
    MKL_INT *cell_idx;                  // indices of cells containing
                                        // interpolation sites
    double *ic;                         // array of internal conditions
    double *bc;                         // array of boundary conditions
    double* scoeff = new double[nscoeff];             // array of spline coefficients

	double* r = new double[NDER*a_nx];	// spline evaluation results
    double *datahint;                   // additional info about structure
                                        // of arrays x and y

    int i, j, errcode = 0;

    int errnums = 0;

    /***** Initializing parameters for Data Fitting task *****/

    sorder     = DF_PP_CUBIC;
    stype      = DF_PP_NATURAL;

    /***** Parameters describing interpolation interval *****/
    xhint      = DF_NO_HINT;

    /***** Parameters describing function *****/
    yhint      = DF_NO_HINT;

    /***** Parameters describing spline coefficients storage *****/
    scoeffhint = DF_NO_HINT;

    /***** Parameters describing boundary conditions type *****/
	bc_type = DF_BC_NOT_A_KNOT;

    /* No additional parameters are provided for boundary conditions */
    bc         = 0;

    /***** Parameters describing internal conditions type *****/
    /* No internal conditions are provided for natural cubic spline */
    ic_type    = DF_NO_IC;
    ic         = 0;

    /***** Parameters describing interpolation sites *****/
    sitehint   = DF_NON_UNIFORM_PARTITION;

    /***** Additional info about structure of arrays x and y *****/
    /* No additional info is provided */
    datahint   = 0;

    /**** Parameter describing interpolation results storage *****/
    rhint      = DF_MATRIX_STORAGE_ROWS;

    /**** Parameter describing array of cell indices *****/
    /* No cell indices are computed */
    cell_idx   = 0;

    /**** Parameter describing size of array for derivative orders *****/
    ndorder    = NDORDER;

    /***** Create Data Fitting task *****/
    errcode = dfdNewTask1D( &task, a_nx, a_x, xhint, ny, a_y, yhint );

    /***** Edit task parameters for natural cubic spline construction *****/

    errcode = dfdEditPPSpline1D( task, sorder, stype, bc_type, 0, ic_type, 0,
                                 scoeff, scoeffhint );
    //CheckDfError(errcode);

    /***** Construct natural cubic spline using STD method *****/
    errcode = dfdConstruct1D( task, DF_PP_SPLINE, DF_METHOD_STD );
    //CheckDfError(errcode);

    /***** Interpolate, and compute derivatives using STD method *****/

    errcode = dfdInterpolate1D( task, DF_INTERP, DF_METHOD_PP,
                                    a_nx, a_x, sitehint, ndorder,
                                    dorder, datahint, r, rhint, cell_idx );


    for ( i = 0; i < a_nx; i++ )
    {
		a_y1[i] = r[i*NDER+0];
		a_y2[i] = r[i*NDER+1];

    }

	errcode = dfDeleteTask( &task );
	delete[] scoeff, r;

    return;

}