#include "stdafx.h"

#include "lin_sys.h"

#include "mkl.h"

inline MKL_Complex16 getMKLCmplx(double ar, double ai){
	MKL_Complex16 r;
	r.real = ar; r.imag = ai;
	return r;}

inline MKL_Complex16 getMKLCmplx(const t_CompVal& val){
	MKL_Complex16 r;
	r.real = val.real();
	r.imag = val.imag();
	return r;
}

t_Complex MKLCmplx2std(const MKL_Complex16& mkl_cmplx){
	return t_Complex(mkl_cmplx.real, mkl_cmplx.imag);

}

// in gzesv we have column-based storage scheme
inline int calcPlainIndByIJ(int col, int row, int ndim){
	return row+col*ndim;
}

void smat::solve_lsys_lu(const t_SqMatCmplx& A, const t_VecCmplx& b, t_VecCmplx& x){

#ifdef _DEBUG
	if (A.nCols()!=b.nRows()) 
		ssuTHROW(matrix::t_Not3D, 
		_T("smat: lin sys LU error: mat-vec size mismatch"));
#endif

	const int N = A.nCols();

	const int nrhs = 1;

	MKL_Complex16 * a_pl = new MKL_Complex16[N*N];
	MKL_Complex16 * b_pl = new MKL_Complex16[N];

	int*  ipiv = new int[N];

	int ind;

	for (int i=0; i<N; i++){

		b_pl[i] = getMKLCmplx(b[i]);

		// TODO: matrix storage scheme ?!
		for (int j=0; j<N; j++){

			ind = calcPlainIndByIJ(i, j, N);
			a_pl[ind] = getMKLCmplx(A[i][j]);

		} 
	}

	int info;

	zgesv(&N, &nrhs, a_pl, &N, ipiv, b_pl, &N, &info);

	for (int i=0; i<N; i++){
		x[i] = MKLCmplx2std(b_pl[i]);
	}


}