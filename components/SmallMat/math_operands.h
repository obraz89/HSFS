///////////////////////////////////////////////////////////////////////////////
// Project:	SmallMat
// Purpose:	Small interface
///////////////////////////////////////////////////////////////////////////////
// File:        math_operands.h
// Purpose:     all exports here, to be used in other projects
// Author:      A.Obraz
///////////////////////////////////////////////////////////////////////////////


#ifndef __SMALL_MAT_OPERANDS
#define __SMALL_MAT_OPERANDS

// matrix operands
#include "small_mat.h"

// try to instantiate needed classes

template class IMPEXP_SMALLMAT t_Matrix<t_Complex>;
template class IMPEXP_SMALLMAT t_Matrix<double>;

template class IMPEXP_SMALLMAT t_Vec<t_Complex>;
template class IMPEXP_SMALLMAT t_Vec<double>;

template class IMPEXP_SMALLMAT t_Vec3<t_Complex>;
template class IMPEXP_SMALLMAT t_Vec3<double>;

template class IMPEXP_SMALLMAT t_SqMatrix<t_Complex>;
template class IMPEXP_SMALLMAT t_SqMatrix<double>;

template class IMPEXP_SMALLMAT t_SqMat3<t_Complex>;
template class IMPEXP_SMALLMAT t_SqMat3<double>;


// gradient search methods
#include "conj_minmax.h"

// finding roots F(x)=0
#include "fun_zero_1D.h"

// interpolation of functions
#include "fun_interpolate.h"

// integration of functions
#include "fun_integrate.h"

#endif // __SMALL_MAT_OPERANDS