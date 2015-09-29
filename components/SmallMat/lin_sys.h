#ifndef __LIN_SYS
#define __LIN_SYS

#include "dll_impexp_smat.h"

#include "small_mat.h"

namespace smat{

	IMPEXP_SMALLMAT void solve_lsys_lu(const t_SqMatCmplx& A, const t_VecCmplx& b, t_VecCmplx& x);

}

#endif  // #ifndef __LIN_SYS