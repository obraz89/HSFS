#include "stdafx.h"
#include "conj_minmax.h"

#include "log.h"

static double ARG_TOL=1.0e-6;
static double GRAD_TOL=1.0e-7;
static double DARG_MIN=ARG_TOL;

//steep params
static int MAX_ITER_SD = 1000;

// conj params
static double GOLD = 1.618034;
// TODO: this should be important in stability part
static double GLIMIT = 10.0;
static double TINY = 1.0e-8;
static int MAX_ITER_CONJ=50;

using namespace smat;

//----------------------------------------------------------------t_GradSrchBase
t_GradSrchBase::t_GradSrchBase(int ndim):_ndim(ndim){}
//---------------------------------------------------------------t_SteepDescSrch
t_SteepDescSrch::t_SteepDescSrch(int ndim):t_GradSrchBase(ndim),
_arg_cur(ndim), _arg_nxt(ndim), _grad_cur(ndim){}

int t_SteepDescSrch::_search_min_desc(t_VecDbl& a_arg){

	_arg_cur = a_arg;

	bool converged = false;
	int n_iter=0;
	double fun_cur, fun_nxt;

	do 
	{

		fun_cur = _calc_fun_desc(_arg_cur);
		_grad_cur = _calc_grad_desc(_arg_cur);
		// IMPORTANT TODO: step?
		double dt = 0.005*_arg_cur.norm()/_grad_cur.norm();

		_arg_nxt = _arg_cur - dt*_grad_cur;
		fun_nxt = _calc_fun_desc(_arg_nxt);

		converged = fun_nxt>fun_cur;
		// TODO: max iterations?
		//if (n_iter++>MAX_ITER_CONJ) return 1;

		_arg_cur = _arg_nxt;

		// debug
		std::wstringstream wsstr;
		wsstr<<"Minimizing{SD} fun="<<std_manip::std_format_sci(fun_nxt);
		log_my::wxLogMessageStd(wsstr.str());

	} while (!converged);

	a_arg = _arg_cur;
	return 0;
}
//----------------------------------------------------------------t_ConjGradSrch
t_ConjGradSrch::t_ConjGradSrch(int ndim):t_GradSrchBase(ndim),
_arg_cur(ndim), _arg_nxt(ndim),
_h_cur(ndim), _h_nxt(ndim),
_g_cur(ndim), _g_nxt(ndim),
_lin_start(ndim), _lin_end(ndim){}

int smat::t_ConjGradSrch::_search_min_desc(t_VecDbl& init_guess){

	_arg_cur = init_guess;

	_g_cur = -1.0*_calc_grad_desc(init_guess);
	_h_cur = _g_cur;

	bool converged;
	int n_iter=0;
	do
	{
		_lin_start = _arg_cur;
		if (_lin_bracket(_lin_start, _lin_end, _h_cur)!=0){
			wxLogMessage(_T("ConjGrad: linmin bracketing error"));
		}
		double lam = _lin_min(_lin_start, _lin_end);
		_arg_nxt = _arg_cur + lam*_h_cur;

		_g_nxt = -1.0*_calc_grad_desc(_arg_nxt);

		double gamma = (vector::dot(_g_nxt - _g_cur,_g_nxt))/(vector::dot(_g_cur, _g_cur));
		// TODO: fast initial convergence with real conjugate grad,
		// but problems near the mimax
		// for now use steepest descent
		gamma=0.0;
		_h_nxt = _g_nxt + gamma*_h_cur;

		double cur_resid = (_arg_nxt - _arg_cur).norm()/(_arg_nxt+_arg_cur).norm() ;

		converged = cur_resid < sqrt(GRAD_TOL);

		_arg_cur = _arg_nxt;
		_h_cur = _h_nxt;
		_g_cur = _g_nxt;

		//  if no convergence
		if (n_iter++>MAX_ITER_CONJ) return 1;

		// debug
		std::wstringstream wsstr;
		wsstr<<"Minimizing{CG} fun="<<std_manip::std_format_sci(_calc_fun_desc(_arg_cur));
		log_my::wxLogMessageStd(wsstr.str());
	} while (!converged);

	std::wstringstream sstr;
	sstr<<_T("ConjGrad: converged niter=")<<n_iter;
	wxLogMessage(&(sstr.str()[0]));

	// lazy to rename init_guess
	init_guess = _arg_cur;

	return 0;

}

int t_ConjGradSrch::_lin_bracket(t_VecDbl& pnt_start, t_VecDbl& pnt_end, const t_VecDbl& dir){

	double cur_norm = pnt_start.norm();
	double lam = DARG_MIN*cur_norm;
	double f_cur, f_nxt;

	double f_start = _calc_fun_desc(pnt_start);
	double f_end;

	int n_iter = 0;

	// init pnt_end
	pnt_end = pnt_start;
	f_end = f_start;

	t_VecDbl dir_n = dir;
	dir_n.normalize();

	do 
	{
		pnt_start = pnt_end;
		f_start = f_end;

		pnt_end = pnt_start + lam*dir_n;
		f_end = _calc_fun_desc(pnt_end);
		lam*=2.0;

		// can't bracket the min
		if (n_iter++>MAX_ITER_CONJ) return 1;
	} while (f_start>=f_end);

	return 0;

}

//========================

int sign(double val){return val>=0.0 ? 1 : -1;}

int t_ConjGradSrch::_lin_bracket_brent(t_VecDbl& a_vec, t_VecDbl& c_vec, const t_VecDbl& a_dir){

	const int ndim = a_vec.size();
	double ax,bx,cx,ulim,u,r,q,fu,dum;
	double fa,fb,fc, f_tmp;
	t_VecDbl dum_vec(ndim);

	const t_VecDbl ref_point = a_vec;
	t_VecDbl dir = a_dir;
	dir.normalize();

	ax=0.0;
	// TODO: empirical 1% initial variation
	cx = std::max(0.01*a_vec.norm(), sqrt(TINY));
	bx = 0.5*(ax+cx);

	fa=_calc_fun_desc(ref_point + ax*dir);
	fb=_calc_fun_desc(ref_point + bx*dir);

#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((a)*sign(b))

	if (fb > fa) {
		SHFT(dum,ax,bx,dum)
		SHFT(dum,fb,fa,dum)
	}
	cx=bx+GOLD*(bx-ax);
	fc=_calc_fun_desc(ref_point + cx*dir);

	int n_iter=0;
	while ((fb > fc)&&(n_iter++<MAX_ITER_CONJ)) {
		//Compute u by parabolic extrapolation
		r=(bx-ax)*(fb - fc);
		q=(bx-cx)*(fb - fa);
		u=(bx)-((bx-cx)*q-(bx-ax)*r)/
			(2.0*SIGN(std::max(abs(q-r),TINY),q-r));
		ulim=bx+GLIMIT*(cx-bx);

		if ((bx-u)*(u-cx) > 0.0) {	//Parabolic u is between b and c: try it
				fu=_calc_fun_desc(ref_point + u*dir);
			if (fu < fc) {			//Got a minimum between b and c
				ax=bx;
				a_vec = ref_point + ax*dir;
				c_vec = ref_point + cx*dir;
				return 0;
			} else if (fu > fb) {	//Got a minimum between between a and u
					cx=u;
					fc=fu;
					a_vec = ref_point + ax*dir;
					c_vec = ref_point + cx*dir;
					return 0;
			}
			u=(cx)+GOLD*(cx-bx);	//Parabolic fit was no use. Use default magnification
			fu=_calc_fun_desc(ref_point + u*dir);
		} else if ((cx-u)*(u-ulim) > 0.0) {	//Parabolic fit is between c and its allowed limit
				fu=_calc_fun_desc(ref_point + u*dir);
				if (fu < fc) {
					SHFT(bx,cx,u,cx+GOLD*(cx-bx))
					SHFT(fb,fc,fu,_calc_fun_desc(ref_point + u*dir));
				}
		} else if ((u-ulim)*(ulim-cx) >= 0.0) { //Limit parabolic u to maximum allowed value
				u=ulim;
				fu=_calc_fun_desc(ref_point + u*dir);
			} else {							//Reject parabolic, use deafult magnification
					u=(cx)+GOLD*(cx-bx);
					fu=_calc_fun_desc(ref_point + u*dir);
		}
		SHFT(ax,bx,cx,u)
		SHFT(fa,fb,fc,fu)
	}	// main loop ends
#undef SHFT
#undef SIGN
	if (n_iter<MAX_ITER_CONJ){
		a_vec = ref_point + ax*dir;
		c_vec = ref_point + cx*dir;
		return 0;
	}else{
		wxLogMessage(_T("LinMin bracketing failure\n"));
		return 1;
	}
}
//========================

double t_ConjGradSrch::_lin_min(const t_VecDbl& pnt_start, const t_VecDbl& pnt_end){

	t_VecDbl lft = pnt_start;
	t_VecDbl rgt = pnt_end;
	t_VecDbl mid = 0.5*(lft+rgt);
	t_VecDbl x(mid);

	double tol;
	do 
	{
		double lft_len = (mid-lft).norm();
		double rgt_len = (rgt-mid).norm();

		if (lft_len>rgt_len){
			// x between left and mid
			x = 0.5*(lft+mid);
			double f_a = _calc_fun_desc(lft);
			double f_x = _calc_fun_desc(x);

			if (f_x<f_a){
				lft = x;
			}else{
				rgt = mid;
				mid = x;
			}
		}else{
			x = 0.5*(mid+rgt);

			double f_b = _calc_fun_desc(mid);
			double f_x = _calc_fun_desc(rgt);

			if(f_b<f_x){
				rgt = x;
			}else{
				lft=mid;
				mid=x;
			}

		}
		tol = (mid-lft).norm();

	} while (tol>ARG_TOL);

	return (mid-pnt_start).norm()/_h_cur.norm();

}


// =====================+DEBUG

smat::t_Conj2D::t_Conj2D(t_pFun2D a_fun, t_pFunGrad2D a_fun_grad):t_ConjGradSrch(2), 
_pFun(a_fun), _pFunGrad(a_fun_grad){}

double t_Conj2D::_calc_fun(const t_VecDbl& arg){
	return _pFun(arg[0], arg[1]);
}

t_VecDbl t_Conj2D::_calc_grad(const t_VecDbl& arg){
	return _pFunGrad(arg[0], arg[1]);
}


class t_Conj2D: public t_ConjGradSrch{
	t_pFun2D _pFun, _pFunGrad;
	double _calc_fun(const t_VecDbl& arg);
	t_VecDbl _calc_grad(const t_VecDbl& arg);
public:
	t_Conj2D(t_pFun2D fun, t_pFun2D fun_grad);
};
// =====================~DEBUG
