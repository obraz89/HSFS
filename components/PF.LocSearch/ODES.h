#include "math_operands.h"

class t_ODES {
protected:

	struct t_OrthPoint{
		t_OrthPoint(const int& a_ind, const int& a_dim);
		t_SqMatCmplx orthMatrix;
		int ind;
	};

	std::vector<t_OrthPoint> _orthStack;
	int _orth_stack_len;

	// cashed vecs and matrices - for performance boost
	// TODO: sizes and initialization!!!
	t_VecCmplx _rk1, _rk2, _rk3, _rk4, _hl, _hr, _h_res, 
			   _sol_vec_tmp1,_sol_vec_tmp2, 
			   _trans_mat_vec_tmp1, _trans_mat_vec_tmp2;

	t_MatCmplx _sol_tmp;
	// _NSolVecs - number of solution vectors to solve
	// _orthMatrix - matrix of transition from initial basis to current 
	//		(after orthogonalization)
	const int _NSolVecs;
	// _SolVecDim - dimension of any solution vectors
	const int _SolVecDim;
	// _nnodes - number of grid nodes interval is splitted into
	int _nnodes;
	// methods

	void stepRK3D(const int& ind, const t_VecCmplx& fun, t_VecCmplx& dest); 
	// needOrtho - check if orthogonalization is needed for current basis
	virtual bool needOrtho(const t_MatCmplx& vecs)=0; 

	t_Complex detGS(const t_MatCmplx& sol, const int& rank) const;
	t_Complex minorGS(const t_MatCmplx& sol, const int& dim, const int& nExcludeCol) const; // rank = dim+1 
	void ortho(const int& nnode);	

	virtual void formRHS3D(const double& a_y, const t_VecCmplx& a_var, t_VecCmplx& rhs) = 0; 	
public:
	// members
	std::vector<double> varRange;

	// solution[i] - vectors of solution at i point of grid  0 < i < _nnodes - 1
	// solution[i][j] - particular vector 0 < j < _ndim - 1
	// solution[i][j][k] - k-component of eigenvector (raw)
	std::vector<t_MatCmplx> solution; 

	t_ODES(const int NSolVecs, const int SolVecDim);
	//getters
	int getTaskDim() const{
		return _NSolVecs;
    }

	int getNNodes() const{
		return _nnodes;
	} 
	virtual void clean();
	virtual void setContext(const int& a_newNnodes);
	virtual void solve();
	std::vector<t_MatCmplx> reconstruct();
	virtual ~t_ODES(){};	
};

//=============================debug
/*
class t_ODESTest : public t_ODES{
public:
	typedef t_SqMatCmplx (*pFunRHS)(const double);
private:
	pFunRHS _rhs;
	t_VecCmplx formRHS3D(const double& a_y, const t_VecCmplx& a_var) const{
		return _rhs(a_y)*a_var;
	}
	t_Complex getResidual3D(){return 0;}
	bool needOrtho(const t_MatCmplx& a_cur_sol){
		//debug
		double max_norm = 0.0;
		double min_norm = 1.0e+12;
		for (int i=0; i<a_cur_sol.nCols(); i++){
			double cur_norm = complex::norm(t_VecCmplx(a_cur_sol[i]).norm());
			if (cur_norm<min_norm){
				min_norm = cur_norm;
			}
			if (cur_norm>max_norm){
				max_norm = cur_norm;
			}
		}

		// TODO: tolerance should be param
		// for alg. solver
		return max_norm>1000.0;
	}
public:
	void init(double y_min, double y_max, int nnodes, t_MatCmplx ddd){
		setContext(nnodes);
		double dy = (y_max-y_min)/double(nnodes-1);
		for (int i=0; i<nnodes; i++){
			varRange[i] = y_min + dy*i;
		}
		solution[0] = ddd;
	}

	t_ODESTest(pFunRHS rhs):_rhs(rhs), t_ODES(){};
};
//
*/ // t_ODESTest