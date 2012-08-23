#include "math_operands.h"
class t_ODES {
protected:
	struct t_OrthPoint{
	t_OrthPoint(const int& a_ind, const int& a_dim);
	t_SqMatCmplx orthMatrix;
	int ind;
};
	std::vector<t_OrthPoint> _orthStack;
protected:
	// _dim - for 2D dim=3 for 3D dim=4
	// _nnodes - number of grid nodes interval is splitted into
	// _pFunRHS - reference to the right hand side function 
	//		(for stability rhs = A*f, f - one of solution vectors)
	// _orthMatrix - matrix of transition from initial basis to current 
	//		(after orthogonalization)
	static const int _dim;
	int _nnodes;
	// methods
	// stepRK - Runge-Kutta step for a particular vector, returns solution vector
	// needOrtho - check if orthogonalization is needed for current basis
	// ortho - makes transition to new basis of solution vectors
	// reconstruct - reconstructs the eigenfunctions

	// solution[i] - vectors of solution at i point of grid  0 < i < _nnodes - 1
	// solution[i][j] - particular vector 0 < j < _ndim - 1
	// solution[i][j][k] - k-component of eigenvector (raw)

	//t_Vec stepRK2D(const int& ind, const t_Vec& fun); 
	t_VecCmplx stepRK3D(const int& ind, const t_VecCmplx& fun); 
	virtual bool needOrtho(const t_MatCmplx& vecs)=0; 
	// Gramm - Shmidt determinant of rank k, basis v0, v1, ...
	// Ermith matrix
	//	( (v0, v0)		(v1, v0)*		 ... (v[k-1], v0)*		)
	//	(  ...													)
	//	(  (v[k-1], v0) (v[k-1], v1)	 ... (v[k-1], v[k-1])*	)

	t_Complex detGS(const t_MatCmplx& sol, const int& rank) const;
	t_Complex minorGS(const t_MatCmplx& sol, const int& dim, const int& nExcludeCol) const; // rank = dim+1 
	void ortho(const int& nnode);	

	//virtual t_Vec formRHS2D(const double& a_y, const t_Vec& a_var) = 0;
	virtual t_VecCmplx formRHS3D(const double& a_y, const t_VecCmplx& a_var) const = 0; 	
	virtual t_Complex getResidual3D() = 0;
	//virtual t_Complex getResidual2D() = 0;
public:
	// members
	std::vector<double> varRange;
	// container for solution vectors - raw, without reconstruction
	// reconstruction is t-consuming - to obtain reconstructed solution
	// reconstruct should be invoked - this is only needed when 
	// eigenfunctions are of interest (not eigenvalues)
	std::vector<t_MatCmplx> solution; 

	t_ODES();
	//getters
	int getTaskDim() const{
		return _dim;
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
