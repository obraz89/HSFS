#include "ODES_operands.h"
class t_ODES {
private:
	// _dim - for 2D dim=3 for 3D dim=4
	// _nnodes - number of grid nodes interval is splitted into
	// _pFunRHS - reference to the right hand side function 
	//		(for stability rhs = A*f, f - one of solution vectors)
	// _orthMatrix - matrix of transition from initial basis to current 
	//		(after orthogonalization)
	const int _dim;
	int _nnodes;
	struct t_OrthPoint{
		t_OrthPoint(const int& a_ind, const int& a_dim);
		t_SqMatrix orthMatrix;
		int ind;
	};
	std::vector<t_OrthPoint> _orthStack;
	// methods
	// stepRK - Runge-Kutta step for a particular vector, returns solution vector
	// needOrtho - check if orthogonalization is needed for current basis
	// ortho - makes transition to new basis of solution vectors
	// reconstruct - reconstructs the eigenfunctions

	// solution[i] - vectors of solution at i point of grid  0 < i < _nnodes - 1
	// solution[i][j] - particular vector 0 < j < _ndim - 1
	// solution[i][j][k] - k-component of eigenvector (raw)
	t_Vec stepRK(const int& ind, const t_Vec& fun); 
	bool needOrtho(const t_Matrix& vecs); 
	// Gramm - Shmidt determinant of rank k, basis v0, v1, ...
	// Ermith matrix
	//	( (v0, v0)		(v1, v0)*		 ... (v[k-1], v0)*		)
	//	(  ...													)
	//	(  (v[k-1], v0) (v[k-1], v1)	 ... (v[k-1], v[k-1])*	)

	t_Complex detGS(const t_Matrix& sol, const int& rank) const;
	t_Complex minorGS(const t_Matrix& sol, const int& dim, const int& nExcludeCol) const; // rank = dim+1 
	void ortho(const int& nnode);	

	virtual t_Vec formRHS(const double& a_y, const t_Vec& a_var) = 0;
	virtual void setInitials() = 0; 	
	virtual t_Complex getResidual3D() = 0;
	virtual t_Complex getResidual2D() = 0;
public:
	// members
	std::vector<double> varRange;
	// container for solution vectors - raw, without reconstruction
	// reconstruction is t-consuming - to obtain reconstructed solution
	// reconstruct should be invoked - this is only needed when 
	// eigenfunctions are of interest (not eigenvalues)
	std::vector<t_Matrix> solution; 

	t_ODES(const int& a_dim, const int& a_nnodes);
	// to change grid size
	void resizeGrid(const int& a_newNnodes);
	virtual void solve();
	std::vector<t_Matrix> reconstruct();
	virtual ~t_ODES(){};	
};
