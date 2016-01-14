#include "stdafx.h"

#include "StabSolver.h"
#include "log.h"

using namespace pf;

// set direct matrix of stability eqs H:
// dphi/dy = H*phi
void t_StabSolver::_setStabMatrix3D(const t_ProfRec& rec){
	t_CompVal imagUnity(0.0, 1.0);

	const double& stabRe = _profStab.scales().ReStab;
	const double& Me = _profStab.scales().Me;

	const mf::t_FldParams& Params = _rFldNS.get_mf_params();

	const double gMaMa = Params.Gamma*Me*Me;
	const double g_1MaMa = (Params.Gamma-1.0)*Me*Me;
	const double Pr = Params.Pr;

	const double inv_t = 1.0/rec.t;
	const double inv_mu = 1.0/rec.mu;

	const double vCoefL = 2.0/3.0*(Params.BulkViscRatio+2.0); 
	const double vCoefS = 2.0/3.0*(Params.BulkViscRatio-1.0); 
	const double vCoefM = 1.0+vCoefS;

	const t_CompVal& alpha = _waveChars.a;
	const t_CompVal& beta = _waveChars.b;
	const t_CompVal& freq = _waveChars.w;

	const t_CompVal dEicon_dt = alpha*rec.u + beta*rec.w - freq;

	const double dMu1 = rec.mu1*rec.t1;
	const double dMu2 = rec.mu2*rec.t1*rec.t1 + rec.mu1*rec.t2;
	const double dMu_U= rec.mu2*rec.t1*rec.u1 + rec.mu1*rec.u2;

	const t_CompVal xi = 1.0/(stabRe*inv_mu+imagUnity*vCoefL*gMaMa*dEicon_dt);

	// compose stab matrix for direct problem first

	// first
	_stab_matrix[1][0]=1.0;
	// second
	_stab_matrix[0][1] = imagUnity*stabRe*inv_t*inv_mu*dEicon_dt +
		pow(alpha,2) + pow(beta,2);

	_stab_matrix[1][1] =  - inv_mu*dMu1;

	_stab_matrix[2][1] = stabRe*inv_t*inv_mu*rec.u1 - 
		imagUnity*alpha*(inv_mu*dMu1+vCoefM*inv_t*rec.t1);

	_stab_matrix[3][1] = alpha*(imagUnity*stabRe*inv_mu -
		vCoefM*gMaMa*dEicon_dt);

	_stab_matrix[4][1] = vCoefM*alpha*inv_t*dEicon_dt - inv_mu*dMu_U;

	_stab_matrix[5][1] = -inv_mu*rec.mu1*rec.u1;
	// third
	_stab_matrix[0][2] = -imagUnity*alpha;

	_stab_matrix[2][2] = inv_t*rec.t1;

	_stab_matrix[3][2] = -imagUnity*gMaMa*dEicon_dt;

	_stab_matrix[4][2] = imagUnity*inv_t*dEicon_dt;

	_stab_matrix[6][2] = -imagUnity*beta;
	// fourth
	_stab_matrix[0][3] = -imagUnity*xi*alpha*
		(2.0*inv_mu*dMu1 + vCoefL*inv_t*rec.t1);

	_stab_matrix[1][3] = -imagUnity*xi*alpha;

	_stab_matrix[2][3] = xi*(-alpha*alpha-beta*beta+
		vCoefL*inv_t*inv_mu*dMu1*rec.t1+
		vCoefL*inv_t*rec.t2-
		imagUnity*stabRe*inv_t*inv_mu*dEicon_dt);

	_stab_matrix[3][3] = -imagUnity*xi*vCoefL*gMaMa*
		(
		(inv_mu*dMu1+inv_t*rec.t1)*dEicon_dt+
		alpha*rec.u1+beta*rec.w1
		);

	_stab_matrix[4][3] = imagUnity*xi*
		(
		(inv_mu*rec.mu1+vCoefL*inv_t)*(alpha*rec.u1+beta*rec.w1)+
		vCoefL*inv_t*inv_mu*dMu1*dEicon_dt
		);

	_stab_matrix[5][3] = imagUnity*xi*vCoefL*inv_t*dEicon_dt;

	_stab_matrix[6][3] = -imagUnity*xi*beta*
		(2.0*inv_mu*dMu1 + vCoefL*inv_t*rec.t1);

	_stab_matrix[7][3] = -imagUnity*xi*beta;
	// fifth
	_stab_matrix[5][4]=1.0;
	// sixth
	_stab_matrix[1][5] = -2.0*Pr*g_1MaMa*rec.u1;

	_stab_matrix[2][5] = Pr*
		(
		stabRe*inv_t*inv_mu*rec.t1 - 
		2.0*imagUnity*g_1MaMa*(alpha*rec.u1+beta*rec.w1)
		);

	_stab_matrix[3][5] = -imagUnity*stabRe*Pr*inv_mu*g_1MaMa*dEicon_dt;

	_stab_matrix[4][5] = imagUnity*stabRe*Pr*inv_t*inv_mu*dEicon_dt +
		alpha*alpha + beta*beta - 
		g_1MaMa*Pr*inv_mu*rec.mu1*(rec.u1*rec.u1+rec.w1*rec.w1)-
		inv_mu*dMu2;

	_stab_matrix[5][5] = -2.0*inv_mu*dMu1;

	_stab_matrix[7][5] = -2.0*Pr*g_1MaMa*rec.w1;
	// seventh
	_stab_matrix[7][6]=1.0;
	// last
	_stab_matrix[2][7] = -imagUnity*beta*(inv_mu*dMu1+vCoefM*inv_t*rec.t1)+
		stabRe*inv_mu*inv_t*rec.w1;

	_stab_matrix[3][7] = imagUnity*stabRe*beta*inv_mu-
		vCoefM*beta*gMaMa*dEicon_dt;

	_stab_matrix[4][7] = vCoefM*beta*inv_t*dEicon_dt-
		inv_mu*(rec.mu2*rec.t1*rec.w1+rec.mu1*rec.w2);

	_stab_matrix[5][7] = -inv_mu*rec.mu1*rec.w1;

	_stab_matrix[6][7] = imagUnity*stabRe*inv_t*inv_mu*dEicon_dt +
		alpha*alpha + beta*beta;

	_stab_matrix[7][7] = -inv_mu*dMu1;

}

void t_StabSolver::_setStabMatrix3D(const double& a_y){

	const t_ProfRec& rec = _profStab.get_rec(a_y);
	_setStabMatrix3D(rec);
};

// Matrix Hx for scalar product of 2 amplitude functions
// 
// Stability matrix corresponds to Nayfeh A.H., Stability of Three-Dimensional
// Boundary Layers," AIAA J., Vol. 18, No. 4, pp. 406-416, 1980.
// In this matrix:
//     m=MM=(e-1)2/3
//     r=RM=(e+2)2/3
//     e=K
//     e=(3/2)*(bulk viscosity)/(viscosity)

// Mack recommends (bulk viscosity)/(viscosity)=0.8, then parameter K=1.2
// TODO: cases where bv/v=0.8 ?
//
//
//	N = number of points in Y-grid of constant step (N should be smaller than parameter NY) 
//	Y1 = upper boundary (real)
//	Basic flow profiles in COMMON/BASIC1/
//	Basic flow parameters in COMMON/BASIC2/
//	DZ(8,N) = vector of conjugate stability problem (complex) (in COMMON/SCAL1/)
//	A,B,W = wavenumbers and frequency of stability problem (complex) (in COMMON/SCAL2/)
//	R = Reynolds number (complex) (in COMMON/SCAL2/)
//	A0(8,N) vector of disturbance (complex) (in COMMON/FIELD/)

void t_StabSolver::_setScalProdMatrix(const t_ProfRec& rec){
	//void set_stab_h(int i,const t_AmpVec<double>& vec_mf, t_SqMatCmplx& h_mat){

		t_CompVal E(0.0,1.0);

		const mf::t_FldParams& Params = _rFldNS.get_mf_params();

		double MM = 2./3.*(Params.BulkViscRatio-1.0);

		double F = MM+1.0;

		double RM = 2./3.*(Params.BulkViscRatio+2.0);

		const double Gamma = Params.Gamma;
		const double Mach = Params.Mach;

		const double Pr = Params.Pr;

		double M2 = Gamma*Mach*Mach;

		double MG = (Gamma-1.0)*Mach*Mach;

		_scal_prod_matrix.setToZero();

		const t_CompVal W = _waveChars.w;
		const t_CompVal A = _waveChars.a;
		const t_CompVal B = _waveChars.b;
		const t_CompVal R = _profStab.scales().ReStab;

		const double U = rec.u;
		const double U1 = rec.u1;
		const double U2 = rec.u2;

		const double T = rec.t;
		const double T1 = rec.t1;
		const double T2 = rec.t2;

		const double MU = rec.mu;
		const double MU1 = rec.mu1;
		const double MU2 = rec.mu2;

		const double WS = rec.w;
		const double WS1 = rec.w1;
		const double WS2 = rec.w2;

		const t_CompVal WA = W - A*U - B*WS;
		const double M1 = 1.0/MU;
		const double MY = MU1*T1;
		const double MYY = MU2*T1*T1+MU1*T2;
		const t_CompVal XI = 1.0/(R*M1-E*RM*M2*WA);
		const t_CompVal TD = 1.0/T;
		const t_CompVal WB = E*(A*U-W);

		const t_CompVal DXI = -E*RM*M2*U*XI*XI;

		// fill non-zero elements

		_scal_prod_matrix[0][1] = U*R*M1*TD-2.*E*A;

		_scal_prod_matrix[2][1]=-F*T1*TD-MY*M1;

		_scal_prod_matrix[3][1]=R*M1-E*F*M2*(WA-A*U);

		_scal_prod_matrix[4][1]=E*F*TD*(WA-A*U);

		_scal_prod_matrix[0][2]=-1.;

		_scal_prod_matrix[3][2]=-M2*U;

		_scal_prod_matrix[4][2]=U*TD;

		_scal_prod_matrix[0][3]=-(XI+A*DXI)*(RM*T1*TD+2.*MY*M1);

		_scal_prod_matrix[1][3]=-XI-A*DXI;

		_scal_prod_matrix[2][3]=-E*(
			DXI*(E*WA*R*M1*TD-A*A-B*B+RM*(T2*TD+MY*T1*M1*TD))
			-XI*(2.*A+E*R*U*M1*TD));

		_scal_prod_matrix[3][3]=-RM*M2*
			(DXI*(A*U1+B*WS1)+XI*U1+(T1*TD+MY*M1)*(XI*U-DXI*WA));

		_scal_prod_matrix[4][3]=DXI*((A*U1+B*WS1)*(RM*TD+MU1*M1)-RM*WA*MY*M1*TD)
			+XI*(RM*U1*TD+MU1*U1*M1+RM*U*MY*M1*TD);

		_scal_prod_matrix[5][3]=RM*TD*(U*XI-WA*DXI);

		_scal_prod_matrix[6][3]=-DXI*B*(RM*T1*TD+2.*MY*M1);

		_scal_prod_matrix[7][3]=-B*DXI;

		_scal_prod_matrix[2][5]=-2.*MG*Pr*U1;

		_scal_prod_matrix[3][5]=-MG*Pr*R*U*M1;

		_scal_prod_matrix[4][5]=R*Pr*U*M1*TD-2.*E*A;

		_scal_prod_matrix[3][7]=E*F*M2*B*U;

		_scal_prod_matrix[4][7]=-E*F*B*U*TD;

		_scal_prod_matrix[6][7]=R*U*M1*TD-2.*E*A;

}

void t_StabSolver::_setScalProdMatrix(const double& a_y){

	const t_ProfRec& rec = _profStab.get_rec(a_y);
	_setScalProdMatrix(rec);

}

// TODO: to speed up - rewrite to "_setAsymptotcs" 
// and write directly to destination formal param

void t_StabSolver::_setAsymptotics(t_MatCmplx& asym_vecs){

	const int dim = getTaskDim();

	t_SqMatCmplx b_coef(dim);
	t_VecCmplx lambda(dim,0.0);

	// TODO: function for simplified asymp: u=1.0, u'=0, u''=0, ... ?
	t_ProfRec out_rec = _profStab.get_last_rec();

	if (_ls_mode.is_flag_on(stab::t_LSMode::ASYM_HOMOGEN)){
		// modify ... a little
		// this is how AVF did it
		out_rec.t1=0.0;
		out_rec.t2=0.0;
		out_rec.u1=0.0;
		out_rec.u2=0.0;
		out_rec.w1=0.0;
		out_rec.w2=0.0;
	}
	_setStabMatrix3D(out_rec);
	// to shorten 
	t_SqMatCmplx& _sm = _stab_matrix;

	// det(H - lambda*E)=0 <=> lambda1, lambda2... lambda8
	// choose l1,l3,l5,l7 : Re(l1,l3,l5,l7)<0

	// eigen vectors h1 for l1 eigenvalue, h3 for l3, h5 for l5, h7 for l7 are used as initials

	b_coef[0][0]=_stab_matrix[0][1];
	b_coef[1][0]=_stab_matrix[3][1];
	b_coef[2][0]=_stab_matrix[4][1];

	b_coef[1][1]=_sm[3][1]*_sm[1][3]+_sm[3][2]*_sm[2][3]+
		_sm[3][5]*_sm[5][3]+_sm[3][7]*_sm[7][3];

	b_coef[2][1]=_sm[4][1]*_sm[1][3]+_sm[4][2]*_sm[2][3]+
		_sm[5][3]*_sm[4][5]+_sm[7][3]*_sm[4][7];

	b_coef[1][2]=_sm[3][5];
	b_coef[2][2]=_sm[4][5];

	b_coef[1][3]=_sm[3][7];
	b_coef[2][3]=_sm[4][7];

	b_coef[3][3] = _stab_matrix[0][1];
	//

	t_CompVal s1 = 0.5*(b_coef[1][1]+b_coef[2][2]);
	t_CompVal s2 = sqrt(
		0.25*std::pow(b_coef[1][1]-b_coef[2][2],2)+
		b_coef[2][1]*b_coef[1][2]
	);
	lambda[0] = -sqrt(b_coef[0][0]);
	lambda[1] = -sqrt(s1+s2);
	lambda[2] = -sqrt(s1-s2);
	lambda[3] = lambda[0];

	b_coef[0][0] = 1.0;
	b_coef[1][0] = 0.0;
	b_coef[2][0] = 0.0;
	b_coef[3][0] = 0.0;
	for (int i=1; i<3; i++){
		t_CompVal L2 = std::pow(lambda[i],2);
		t_CompVal denom = _stab_matrix[0][1] - L2;
		b_coef[0][i] = (
			(L2 - _stab_matrix[4][5])*_stab_matrix[3][1]+
			_stab_matrix[4][1]*_stab_matrix[3][5]
		)/denom;

		b_coef[1][i] = _stab_matrix[4][5] - L2;
		b_coef[2][i] = -_stab_matrix[3][5];
		b_coef[3][i] = (
			_stab_matrix[3][5]*_stab_matrix[4][7]+
			(L2 - _stab_matrix[4][5])*_stab_matrix[3][7]
		)/denom;
	};

	b_coef[0][3]=0.0;
	b_coef[1][3]=0.0;
	b_coef[2][3]=0.0;
	b_coef[3][3]=1.0;

	for (int i=0; i<4; i++){	
		asym_vecs[i][0] = b_coef[0][i];
		asym_vecs[i][1] = lambda[i]*b_coef[0][i];
		asym_vecs[i][2] = (
			_stab_matrix[0][2]*b_coef[0][i]+
			_stab_matrix[3][2]*b_coef[1][i]+
			_stab_matrix[4][2]*b_coef[2][i]+
			_stab_matrix[6][2]*b_coef[3][i]
		)/lambda[i];

		asym_vecs[i][3] = b_coef[1][i];
		asym_vecs[i][4] = b_coef[2][i];
		asym_vecs[i][5] = lambda[i]*b_coef[2][i];
		asym_vecs[i][6] = b_coef[3][i];
		asym_vecs[i][7] = (
			_stab_matrix[3][7]*b_coef[1][i]+
			_stab_matrix[4][7]*b_coef[2][i]+
			_stab_matrix[6][7]*b_coef[3][i]
		)/lambda[i];
	}

	// verification of asymptotics, enable when debug needed
	
	t_VecCmplx asym_resid_v(2*dim, 0.0), init_vec(2*dim, 0.0);
	t_VecCmplx v1(2*dim, 0.0), v2(2*dim, 0.0), v3(2*dim, 0.0);
	double resid=0.0;

	for (int i=0; i<dim; i++){

			asym_vecs.col_to_vec(i, init_vec);
			//asym_resid_v = _stab_matrix*init_vec - lambda[i]*init_vec;		
			matrix::base::mat_mul<t_Complex, t_Complex>(_stab_matrix, init_vec, v1);
			matrix::base::mul<t_Complex, t_Complex>(lambda[i], init_vec, v2);
			matrix::base::minus<t_Complex, t_Complex>(v1, v2, asym_resid_v);

			resid = resid + asym_resid_v.norm().real();
	}

	//wxLogMessage(_T("Verify asymptotics v1: resid = %f"), resid);
	//if (resid>ASYM_TOL_DEFAULT)
	//	ssuGENTHROW(_T("StabSolver Error: Verification of Asymptotics failed"));
	
	// asym_vecs are set
}

void t_StabSolver::_setAsymptotics_v2(t_MatCmplx& asym_vecs){

	const int dim = getTaskDim();
	t_SqMatCmplx b_coef(dim);

	// matrix of fundamental solutions of direct problem
	// first 4 vecs are initials with Re[lambda]<0
	// others are with Re[lambda]>0 : lambda[i+4] = -lambda[i], i=0..3
	t_SqMatCmplx Z(STAB_MATRIX_DIM);

	t_VecCmplx lambda(STAB_MATRIX_DIM,0.0);
		
	// TODO: function for simplified asymp: u=1.0, u'=0, u''=0, ... ?
	t_ProfRec out_rec = _profStab.get_last_rec();

	if (_ls_mode.is_flag_on(stab::t_LSMode::ASYM_HOMOGEN)){
		// modify ... a little
		// this is how AVF did it
		out_rec.t1=0.0;
		out_rec.t2=0.0;
		out_rec.u1=0.0;
		out_rec.u2=0.0;
		out_rec.w1=0.0;
		out_rec.w2=0.0;
	}

	_setStabMatrix3D(out_rec);

	t_SqMatCmplx& _sm = _stab_matrix;


	// AVF method of computing Z
	{

		b_coef[0][0]=_stab_matrix[0][1];
		b_coef[1][0]=_stab_matrix[3][1];
		b_coef[2][0]=_stab_matrix[4][1];

		b_coef[1][1]=_sm[3][1]*_sm[1][3]+_sm[3][2]*_sm[2][3]+
			_sm[3][5]*_sm[5][3]+_sm[3][7]*_sm[7][3];

		b_coef[2][1]=_sm[4][1]*_sm[1][3]+_sm[4][2]*_sm[2][3]+
			_sm[5][3]*_sm[4][5]+_sm[7][3]*_sm[4][7];

		b_coef[1][2]=_sm[3][5];
		b_coef[2][2]=_sm[4][5];

		b_coef[1][3]=_sm[3][7];
		b_coef[2][3]=_sm[4][7];

		b_coef[3][3] = _stab_matrix[0][1];
		//

		t_CompVal s1 = 0.5*(b_coef[1][1]+b_coef[2][2]);
		t_CompVal s2 = sqrt(
			0.25*std::pow(b_coef[1][1]-b_coef[2][2],2)+
			b_coef[2][1]*b_coef[1][2]
		);
		lambda[0] = -sqrt(b_coef[0][0]);
		lambda[1] = -sqrt(s1+s2);
		lambda[2] = -sqrt(s1-s2);
		lambda[3] = lambda[0];

		for (int i=0; i<dim; i++) lambda[i+dim] = -lambda[i];

		b_coef[0][0] = 1.0;
		b_coef[1][0] = 0.0;
		b_coef[2][0] = 0.0;
		b_coef[3][0] = 0.0;
		for (int i=1; i<3; i++){
			t_CompVal L2 = std::pow(lambda[i],2);
			t_CompVal denom = _stab_matrix[0][1] - L2;
			b_coef[0][i] = (
				(L2 - _stab_matrix[4][5])*_stab_matrix[3][1]+
				_stab_matrix[4][1]*_stab_matrix[3][5]
			)/denom;

			b_coef[1][i] = _stab_matrix[4][5] - L2;
			b_coef[2][i] = -_stab_matrix[3][5];
			b_coef[3][i] = (
				_stab_matrix[3][5]*_stab_matrix[4][7]+
				(L2 - _stab_matrix[4][5])*_stab_matrix[3][7]
			)/denom;
		};

		b_coef[0][3]=0.0;
		b_coef[1][3]=0.0;
		b_coef[2][3]=0.0;
		b_coef[3][3]=1.0;

		for (int k=0; k<2; k++)
		for (int i=0; i<4; i++){	
			int ind = k*dim + i;
			Z[ind][0] = b_coef[0][i];
			Z[ind][1] = lambda[ind]*b_coef[0][i];
			Z[ind][2] = (
				_stab_matrix[0][2]*b_coef[0][i]+
				_stab_matrix[3][2]*b_coef[1][i]+
				_stab_matrix[4][2]*b_coef[2][i]+
				_stab_matrix[6][2]*b_coef[3][i]
			)/lambda[ind];

			Z[ind][3] = b_coef[1][i];
			Z[ind][4] = b_coef[2][i];
			Z[ind][5] = lambda[ind]*b_coef[2][i];
			Z[ind][6] = b_coef[3][i];
			Z[ind][7] = (
				_stab_matrix[3][7]*b_coef[1][i]+
				_stab_matrix[4][7]*b_coef[2][i]+
				_stab_matrix[6][7]*b_coef[3][i]
			)/lambda[ind];
		}
	
	} // Z calcs by AVF done

	// compute Z following Forgoston diss
	/*
	{

	b_coef[0][0]=_stab_matrix[0][1];

	b_coef[1][1]=_sm[3][1]*_sm[1][3]+_sm[3][2]*_sm[2][3]+
			_sm[3][5]*_sm[5][3]+_sm[3][7]*_sm[7][3];

	b_coef[2][1]=_sm[4][1]*_sm[1][3]+_sm[4][2]*_sm[2][3]+
			_sm[5][3]*_sm[4][5]+_sm[7][3]*_sm[4][7];

	b_coef[1][2]=_sm[3][5];

	b_coef[2][2]=_sm[4][5];

	t_CompVal s1 = 0.5*(b_coef[1][1]+b_coef[2][2]);
	t_CompVal s2 = sqrt(
			0.25*std::pow(b_coef[1][1]-b_coef[2][2],2)+
			b_coef[2][1]*b_coef[1][2]
		);

	// eigenvalues 
	lambda[0] = -sqrt(b_coef[0][0]);
	lambda[1] = -sqrt(s1+s2);
	lambda[2] = -sqrt(s1-s2);
	lambda[3] = lambda[0];
	
	for (int i=0; i<dim; i++) lambda[i+dim] = -lambda[i];

	{

		// asymptotics for direct problem

		for (int k=0; k<2; k++){

			{

				// z0, z3
				const int ind = 0+dim*k;
				Z[ind][0] = 1.0;
				Z[ind][1] = lambda[ind];
				Z[ind][2] = _sm[0][2]/lambda[ind];
				Z[ind][3] = 0.0;
				Z[ind][4] = 0.0;
				Z[ind][5] = 0.0;
				Z[ind][6] = 0.0;
				Z[ind][7] = 0.0;

			}

			for (int j=1; j<3; j++){

				//z1, z5, z2, z6

				// index of current fundamental vector z_ind
				const int ind = j+dim*k;

				t_Complex z2j, z3j, z4j, z6j;
				t_Complex l2 = lambda[ind]*lambda[ind];

				b_coef[1][0] = _sm[3][1]*b_coef[2][1] - _sm[4][1]*(b_coef[1][1]-l2);

				z3j = (l2 - _sm[0][1])*b_coef[2][1]/b_coef[1][0];

				z4j = -(b_coef[1][1] - l2)*(l2 - _sm[0][1])/b_coef[1][0];

				z6j = (_sm[3][7]*z3j + _sm[4][7]*z4j)/(l2 - _sm[6][7]);

				z2j = (_sm[0][2]+_sm[3][2]*z3j + _sm[4][2]*z4j + _sm[6][2]*z6j)/lambda[ind];

				Z[ind][0]= 1.0;
				Z[ind][1]= lambda[ind];
				Z[ind][2]= z2j;
				Z[ind][3] = z3j;
				Z[ind][4] = z4j;
				Z[ind][5] = lambda[ind]*z4j;
				Z[ind][6] = z6j;
				Z[ind][7] = lambda[ind]*z6j;

			}

			{

				// z3, z7
				const int ind = 3+dim*k;

				Z[ind][0] = 0.0;
				Z[ind][1] = 0.0;
				Z[ind][2] = _sm[6][2]/lambda[ind];
				Z[ind][3] = 0.0;
				Z[ind][4] = 0.0;
				Z[ind][5] = 0.0;
				Z[ind][6] = 1.0;
				Z[ind][7] = lambda[ind];

			}


		}

	}
	} */	//~ fundamental matrix Z of direct problem is set

	if (_ls_mode.is_flag_on(stab::t_LSMode::DIRECT)){

		// direct problem, take z1, z2, z3, z4 as initials

		for (int i=0; i<dim; i++)
			for (int j=0; j<STAB_MATRIX_DIM; j++) asym_vecs[i][j] = Z[i][j];


	}
	if(_ls_mode.is_flag_on(stab::t_LSMode::CONJUGATE)){

		// asymptotics for conjugate problem

		t_VecCmplx rhs(STAB_MATRIX_DIM);
		t_VecCmplx conj_vec(STAB_MATRIX_DIM);

		t_SqMatCmplx Z_tr(STAB_MATRIX_DIM);
		Z_tr = Z; Z_tr.transpose();

		// use 4 last vecs as initials
		for (int i=0; i<dim; i++){

			rhs.setToZero();
			rhs[i+dim] = 1.0;

			smat::solve_lsys_lu(Z_tr, rhs, conj_vec);

			// to solve "honest" conjugate problem
			//for (int j=0; j<STAB_MATRIX_DIM; j++) asym_vecs[i][j] = std::conj(conj_vec[j]);

			// to solve "conjugated" conjugate problem
			asym_vecs.set_col(i, conj_vec);

		}

	}

	// attempts of explicit conjugate vectors calcs
	// z~2
	/*
	{

			t_Complex l2 = std::conj(lambda[0])*std::conj(lambda[0]);

		initial_vectors[0][0] = 1.0;
		initial_vectors[0][1] = std::conj(lambda[0])/_sm[1][0];
		initial_vectors[0][2] = 0.0;
		initial_vectors[0][3] = 0.0;

		t_Complex ksi52 = (_sm[1][3]*_sm[7][4] - _sm[1][4]*_sm[7][3])/
			(_sm[7][3]*(_sm[5][4]+l2)-_sm[5][3]*_sm[7][4]);

		initial_vectors[0][4] = ksi52;

		initial_vectors[0][5] = -ksi52/std::conj(lambda[0]);

		t_Complex ksi72 = (_sm[5][3]*_sm[1][4]-_sm[1][3]*(_sm[5][4]+l2))/
			(_sm[7][3]*(_sm[5][4]+l2)-_sm[5][3]*_sm[7][4]);

		initial_vectors[0][6] = ksi72;
		initial_vectors[0][7] = -ksi72/std::conj(lambda[0]);
		}
	}
	*/

	
	// verification of asymptotics, enable when debug needed
	
	t_VecCmplx asym_resid_v(STAB_MATRIX_DIM, 0.0), init_vec(STAB_MATRIX_DIM, 0.0);
	t_VecCmplx v1(STAB_MATRIX_DIM, 0.0), v2(STAB_MATRIX_DIM, 0.0), v3(STAB_MATRIX_DIM, 0.0);
	double resid=0.0;

	// verify Z
	for (int i=0; i<STAB_MATRIX_DIM; i++){

			Z.col_to_vec(i, init_vec);
			//asym_resid_v = _stab_matrix*init_vec - lambda[i]*init_vec;		
			matrix::base::mat_mul<t_Complex, t_Complex>(_stab_matrix, init_vec, v1);
			matrix::base::mul<t_Complex, t_Complex>(lambda[i], init_vec, v2);
			matrix::base::minus<t_Complex, t_Complex>(v1, v2, asym_resid_v);

			resid = resid + asym_resid_v.norm().real();
	}

	//wxLogMessage(_T("Verify asymptotics v2 : resid_direct = %f"), resid);

	if (_ls_mode.is_flag_on(stab::t_LSMode::CONJUGATE)){

		t_SqMatCmplx H1(STAB_MATRIX_DIM);

		// H1 = -H_tr
		H1 = _stab_matrix;
		H1.mul_by_factor(-1.0);
		H1.transpose();

		resid = 0.0;

		for (int i=0; i<dim; i++){

			asym_vecs.col_to_vec(i, init_vec);

			//asym_resid_v = -1.0*_stab_matrix_tr*init_vec - lambda[i]*init_vec;		
			matrix::base::mat_mul<t_Complex, t_Complex>(H1, init_vec, v1);

			matrix::base::mul<t_Complex, t_Complex>(-lambda[i+dim], init_vec, v2);
			matrix::base::minus<t_Complex, t_Complex>(v1, v2, asym_resid_v);

			double dres = smat::norm(asym_resid_v.norm());
			resid = resid + dres;
		}

		//wxLogMessage(_T("Verify asymptotics v2 : resid_conjugate = %f"), resid);
		//std::wcout<<_T("Asym vecs:\n")<<asym_vecs<<_T("\n");


	}


}

t_Complex t_StabSolver::_calcResidual(t_Complex* out_resid_coefs) const{

	const std::vector<t_MatCmplx>& solution = _math_solver.solution;

	const t_MatCmplx& wall_func = solution[_math_solver.getNNodes()-1];

	t_VecCmplx rhs(4);
	t_VecCmplx resid_coefs(4);
	t_SqMatCmplx mat(4);

	t_Complex resid(0.0);

	if (!_ls_mode.is_flag_on(stab::t_LSMode::CONJUGATE)){

		// direct problem
		// get resid by temperature residual
		// construct matrix of 4 cols (u,u',v,w)
		// of solutions at wall
		// solve with rhs (0,1,0,0)
		// and construct resid of temperature

		rhs[1]=1.0;
		for (int i=0; i<4; i++){
			mat[i][0] = wall_func[i][0];
			mat[i][1] = wall_func[i][1];
			mat[i][2] = wall_func[i][2];
			mat[i][3] = wall_func[i][6];
		}
		matrix::base::mat_mul<t_CompVal, t_CompVal>(mat.inverse(),rhs,resid_coefs);

		for (int i=0; i<4; i++){
			resid+=resid_coefs[i]*wall_func[i][4];
		}

	}else{

		// conjugate problem
		// get resid by pressure residual
		// construct matrix of 4 cols (v,u',r,w')
		// of solutions at wall
		// solve with rhs (1,0,0,0)
		// and construct resid of pressure

		rhs[0]=1.0;
		for (int i=0; i<4; i++){
			mat[i][0] = wall_func[i][2];
			mat[i][1] = wall_func[i][1];
			mat[i][2] = wall_func[i][5];
			mat[i][3] = wall_func[i][7];
		}

		matrix::base::mat_mul<t_CompVal, t_CompVal>(mat.inverse(),rhs,resid_coefs);

		for (int i=0; i<4; i++){
			resid+=resid_coefs[i]*wall_func[i][3];
		}

	}

	// when needed write out residual coefs
	if (out_resid_coefs!=NULL)
		for (int i=0; i<4; i++)
			out_resid_coefs[i] = resid_coefs[i];

	return resid;

}

t_Complex t_StabSolver::solve(t_WCharsLoc& stab_point){
	this->_waveChars = stab_point;
	this->_math_solver.clean();

	_setAsymptotics_v2(_math_solver.solution[0]);

	_math_solver.solve();
	t_Complex resid = _calcResidual();
	_waveChars.resid = resid;
	stab_point.resid = resid;
	return resid;
}

void t_StabSolver::dumpEigenFuctions(const std::string& fname){

	std::vector<t_VecCmplx> amp_funcs(getNNodes(), t_VecCmplx(STAB_MATRIX_DIM));

	std::wofstream fstr(&fname[0]);
	fstr<<_T("u_re\tu_im\tu'_re\tu'_im\tv_re\tv_im\tp_re\tp_im\tt_re\tt_im\tt'_re\tt'_im\tw_re\tw_im\tw'_re\tw'_im\tY\n");

	getAmpFuncs(amp_funcs);

	for (int j=0; j<getNNodes(); j++){

		for (int k=0; k<STAB_MATRIX_DIM; k++){
			fstr<<std_manip::std_format_sci<double>(amp_funcs[j][k].real())
				<<_T("\t")
				<<std_manip::std_format_sci<double>(amp_funcs[j][k].imag())
				<<_T("\t");
		}
		fstr<<_math_solver.varRange[j];
		fstr<<_T("\n");
	}

}

void t_StabSolver::dumpProfileStab(const std::string& fname) const{
	_profStab.dump(fname);
}


