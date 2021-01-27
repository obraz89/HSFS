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

	const double m13 = _curv_coefs.m13;
	const double m31 = _curv_coefs.m31;
	const double m12 = _curv_coefs.m12;
	const double m32 = _curv_coefs.m32;

	const t_CompVal xi = 1.0 / (stabRe*inv_mu + imagUnity*vCoefL*gMaMa*dEicon_dt);
	//const t_CompVal xi_m = 1.0 / (stabRe*inv_mu + imagUnity*vCoefL*gMaMa*dEicon_dt + vCoefL*gMaMa*m_e3);

	//const t_CompVal xi_dd = xi_m - xi;

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

	_stab_matrix[6][1] = 0.0;

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
	_stab_matrix[0][7] = 0.0;

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

	const bool curv_on = _params.CurvTermsOn;

	// compressible addition to curv terms
	/*
	if (curv_on) {

		const double m_e1 = m12 + m32;

		const double m_e3 = m31*rec.u + m13*rec.w;

		const double m_e4 = m12*rec.u*rec.u + m32*rec.w*rec.w;

		// second row
		_stab_matrix[0][1] += stabRe*inv_mu*inv_t*rec.w*m12 +
			vCoefM*imagUnity*alpha*m21;

		_stab_matrix[2][1] += stabRe*inv_mu*inv_t*rec.u*m13 +
			vCoefM*imagUnity*alpha*m_e1;

		const double u_e2 = rec.u*rec.w*m12 - rec.w*rec.w*m21;

		_stab_matrix[3][1] += stabRe*inv_mu*inv_t*gMaMa*u_e2 +
			vCoefM*imagUnity*alpha*gMaMa*m_e3;

		_stab_matrix[4][1] += -stabRe*inv_mu*inv_t*inv_t*u_e2 -
			vCoefM*imagUnity*alpha*inv_t*m_e3;

		_stab_matrix[6][1] += stabRe*inv_mu*inv_t*(rec.u*m12 - 2 * rec.w*m21) +
			vCoefM*imagUnity*alpha*m12;

		// thrid row

		_stab_matrix[0][2] += -m21;

		_stab_matrix[2][2] += -m_e1;

		_stab_matrix[3][2] += -gMaMa*m_e3;

		_stab_matrix[4][2] += inv_t*m_e3;

		_stab_matrix[6][2] += -m12;

		// fourth row
		// IMPORTANT TODO: exact expressions for curvature terms
		// for now assume xi = xi_m, xi_dd = 0, 
		// this is the case when m12=m21=0

		_stab_matrix[0][3] += xi*(-2.0*stabRe*inv_mu*inv_t*rec.u*m13 +
			vCoefL*imagUnity*alpha*m_e1);

		_stab_matrix[2][3] +=
			xi*vCoefL*m_e1*(-2.0*inv_t*rec.t1 + m_e1 - inv_mu*dMu1);

		_stab_matrix[3][3] += xi*gMaMa*(stabRe*inv_mu*inv_t*m_e4 +
			vCoefL*m_e1*imagUnity*dEicon_dt);

		_stab_matrix[4][3] += -xi*inv_t*(stabRe*inv_mu*inv_t*m_e4 +
			vCoefL*m_e1*imagUnity*dEicon_dt);

		_stab_matrix[6][3] += xi*(2.0*stabRe*inv_mu*inv_t*rec.w*m23 +
			vCoefL*imagUnity*beta*m_e1);

		// eight row
		_stab_matrix[0][7] += stabRe*inv_mu*inv_t*(rec.w*m21 - 2.0*rec.u*m12) +
			vCoefM*imagUnity*beta*m21;

		_stab_matrix[2][7] += stabRe*inv_mu*inv_t*rec.w*m23 +
			vCoefM*imagUnity*beta*m_e1;

		const double w_e2 = rec.u*rec.w*m21 - rec.u*rec.u*m12;

		_stab_matrix[3][7] += stabRe*inv_mu*inv_t*gMaMa*w_e2 +
			vCoefM*imagUnity*beta*gMaMa*m_e3;

		_stab_matrix[4][7] += -stabRe*inv_mu*inv_t*inv_t*w_e2 -
			vCoefM*imagUnity*beta*inv_t*m_e3;

		_stab_matrix[6][7] += stabRe*inv_mu*inv_t*rec.u*m21 +
			vCoefM*imagUnity*beta*m12;

	} */ // Compressible curvature additions to stability matrix 

	// incompressible additions after Cebeci
	if (curv_on) {

		//wxLogMessage(_T("LS: Curv terms on!"));

		// second row, 'du/dy'-eq
		_stab_matrix[0][1] += stabRe*rec.w*m13;

		_stab_matrix[2][1] += stabRe*rec.u*m12;

		_stab_matrix[6][1] += stabRe*(rec.u*m13 - 2 * rec.w*m31);

		// third row, 'v'-eq

		_stab_matrix[0][2] += -m31;

		_stab_matrix[2][2] += -(m12+m32);

		_stab_matrix[6][2] += -m13;

		// fourth row, 'p'-eq or 'dv/dz'-eq

		double m2 = m12 + m32;

		_stab_matrix[0][3] += xi*(2.0*stabRe*inv_mu*inv_t*rec.u*m12 +
			vCoefL*imagUnity*alpha*m2);

		_stab_matrix[2][3] +=
			xi*vCoefL*m2*(m2-inv_t*rec.t1 - inv_mu*dMu1);

		const double m_e4 = m12*rec.u*rec.u + m32*rec.w*rec.w;

		_stab_matrix[3][3] += xi*gMaMa*(stabRe*inv_mu*inv_t*m_e4 +
			vCoefL*m2*imagUnity*dEicon_dt);

		_stab_matrix[4][3] += -xi*inv_t*(stabRe*inv_mu*inv_t*m_e4 +
			vCoefL*m2*imagUnity*dEicon_dt);

		_stab_matrix[6][3] += xi*(2.0*stabRe*inv_mu*inv_t*rec.w*m32 +
			vCoefL*imagUnity*beta*m2);

		// eighth row, 'dw/dy'-eq

		// eight row
		_stab_matrix[0][7] += stabRe*(rec.w*m31 - 2.0*rec.u*m13);

		_stab_matrix[2][7] += stabRe*rec.w*m32;

		_stab_matrix[6][7] += stabRe*rec.u*m31;

	}

}

void t_StabSolver::_setCurvCoefs(const mf::t_GeomPoint& a_xyz) {

	// sphere case, dimensional raduis
	const double R_dim = 0.089;

	_curv_coefs.set_coefs_sphere(_profStab.scales(), R_dim);

};

void t_StabSolver::_setStabMatrix3D(const double& a_y){

	const t_ProfRec& rec = _profStab.get_rec(a_y);
	_setStabMatrix3D(rec);
};

// Matrix Hx for scalar product of 2 amplitude functions

// The matrix of scalar product is H1 = -i*dH0/da
// H0 is the basic stability matrix
// so this is to prepare calculations of <H1*x, z>
// x - direct amplitude vector, z - conjugate amplitude vector
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
void t_StabSolver::_setScalProdMatrix_H1(const t_ProfRec& rec){

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

		_scal_prod_matrix_H1.setToZero();

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

		_scal_prod_matrix_H1[0][1] = U*R*M1*TD-2.*E*A;

		_scal_prod_matrix_H1[2][1]=-F*T1*TD-MY*M1;

		_scal_prod_matrix_H1[3][1]=R*M1-E*F*M2*(WA-A*U);

		_scal_prod_matrix_H1[4][1]=E*F*TD*(WA-A*U);

		_scal_prod_matrix_H1[0][2]=-1.;

		_scal_prod_matrix_H1[3][2]=-M2*U;

		_scal_prod_matrix_H1[4][2]=U*TD;

		_scal_prod_matrix_H1[0][3]=-(XI+A*DXI)*(RM*T1*TD+2.*MY*M1);

		_scal_prod_matrix_H1[1][3]=-XI-A*DXI;

		_scal_prod_matrix_H1[2][3]=-E*(
			DXI*(E*WA*R*M1*TD-A*A-B*B+RM*(T2*TD+MY*T1*M1*TD))
			-XI*(2.*A+E*R*U*M1*TD));

		_scal_prod_matrix_H1[3][3]=-RM*M2*
			(DXI*(A*U1+B*WS1)+XI*U1+(T1*TD+MY*M1)*(XI*U-DXI*WA));

		_scal_prod_matrix_H1[4][3]=DXI*((A*U1+B*WS1)*(RM*TD+MU1*M1)-RM*WA*MY*M1*TD)
			+XI*(RM*U1*TD+MU1*U1*M1+RM*U*MY*M1*TD);

		_scal_prod_matrix_H1[5][3]=RM*TD*(U*XI-WA*DXI);

		_scal_prod_matrix_H1[6][3]=-DXI*B*(RM*T1*TD+2.*MY*M1);

		_scal_prod_matrix_H1[7][3]=-B*DXI;

		_scal_prod_matrix_H1[2][5]=-2.*MG*Pr*U1;

		_scal_prod_matrix_H1[3][5]=-MG*Pr*R*U*M1;

		_scal_prod_matrix_H1[4][5]=R*Pr*U*M1*TD-2.*E*A;

		_scal_prod_matrix_H1[3][7]=E*F*M2*B*U;

		_scal_prod_matrix_H1[4][7]=-E*F*B*U*TD;

		_scal_prod_matrix_H1[6][7]=R*U*M1*TD-2.*E*A;

}

void t_StabSolver::_setScalProdMatrix_H1(const double& a_y){

	const t_ProfRec& rec = _profStab.get_rec(a_y);
	_setScalProdMatrix_H1(rec);

}

// calculate H2*z = {G1, ..., G8} - vector of non-par effects of mean flow
// H2*z must be computed explicitly as {G1, ..., G8} as it contains non-homogen part (dz3/dy and dz4/dy terms)
// NB : dz1/dy = z2, dz5/dy = z6, dz7/dy = z8
// dz3/dy and dz4/dy must be computed
// NB : rec_grad.ug[1] and rec.u1 are the same value : dU_dy, check for verification
// NB : currently terms d_dz not included (non-par effects for 2D flow)
void t_StabSolver::_calc_H2z(int i, std::vector<t_VecCmplx>& fun_direct, t_VecCmplx& H2z) {

	// amp funcs in reverse order, get
	int nnodes = _math_solver.getNNodes();
	int ind_r = nnodes - 1 - i;

	const std::vector<double>& y_distrib = get_y_distrib();

	const t_ProfRec& rec = _profStab.get_rec(i);

	const mf::t_RecGrad& rec_grad = _profStab.get_rec_grad(i);

	const t_CompVal R = _profStab.scales().ReStab;
	
	t_CompVal E(0.0, 1.0);

	t_CompVal z1 = fun_direct[ind_r][0];
	t_CompVal z2 = fun_direct[ind_r][1];
	t_CompVal z3 = fun_direct[ind_r][2];
	t_CompVal z4 = fun_direct[ind_r][3];
	t_CompVal z5 = fun_direct[ind_r][4];
	t_CompVal z6 = fun_direct[ind_r][5];
	t_CompVal z7 = fun_direct[ind_r][6];
	t_CompVal z8 = fun_direct[ind_r][7];

	t_CompVal z1_y = z2;

	t_CompVal z3_y, z4_y;
	_calc_amp_fun_dv_dy_dp_dy(ind_r, fun_direct, z3_y, z4_y);

	t_CompVal z5_y = z6;
	t_CompVal z7_y = z8;

	const mf::t_FldParams& Params = _rFldNS.get_mf_params();

	const double Gamma = Params.Gamma;
	const double Mach = Params.Mach;

	const double Pr = Params.Pr;

	double M2 = Gamma*Mach*Mach;

	double MG = (Gamma - 1.0)*Mach*Mach;

	double U = rec.u;
	double V = rec.v;
	double W = rec.w;

	double T = rec.t;
	double T_inv = 1.0 / T;
	double T_inv2 = T_inv*T_inv;

	double Mu = rec.mu;
	double Mu_inv = 1.0 / Mu;

	double dU_dx = rec_grad.ug[0];
	double dU_dy = rec_grad.ug[1];

	double dV_dy = rec_grad.vg[1];

	double dP_dx = rec_grad.pg[0];

	double dT_dx = rec_grad.tg[0];
	double dT_dy = rec_grad.tg[1];

	t_Complex Fun1 = M2*T_inv*z4 - T_inv2*z5;

	H2z.setToZero();
	// first row - zero
	 
	H2z[1] = R*Mu_inv*(z1*T_inv*dU_dx + V*T_inv*z2 + Fun1*(U*dU_dx + V*dU_dy));
	 
	H2z[2] = z1*T_inv*dT_dx + U*M2*z4*T_inv*dT_dx - 2*U*z5*T_inv2*dT_dx - (dU_dx + dV_dy)*Fun1 
		- V*(M2*z4_y - M2*T_inv*z4*dT_dy - T_inv*z5_y + 2.0 * z5*T_inv2*dT_dy);
	
	H2z[3] = -T_inv*(V*z3_y + z3*dV_dy);

	H2z[5] = -R*Pr*Mu_inv*(-V*T_inv*z5_y - Fun1*(U*dT_dx + V*dT_dy) + MG*(dP_dx*z1 + V*z4_y - z1*T_inv*dT_dx));

	H2z[7] = R*Mu_inv*V*T_inv*z7_y;

}

void t_StabSolver::setAsymptotics(t_MatCmplx& asym_vecs, t_CompVal* a_lambdas /*=NULL*/){

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

		if (a_lambdas != NULL) {
			for (int i = 0; i < 4; i++) a_lambdas[i] = lambda[i];
		}

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

	setAsymptotics(_math_solver.solution[0]);

	_math_solver.solve();
	t_Complex resid = _calcResidual();
	_waveChars.resid = resid;
	stab_point.resid = resid;
	return resid;
}

void t_StabSolver::dumpEigenFuctions(const std::string& fname){

	std::vector<t_VecCmplx> amp_funcs(getNNodes(), t_VecCmplx(STAB_MATRIX_DIM));

	std::wofstream fstr(&fname[0]);
	fstr<<_T("Y\tu_re\tu_im\tu'_re\tu'_im\tv_re\tv_im\tp_re\tp_im\tt_re\tt_im\tt'_re\tt'_im\tw_re\tw_im\tw'_re\tw'_im\n");

	getAmpFuncs(amp_funcs);

	// write out in reverse order (from wall to outer)
	for (int j=getNNodes()-1; j>=0; j--){

		fstr<<_math_solver.varRange[j]<<_T("\t");

		for (int k=0; k<STAB_MATRIX_DIM; k++){
			fstr<<std_manip::std_format_sci<double>(amp_funcs[j][k].real())
				<<_T("\t")
				<<std_manip::std_format_sci<double>(amp_funcs[j][k].imag())
				<<_T("\t");
		}
		fstr<<_T("\n");
	}

}

void t_StabSolver::dumpProfileStab(const std::string& fname) const{
	_profStab.dump(fname);
}


