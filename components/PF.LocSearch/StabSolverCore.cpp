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

// Matrices HA and HW for scalar product of 2 amplitude functions

// The first matrix of scalar product is H1 = -i*dH0/da
// H0 is the basic stability matrix
// so this is to prepare calculations of <H1*x, z>
// x - direct amplitude vector, z - conjugate amplitude vector

// Addition matrix of scalar product is HW = +i*dH0/dw
// it is used in group velo calcs etc 

// Stability matrix corresponds to Nayfeh A.H., Stability of Three-Dimensional
// Boundary Layers," AIAA J., Vol. 18, No. 4, pp. 406-416, 1980.
// In this matrix:
//     m=MM=(e-1)2/3
//     r=RM=(e+2)2/3
//     e=K
//     e=(3/2)*(bulk viscosity)/(viscosity)

// Mack recommends (bulk viscosity)/(viscosity)=0.8, then parameter K=1.2
// TODO: cases where bv/v=0.8 ?
void t_StabSolver::_setScalProdMatrix_H1_HW(const t_ProfRec& rec){

	t_CompVal i(0.0, 1.0);

	const double& R = _profStab.scales().ReStab;
	const double& Me = _profStab.scales().Me;

	const mf::t_FldParams& Params = _rFldNS.get_mf_params();

	const double gMaMa = Params.Gamma*Me*Me;
	const double g_1MaMa = (Params.Gamma - 1.0)*Me*Me;
	const double Pr = Params.Pr;

	const double inv_t = 1.0 / rec.t;
	const double inv_mu = 1.0 / rec.mu;

	const double vCoefL = 2.0 / 3.0*(Params.BulkViscRatio + 2.0);
	const double vCoefS = 2.0 / 3.0*(Params.BulkViscRatio - 1.0);
	const double vCoefM = 1.0 + vCoefS;

	const t_CompVal& alpha = _waveChars.a;
	const t_CompVal& beta = _waveChars.b;
	const t_CompVal& freq = _waveChars.w;

	const double& u = rec.u;
	const double& u1 = rec.u1;

	const double& t = rec.t;
	const double& t1 = rec.t1;

	const double& w = rec.w;
	const double& w1 = rec.w1;

	const t_CompVal E = alpha*u + beta*w - freq;

	const double dMu1 = rec.mu1*t1;
	const double dMu2 = rec.mu2*t1*t1 + rec.mu1*rec.t2;
	const double dMu_U = rec.mu2*t1*u1 + rec.mu1*rec.u2;

	const t_CompVal xi = 1.0 / (R*inv_mu + i*vCoefL*gMaMa*E);

	const t_CompVal dxi_da = -i*xi*xi*vCoefL*gMaMa*u;
	const t_CompVal dxi_dw = i*xi*xi*vCoefL*gMaMa;

	_scal_prod_matrix_H1.setToZero();

	_scal_prod_matrix_H1[0][1] = R*u*inv_t*inv_mu - 2.0*i*alpha;

	_scal_prod_matrix_H1[2][1] = -inv_mu*dMu1 - vCoefM*inv_t*rec.t1;

	_scal_prod_matrix_H1[3][1] = R*inv_mu + i*vCoefM*gMaMa*(E + alpha*u);

	_scal_prod_matrix_H1[4][1] = -i*vCoefM*inv_t*(E + alpha*u);

	_scal_prod_matrix_H1[0][2] = -1.0;

	_scal_prod_matrix_H1[3][2] = -gMaMa*u;

	_scal_prod_matrix_H1[4][2] = inv_t*u;

	_scal_prod_matrix_H1[0][3] = -(2.0*inv_mu*dMu1 + vCoefL*inv_t*rec.t1)*(xi+ alpha*dxi_da);

	_scal_prod_matrix_H1[1][3] = -xi - alpha*dxi_da;

	_scal_prod_matrix_H1[2][3] = -i*(dxi_da*(-alpha*alpha - beta*beta +
		vCoefL*inv_t*(inv_mu*dMu1*rec.t1 + rec.t2) -
		i*R*inv_t*inv_mu*E) - xi*(2.0 * alpha + i*R*inv_mu*inv_t*u));

	_scal_prod_matrix_H1[3][3] = -vCoefL*gMaMa*(
		dxi_da*(alpha*u1 + beta*w1) + xi*u1 +
		(inv_mu*dMu1 + inv_t*rec.t1)*(dxi_da*E + xi*u));

	{
		const t_Complex b1 = inv_mu*rec.mu1 + vCoefL*inv_t;
		const t_Complex c1 = vCoefL*inv_t*inv_mu*dMu1;
		_scal_prod_matrix_H1[4][3] = dxi_da*(b1*(alpha*rec.u1 + beta*rec.w1)+c1*E) +
			xi*(b1*rec.u1 + c1*rec.u);
	}

	_scal_prod_matrix_H1[5][3] = vCoefL*inv_t*(dxi_da*E + xi*u);

	_scal_prod_matrix_H1[6][3] = -dxi_da*beta*(2.0*inv_mu*dMu1 + vCoefL*inv_t*t1);

	_scal_prod_matrix_H1[7][3] = -dxi_da*beta;

	_scal_prod_matrix_H1[2][5] = -2 * Pr*g_1MaMa*u1;

	_scal_prod_matrix_H1[3][5] = -R*Pr*inv_mu*g_1MaMa*u;

	_scal_prod_matrix_H1[4][5] = R*Pr*inv_mu*inv_t*u - 2.0*i*alpha;

	_scal_prod_matrix_H1[3][7] = i*vCoefM*beta*gMaMa*u;

	_scal_prod_matrix_H1[4][7] = -i*vCoefM*beta*inv_t*u;

	_scal_prod_matrix_H1[6][7] = R*inv_mu*inv_t*u - 2.0*i*alpha;

	// set HW

	_scal_prod_matrix_HW.setToZero();

	_scal_prod_matrix_HW[0][1] = R*inv_mu*inv_t;

	_scal_prod_matrix_HW[3][1] = i*alpha*vCoefM*gMaMa;

	_scal_prod_matrix_HW[4][1] = -i*vCoefM*inv_t*alpha;

	_scal_prod_matrix_HW[3][2] = -gMaMa;

	_scal_prod_matrix_HW[4][2] = inv_t;

	_scal_prod_matrix_HW[0][3] = alpha*dxi_dw*(2.0*inv_mu*dMu1 + vCoefL*inv_t*rec.t1);

	_scal_prod_matrix_HW[1][3] = alpha*dxi_dw;

	_scal_prod_matrix_HW[2][3] = i*dxi_dw*(-alpha*alpha - beta*beta +
		vCoefL*inv_t*inv_mu*dMu1*rec.t1 + 
		vCoefL*inv_t*rec.t2 - i*R*inv_t*inv_mu*E) - xi*R*inv_mu*inv_t;


	{
		t_Complex b1 = inv_mu*dMu1 + inv_t*rec.t1;
		t_Complex c1 = alpha*rec.u1 + beta*rec.w1;
		_scal_prod_matrix_HW[3][3] = vCoefL*gMaMa*(dxi_dw*(b1*E + c1) - xi*b1);

	}
	{
		t_Complex c1 = (inv_mu*rec.mu1 + vCoefL*inv_t)*(alpha*rec.u1 + beta*rec.w1);
		t_Complex b1 = vCoefL*inv_t*inv_mu*dMu1;
		_scal_prod_matrix_HW[4][3] = -dxi_dw*c1 - b1*(dxi_dw*E - xi);
	}

	_scal_prod_matrix_HW[5][3] = vCoefL*inv_t*(xi - dxi_dw*E);

	_scal_prod_matrix_HW[6][3] = beta*dxi_dw*(2.0*inv_mu*dMu1 + vCoefL*inv_t*rec.t1);

	_scal_prod_matrix_HW[7][3] = dxi_dw*beta;

	_scal_prod_matrix_HW[3][5] = -R*Pr*inv_mu*g_1MaMa;

	_scal_prod_matrix_HW[4][5] = R*Pr*inv_t*inv_mu;

	_scal_prod_matrix_HW[3][7] = i*vCoefM*beta*gMaMa;

	_scal_prod_matrix_HW[4][7] = -i*vCoefM*beta*inv_t;

	_scal_prod_matrix_HW[6][7] = R*inv_t*inv_mu;
		//*******************
		//*******************
}

void t_StabSolver::_setScalProdMatrix_H1_HW(const double& a_y){

	const t_ProfRec& rec = _profStab.get_rec(a_y);
	_setScalProdMatrix_H1_HW(rec);

}

// calculate H2 - matrix of non-parallel effects
// H2*z = H2_h*z + H2_i*dz/dy
// dz/dy = H0*z => H2 = H2_h + H2_i*H0

// here we use notation of Fedorov, AIAA-2002-2846
// x1,y1, z1 - global reference frame, nondim by Lref
// x,y,z - local reference frame, nondim by Dels

// H2*f gives vector in global reference frame;
// to get matrix in local rf, result is multiplied by dx/dx1 = L/Dels;
// in this case <H2*dze, ksi> gives addition in local reference frame
void t_StabSolver::_setScalProdMatrix_H2(const t_ProfRec& rec, const mf::t_RecGrad& rec_grad) {

	_setStabMatrix3D(rec);

	const mf::t_FldParams& Params = _rFldNS.get_mf_params();

	const double Gamma = Params.Gamma;
	const double Mach = Params.Mach;

	const double Pr = Params.Pr;

	double M2 = Gamma*Mach*Mach;

	double MG = (Gamma - 1.0)*Mach*Mach;

	const double R = _profStab.scales().ReStab;

	const double dx_dx1 = Params.L_ref / _profStab.scales().Dels;

	double U = rec.u;
	double V0 = rec.v *dx_dx1;
	double W = rec.w;

	double T = rec.t;
	double T_inv = 1.0 / T;
	double T_inv2 = T_inv*T_inv;

	double Mu = rec.mu;
	double Mu_inv = 1.0 / Mu;

	double dU_dx1 = rec_grad.ug[0] * dx_dx1;
	double dU_dy = rec_grad.ug[1];
	double dU_dz1 = rec_grad.ug[2] * dx_dx1;

	double dV_dy = rec_grad.vg[1];

	double dP_dx1 = rec_grad.pg[0] * dx_dx1;
	double dP_dz1 = rec_grad.pg[2] * dx_dx1;

	double dT_dx1 = rec_grad.tg[0] * dx_dx1;
	double dT_dy = rec_grad.tg[1];
	double dT_dz1 = rec_grad.tg[2] * dx_dx1;

	double dW_dx1 = rec_grad.wg[0] * dx_dx1;
	double dW_dy = rec_grad.wg[1];
	double dW_dz1 = rec_grad.wg[2] * dx_dx1;

	t_SqMatCmplx& H2 = _scal_prod_matrix_H2;

	t_SqMatCmplx& H2hom = _mat_tmp1;
	t_SqMatCmplx& H2inh = _mat_tmp2;

	H2.setToZero();
	H2hom.setToZero();
	H2inh.setToZero();

	double xii = U*dT_dx1 + V0*dT_dy + W*dT_dz1;

	// second row of homogen & non-homogen parts
	{
		double psi = U*dU_dx1 + V0*dU_dy + W*dU_dz1;

		double rmt = R*Mu_inv*T_inv;

		H2hom[0][1] = rmt*dU_dx1;
		H2hom[3][1] = rmt*M2*psi;
		H2hom[4][1] = -1.0*rmt*T_inv*psi;
		H2hom[6][1] = rmt*dU_dz1;

		H2inh[0][1] = rmt*V0;
	}

	// third
	{
		double div = dU_dx1 + dV_dy + dW_dz1;

		H2hom[0][2] = T_inv*dT_dx1;
		H2hom[3][2] = M2*(U*T_inv*dT_dx1 - div + V0*T_inv*dT_dy + W*T_inv*dT_dz1);
		H2hom[4][2] = T_inv2*(-2.0*xii + T*div);
		H2hom[6][2] = T_inv*dT_dz1;

		H2inh[3][2] = -1.0*V0*M2;
		H2inh[4][2] = T_inv*V0;
	}

	// fourth
	H2hom[2][3] = -1.0*T_inv*dV_dy;
	H2inh[2][3] = -1.0*T_inv*V0;

	// sixth
	{
		double rpm = R*Pr*Mu_inv;

		H2hom[0][5] = -1.0*rpm*MG*(dP_dx1 - T_inv*dT_dx1);
		H2hom[3][5] = rpm*M2*T_inv*xii;
		H2hom[4][5] = -1.0*rpm*T_inv2*xii;
		H2hom[6][5] = -1.0*rpm*MG*(dP_dz1 - T_inv*dT_dz1);

		H2inh[3][5] = -1.0*rpm*MG*V0;
		H2inh[4][5] = rpm*T_inv*V0;
	}
	// 8
	{
		double rmt = R*Mu_inv*T_inv;

		double eta = U*dW_dx1 + V0*dW_dy + W*dW_dz1;

		H2hom[0][7] = rmt*dW_dx1;
		H2hom[3][7] = rmt*M2*T_inv*eta;
		H2hom[4][7] = -1.0*rmt*T_inv*eta;
		H2hom[6][7] = rmt*dW_dz1;

		H2inh[6][7] = rmt*V0;
	}
	// m3 = H2inh*H0
	matrix::base::mat_mul<t_Complex, t_Complex>(H2inh, _stab_matrix, _mat_tmp3);
	// H2 = H2hom + m3
	matrix::base::plus<t_Complex, t_Complex>(H2hom, _mat_tmp3, _scal_prod_matrix_H2);

	// computed matrix is for x1 reference frame, return to local (x) reference
	_scal_prod_matrix_H2.mul_by_factor(1.0/dx_dx1);

	// debug, using only homogen part
	//_scal_prod_matrix_H2 = H2hom;

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
		// construct matrix of 4 cols (f1,f2,f3,f4)
		// of solutions at wall
		// solve with rhs (0,1,0,0)
		// and construct resid of temperature

		rhs[1]=1.0;

		// fill in matrix of values f[i][j]
		{
			// homogen BC, f1 = u, f2 = u', f3 = v, f4 = w
			if (_params.WallBC == pf::t_StabSolverParams::WALL_HOMOGEN) {
				for (int i = 0; i < 4; i++) {
					mat[i][0] = wall_func[i][0];
					mat[i][1] = wall_func[i][1];
					mat[i][2] = wall_func[i][2];
					mat[i][3] = wall_func[i][6];
				}
			}

			// slip BC, f1 = u - eta_u*u', f2 = u', f3 = v, f4 = w - eta_w*w'
			if (_params.WallBC == pf::t_StabSolverParams::WALL_SLIP) {
				// WallBC_EtaU, WallBC_EtaW are in global HSFlow nondim scale
				// in stability reference frame (length scale is Dels)
				// coefs are multiplied by Lref/Dels
				
				double l_coef = _rFldNS.get_mf_params().L_ref / get_stab_scales().Dels;
				double eu = _params.WallBC_EtaU * l_coef;
				double ew = _params.WallBC_EtaW * l_coef;
				for (int i = 0; i < 4; i++) {
					mat[i][0] = wall_func[i][0] - eu*wall_func[i][1];
					mat[i][1] = wall_func[i][1];
					mat[i][2] = wall_func[i][2];
					mat[i][3] = wall_func[i][6] - ew*wall_func[i][7];
				}
			}
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

		// IMPORTANT TODO: conditions for slip bc
		if (_params.WallBC != pf::t_StabSolverParams::WALL_HOMOGEN) {
			wxLogError(_T("t_StabSolver::_calcResidual not implemented for conjugate with this wall bc type"));
		}
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

#ifdef __PERFORMANCE_MEASURE
	t_TimeInterval::log(_T("t_StabSolver::solve:start"));
#endif // __PERFORMANCE_MEASURE

	this->_waveChars = stab_point;
	this->_math_solver.clean();

	setAsymptotics(_math_solver.solution[0]);

#ifdef __PERFORMANCE_MEASURE
	t_TimeInterval::log(_T("t_StabSolver::solve:asymptotic set"));
#endif // __PERFORMANCE_MEASURE

	_math_solver.solve();

#ifdef __PERFORMANCE_MEASURE
	//t_TimeInterval::log(_T("t_StabSolver::solve:_solve done"));
#endif // __PERFORMANCE_MEASURE

	t_Complex resid = _calcResidual();

#ifdef __PERFORMANCE_MEASURE
	//t_TimeInterval::log(_T("t_StabSolver::solve:resid computed"));
#endif // __PERFORMANCE_MEASURE

	_waveChars.resid = resid;
	stab_point.resid = resid;
	return resid;
}

void t_StabSolver::dumpEigenFuctions(const std::string& fname){

	std::vector<t_VecCmplx> amp_funcs(getNNodes(), t_VecCmplx(STAB_MATRIX_DIM));

	getAmpFuncs(amp_funcs);

	stab::dumpEigenFuncs(fname, getNNodes(), _math_solver.varRange, amp_funcs);

}

void t_StabSolver::dumpProfileStab(const std::string& fname) const{
	_profStab.dump(fname);
}


