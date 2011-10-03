#include "StabSolver.h"
//#include "SolverCore.h"
t_StabSolver::t_StabSolver(const MF_Field& a_rFldNS, t_StabField& a_rFldStab):
_rFldNS(a_rFldNS), _rFldStab(a_rFldStab), _profStab(0), 
// default is 3D - 8x8 stab matrix
_math_solver(), _stab_matrix(8){
	_math_solver._pStab_solver = this;
};
// obsolete
/*
void t_StabSolver::setParameters(t_ProfileStab& a_profStab){
	_profile = a_profStab;
	// globals -flow
	DNS.AMINF = MF_Field::Mach;
	DNS.REINF = MF_Field::Re;
	DNS.TINF = MF_Field::T_inf;
	DNS.XLL = MF_Field::L_ref;
	// -gas parameters
	OUT.G = 1.4;	// specific heat ratio
	OUT.SIG = 0.72;	// Prandtl number
	OUT.K = 0.0;	  // second viscosity
	OUT.VISC_FLAG = 1; // for viscosity: iupt=0 - Sutherland, 1 - power
	OUT.VISC_POWER = 0.75; // if VISC_FLAG=1, VISC_POWER - viscosity power coef  
	// locals
	NS.XST[0] = fld_ref.fld[i_ind][0][k_ind].x;
	CF_PROFILE_DATA.MF_WE = profile.w_e;
	CF_PROFILE_DATA.MF_ANG = fld_ref.get_cf_wave_dir(i_ind, k_ind);
	ADDITIONAL_NS.Y_DIM = profile.ny;
	MACKY1.NO = 150;	// number of grid points for stability computations
	HADY1.NS = 1;		// coef increasing number of points for stability 
	//computations (adaptation by NAVSTOK on uniform grid) 
	OUT.XI = 0.0;	// notused?
	OUT.PB = 0.0;	// ??
	ABOCT.MAB = 1;	// direct problem
};*/

// private, to be used in 3D context

void t_StabSolver::_setStabMatrix3D(const double& a_y){

	//t_SqMatrix stab_matrix(8);
	t_CompVal imagUnity(0.0, 1.0);

	const double& stabRe = _profStab.stabRe;
	const double& Me = _profStab.Me;
	const double gMaMa = MF_Field::Gamma*Me*Me;
	const double g_1MaMa = (MF_Field::Gamma-1.0)*Me*Me;

	const double u = _profStab.getValue(a_y, _profStab.u);
	const double u1 = _profStab.getValue(a_y, _profStab.u1);
	const double u2 = _profStab.getValue(a_y, _profStab.u2);

	const double w = _profStab.getValue(a_y, _profStab.w);
	const double w1 = _profStab.getValue(a_y, _profStab.w1);
	const double w2 = _profStab.getValue(a_y, _profStab.w2);

	const double t = _profStab.getValue(a_y, _profStab.t);
	const double t1 = _profStab.getValue(a_y, _profStab.t1);
	const double t2 = _profStab.getValue(a_y, _profStab.t2);

	const double mu = _profStab.getValue(a_y, _profStab.mu);
	const double mu1 = _profStab.getValue(a_y, _profStab.mu1);
	const double mu2 = _profStab.getValue(a_y, _profStab.mu2);

	const double inv_t = 1.0/t;
	const double inv_mu = 1.0/mu;

	const double vCoefL = 2.0/3.0*(0.0+2.0); //TODO: second Visc coef instead of 0.0
	const double vCoefS = 2.0/3.0*(0.0-1.0); // 
	const double vCoefM = 1.0+vCoefS;

	const t_CompVal& alpha = _waveChars.a;
	const t_CompVal& beta = _waveChars.b;
	const t_CompVal& freq = _waveChars.w;

	const t_CompVal dEicon_dt = alpha*u + beta*w - freq;
	const double dMu2 = mu2*t1*u1 + mu1*u2;
	const t_CompVal xi = 1.0/(stabRe*inv_mu+imagUnity*vCoefL+gMaMa*dEicon_dt);

	// first
	_stab_matrix[1][0]=1.0;
	// second
	_stab_matrix[0][1] = imagUnity*stabRe*inv_t*inv_mu*dEicon_dt +
		pow(alpha,2) + pow(beta,2);

	_stab_matrix[1][1] =  - inv_mu*mu1*t1;

	_stab_matrix[2][1] = stabRe*inv_t*inv_mu*u1 - 
		imagUnity*alpha*(inv_mu*mu1*t1+vCoefM*inv_t*t1);

	_stab_matrix[3][1] = alpha*(imagUnity*stabRe*inv_mu -
		vCoefM*gMaMa*dEicon_dt);

	_stab_matrix[4][1] = vCoefM*alpha*inv_t*dEicon_dt - inv_mu*dMu2;

	_stab_matrix[5][1] = -inv_mu*mu1*u1;
	// third
	_stab_matrix[0][2] = -imagUnity*alpha;

	_stab_matrix[2][2] = inv_t*t1;

	_stab_matrix[3][2] = -imagUnity*gMaMa*dEicon_dt;

	_stab_matrix[4][2] = imagUnity*inv_t*dEicon_dt;

	_stab_matrix[6][2] = -imagUnity*beta;
	// fourth
	_stab_matrix[0][3] = -imagUnity*xi*alpha*t1*
		(2.0*inv_mu*mu1 + vCoefL*inv_t);

	_stab_matrix[1][3] = -imagUnity*xi*alpha;

	_stab_matrix[2][3] = xi*(-alpha*alpha-beta*beta+
		vCoefL*inv_t*inv_mu*mu1*t1*t1+
		vCoefL*inv_t*t2-
		imagUnity*stabRe*inv_t*inv_mu*dEicon_dt);

	_stab_matrix[3][3] = -imagUnity*xi*vCoefL*gMaMa*
		(
		(inv_mu*mu1*t1+inv_t*t1)*dEicon_dt+
		alpha*u1+beta*w1
		);

	_stab_matrix[4][3] = imagUnity*xi*
		(
		(inv_mu*mu1+vCoefL*inv_t)*(alpha*u1+beta*w1)+
		vCoefL*inv_t*inv_mu*mu1*t1*dEicon_dt
		);

	_stab_matrix[5][3] = imagUnity*xi*vCoefL*inv_t*dEicon_dt;

	_stab_matrix[6][3] = -imagUnity*xi*beta*
		(2.0*inv_mu*mu1*t1 + vCoefL*inv_t*t1);

	_stab_matrix[7][3] = -imagUnity*xi*beta;
	// fifth
	_stab_matrix[5][4]=1.0;
	// sixth
	_stab_matrix[1][5] = -2.0*MF_Field::Pr*g_1MaMa*u1;

	_stab_matrix[2][5] = MF_Field::Pr*
		(
		stabRe*inv_t*inv_mu*t1 - 
		2.0*imagUnity*g_1MaMa*(alpha*u1+beta*w1)
		);

	_stab_matrix[3][5] = -imagUnity*stabRe*MF_Field::Pr*inv_mu*g_1MaMa*dEicon_dt;

	_stab_matrix[4][5] = imagUnity*stabRe*MF_Field::Pr*inv_t*inv_mu*dEicon_dt +
		alpha*alpha + beta*beta - 
		g_1MaMa*MF_Field::Pr*inv_mu*mu1*(u1*u1+w1*w1)-
		inv_mu*(mu2*t1*t1+mu1*t2);
	_stab_matrix[5][5] = -2.0*inv_mu*mu1*t1;

	_stab_matrix[7][5] = -2.0*MF_Field::Pr*g_1MaMa*w1;
	// seventh
	_stab_matrix[7][6]=1.0;
	// last
	_stab_matrix[2][7] = -imagUnity*beta*(inv_mu*mu1*t1+vCoefM*inv_t*t1)+
		stabRe*inv_mu*inv_t*w1;

	_stab_matrix[3][7] = imagUnity*stabRe*beta*inv_mu-
		vCoefM*beta*gMaMa*dEicon_dt;

	_stab_matrix[4][7] = vCoefM*beta*inv_t*dEicon_dt-
		inv_mu*(mu2*t1*w1+mu1*w2);

	_stab_matrix[5][7] = -inv_mu*mu1*w1;

	_stab_matrix[6][7] = imagUnity*stabRe*inv_t*inv_mu*dEicon_dt +
		alpha*alpha + beta*beta;

	_stab_matrix[7][7] = -inv_mu*mu1*t1;
};

/*const t_SqMatrix& t_StabSolver::_getStabMatrix3D() const{
	return _stab_matrix;
}*/

// old - can be used only in local spectral procedures
/*void t_StabSolver::setMatSplitByW3D_88(const double& a_y, t_SqMatrix& mat_no_w, t_SqMatrix& mat_w) const{
	t_CompVal imagUnity(0.0, 1.0);

	const double& stabRe = _profStab.stabRe;
	const double& Me = _profStab.Me;
	const double gMaMa = MF_Field::Gamma*Me*Me;
	const double g_1MaMa = (MF_Field::Gamma-1.0)*Me*Me;

	const double u = _profStab.getValue(a_y, _profStab.u);
	const double u1 = _profStab.getValue(a_y, _profStab.u1);
	const double u2 = _profStab.getValue(a_y, _profStab.u2);

	const double w = _profStab.getValue(a_y, _profStab.w);
	const double w1 = _profStab.getValue(a_y, _profStab.w1);
	const double w2 = _profStab.getValue(a_y, _profStab.w2);

	const double t = _profStab.getValue(a_y, _profStab.t);
	const double t1 = _profStab.getValue(a_y, _profStab.t1);
	const double t2 = _profStab.getValue(a_y, _profStab.t2);

	const double mu = _profStab.getValue(a_y, _profStab.mu);
	const double mu1 = _profStab.getValue(a_y, _profStab.mu1);
	const double mu2 = _profStab.getValue(a_y, _profStab.mu2);

	const double inv_t = 1.0/t;
	const double inv_mu = 1.0/mu;

	const double vCoefL = 2.0/3.0*(0.0+2.0); //TODO: second Visc coef instead of 0.0
	const double vCoefS = 2.0/3.0*(0.0-1.0); // 
	const double vCoefM = 1.0+vCoefS;

	const t_CompVal& alpha = _waveChars.a;
	const t_CompVal& beta = _waveChars.b;

	const t_CompVal dEicon_no_W = alpha*u + beta*w ;
	const t_CompVal dEicon_W = -1.0;
	const double dMu2 = mu2*t1*u1 + mu1*u2;
	const t_CompVal xi = 1.0/(stabRe*inv_mu+imagUnity*vCoefL+gMaMa*dEicon_no_W);
// set no_W matrix
// of size 5rowsx8cols
	// first <--> "2"
	mat_no_w[0][0] = imagUnity*stabRe*inv_t*inv_mu*dEicon_no_W +
		pow(alpha,2) + pow(beta,2);

	mat_no_w[1][0] =  - inv_mu*mu1*t1;

	mat_no_w[2][0] = stabRe*inv_t*inv_mu*u1 - 
		imagUnity*alpha*(inv_mu*mu1*t1+vCoefM*inv_t*t1);

	mat_no_w[3][0] = alpha*(imagUnity*stabRe*inv_mu -
		vCoefM*gMaMa*dEicon_no_W);

	mat_no_w[4][0] = vCoefM*alpha*inv_t*dEicon_no_W - inv_mu*dMu2;

	mat_no_w[5][0] = -inv_mu*mu1*u1;
	// second <--> "3"
	mat_no_w[0][1] = -imagUnity*alpha;

	mat_no_w[2][1] = inv_t*t1;

	mat_no_w[3][1] = -imagUnity*gMaMa*dEicon_no_W;

	mat_no_w[4][1] = imagUnity*inv_t*dEicon_no_W;

	mat_no_w[6][1] = -imagUnity*beta;
	// third <--> "4"
	mat_no_w[0][2] = -imagUnity*xi*alpha*t1*
		(2.0*inv_mu*mu1 + vCoefL*inv_t);

	mat_no_w[1][2] = -imagUnity*xi*alpha;

	mat_no_w[2][2] = xi*(-alpha*alpha-beta*beta+
		vCoefL*inv_t*inv_mu*mu1*t1*t1+
		vCoefL*inv_t*t2-
		imagUnity*stabRe*inv_t*inv_mu*dEicon_no_W);

	mat_no_w[3][2] = -imagUnity*xi*vCoefL*gMaMa*
		(
		(inv_mu*mu1*t1+inv_t*t1)*dEicon_no_W+
		alpha*u1+beta*w1
		);

	mat_no_w[4][2] = imagUnity*xi*
		(
		(inv_mu*mu1+vCoefL*inv_t)*(alpha*u1+beta*w1)+
		vCoefL*inv_t*inv_mu*mu1*t1*dEicon_no_W
		);

	mat_no_w[5][2] = imagUnity*xi*vCoefL*inv_t*dEicon_no_W;

	mat_no_w[6][2] = -imagUnity*xi*beta*
		(2.0*inv_mu*mu1*t1 + vCoefL*inv_t*t1);

	mat_no_w[7][2] = -imagUnity*xi*beta;
	// fourth <--> "6"
	mat_no_w[1][3] = -2.0*MF_Field::Pr*g_1MaMa*u1;

	mat_no_w[2][3] = MF_Field::Pr*
		(
		stabRe*inv_t*inv_mu*t1 - 
		2.0*imagUnity*g_1MaMa*(alpha*u1+beta*w1)
		);

	mat_no_w[3][3] = -imagUnity*stabRe*MF_Field::Pr*inv_mu*g_1MaMa*dEicon_no_W;

	mat_no_w[4][3] = imagUnity*stabRe*MF_Field::Pr*inv_t*inv_mu*dEicon_no_W +
		alpha*alpha + beta*beta - 
		g_1MaMa*MF_Field::Pr*inv_mu*mu1*(u1*u1+w1*w1)-
		inv_mu*(mu2*t1*t1+mu1*t2);
	mat_no_w[5][3] = -2.0*inv_mu*mu1*t1;

	mat_no_w[7][3] = -2.0*MF_Field::Pr*g_1MaMa*w1;
	// fifth <--> "8"
	mat_no_w[2][4] = -imagUnity*beta*(inv_mu*mu1*t1+vCoefM*inv_t*t1)+
		stabRe*inv_mu*inv_t*w1;

	mat_no_w[3][4] = imagUnity*stabRe*beta*inv_mu-
		vCoefM*beta*gMaMa*dEicon_no_W;

	mat_no_w[4][4] = vCoefM*beta*inv_t*dEicon_no_W-
		inv_mu*(mu2*t1*w1+mu1*w2);

	mat_no_w[5][4] = -inv_mu*mu1*w1;

	mat_no_w[6][4] = imagUnity*stabRe*inv_t*inv_mu*dEicon_no_W +
		alpha*alpha + beta*beta;

	mat_no_w[7][4] = -inv_mu*mu1*t1;
// set W matrix
// of the same sizes
// 5rows x 8 cols
	// first <--> part of "2", 3 non-zeros
	mat_w[0][0] = imagUnity*stabRe*inv_t*inv_mu*dEicon_W;

	mat_w[3][0] = alpha*( -	vCoefM*gMaMa*dEicon_W);

	mat_w[4][0] = vCoefM*alpha*inv_t*dEicon_W;

	// second <--> part of "3", 2 non-zeros

	mat_w[3][1] = -imagUnity*gMaMa*dEicon_W;

	mat_w[4][1] = imagUnity*inv_t*dEicon_W;

	// third <--> "4", 4 non-zeros, a lot)

	mat_w[2][2] = xi*(-imagUnity*stabRe*inv_t*inv_mu*dEicon_W);

	mat_w[3][2] = -imagUnity*xi*vCoefL*gMaMa*
		(
		(inv_mu*mu1*t1+inv_t*t1)*dEicon_W
		);

	mat_w[4][2] = imagUnity*xi*
		(
		vCoefL*inv_t*inv_mu*mu1*t1*dEicon_W
		);

	mat_w[5][2] = imagUnity*xi*vCoefL*inv_t*dEicon_W;

	// fourth <--> part of "6", 2 non-zeros

	mat_w[3][3] = -imagUnity*stabRe*MF_Field::Pr*inv_mu*g_1MaMa*dEicon_W;

	mat_w[4][3] = imagUnity*stabRe*MF_Field::Pr*inv_t*inv_mu*dEicon_W;

	// fifth <--> "8", 3 non-zeros

	mat_w[3][4] = -vCoefM*beta*gMaMa*dEicon_W;

	mat_w[4][4] = vCoefM*beta*inv_t*dEicon_W;

	mat_w[6][4] = imagUnity*stabRe*inv_t*inv_mu*dEicon_W;

	return;
};*/

t_Vec t_StabSolver::_formRHS3D(const double& a_y, const t_Vec& a_vars){
	// TODO: check input vector dimension and
	// if it is not of 8 elems throw mega-exception)
	_setStabMatrix3D(a_y);
	t_Matrix input(1, 8);
	for (int i=0; i<8; i++){
		input[0][i] = a_vars[i];
	};
	// vars used in stability matrix 
	
	// after multiplication we have a rhs vector - matrix 1x8
	t_Matrix output = _stab_matrix.mul(input); 
	return output[0];
};

t_Matrix t_StabSolver::_getAsymptotics3D(const t_WaveChars& a_waveChars){
	t_Matrix initial_vectors(4,8);
	t_SqMatrix b_coef(4);
	t_Vec lambda(4,0.0);
	const double& y_e = _profStab.y.back();
	// TODO: function for simplified asymp: u=1.0, u'=0, u''=0, ... ?
	_setStabMatrix3D(y_e);	

	// TODO: what is all this about?

	b_coef[0][0]=_stab_matrix[0][1];
	b_coef[1][0]=_stab_matrix[3][1];
	b_coef[2][0]=_stab_matrix[4][1];
	for (int i=1; i<3; i++){
		b_coef[i][1] =  _stab_matrix[i+2][1]*_stab_matrix[1][3]+
						_stab_matrix[i+2][2]*_stab_matrix[2][3]+
						_stab_matrix[i+2][5]*_stab_matrix[5][3]+
						_stab_matrix[i+2][7]*_stab_matrix[7][3];
		b_coef[i][2] = _stab_matrix[i+2][5];
		b_coef[i][3] = _stab_matrix[i+2][7];
	};
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
		t_CompVal l2 = std::pow(lambda[i],2);
		t_CompVal denom = _stab_matrix[0][1] - l2;
		b_coef[0][i] = (
							(l2 - _stab_matrix[4][5])*_stab_matrix[3][1]+
							_stab_matrix[4][1]*_stab_matrix[3][5]
					   )/denom;

		b_coef[1][i] = _stab_matrix[4][5] - l2;
		b_coef[2][i] = -_stab_matrix[3][5];
		b_coef[3][i] = (
							_stab_matrix[3][5]*_stab_matrix[4][7]+
							(l2 - _stab_matrix[4][5])*_stab_matrix[3][7]
						)/denom;
	};

	b_coef[0][3]=0.0;
	b_coef[1][3]=0.0;
	b_coef[2][3]=0.0;
	b_coef[3][3]=1.0;
	
	for (int i=0; i<4; i++){	
		initial_vectors[i][0] = b_coef[0][i];
		initial_vectors[i][1] = lambda[i]*b_coef[0][i];
		initial_vectors[i][2] = (
									_stab_matrix[0][2]*b_coef[0][i]+
									_stab_matrix[3][2]*b_coef[1][i]+
									_stab_matrix[4][2]*b_coef[2][i]+
									_stab_matrix[6][2]*b_coef[3][i]
								)/lambda[i];

		initial_vectors[i][3] = b_coef[1][i];
		initial_vectors[i][4] = b_coef[2][i];
		initial_vectors[i][5] = lambda[i]*b_coef[2][i];
		initial_vectors[i][6] = b_coef[3][i];
		initial_vectors[i][7] = (
									_stab_matrix[3][7]*b_coef[1][i]+
									_stab_matrix[4][7]*b_coef[2][i]+
									_stab_matrix[6][7]*b_coef[3][i]
								)/lambda[i];
	}
	return initial_vectors;
}


void t_StabSolver::set3DContext(const int& i_ind, const int& k_ind, const int& a_nnodesStab){
	_math_solver.set3DContext();
	_math_solver.resizeGrid(a_nnodesStab);

	t_ProfileNS profNS(_rFldNS);
	profNS.setProfiles(i_ind, k_ind);
	
	if (_profStab.size()!=a_nnodesStab){
		_profStab.resize(a_nnodesStab);
	}
	_profStab.setProfiles(profNS);

	if (_stab_matrix.nCols!=8){
		_stab_matrix.resize(8);
	}

	for (int j=0; j<a_nnodesStab; j++){
		_math_solver.varRange[j] = _profStab.y[a_nnodesStab-1-j];
	}

	return;
}

/*   COMPUTATION OF GROUP VELOCITY (VA,VB)=(DW/DA,DW/DB) */        
void t_StabSolver::_calcGroupVelocity(t_WaveChars &a_wave_chars){
	// empiric constant - tolerance
	const double dd = 1.0e-6;
	// ensure that we are in eigen 
	adjustLocal(a_wave_chars, W_MODE);

	t_WaveChars left_chars = a_wave_chars;
	t_WaveChars rght_chars = a_wave_chars;
	// dalpha
	left_chars.a-=dd;
	rght_chars.a+=dd;
    adjustLocal(left_chars, W_MODE);
	adjustLocal(rght_chars, W_MODE);
	a_wave_chars.vga = 0.5*(rght_chars.w - left_chars.w)/dd;

	// reset left, right
	left_chars = a_wave_chars;
	rght_chars = a_wave_chars;

	left_chars.b-=dd;
	rght_chars.b+=dd;
    adjustLocal(left_chars, W_MODE);
	adjustLocal(rght_chars, W_MODE);
	a_wave_chars.vgb = 0.5*(rght_chars.w - left_chars.w)/dd;
}

t_WaveChars t_StabSolver::_getStationaryMaxInstabTime(const t_WaveChars& initial_guess){
	t_WaveChars cur_wave = initial_guess;
	t_WaveChars def_wave = initial_guess;
	t_Complex gv;
	double da, db;     

	do{
		adjustLocal(def_wave, W_MODE);   
		_calcGroupVelocity(def_wave);
		db = 0.01*def_wave.a.real();
		double coef = def_wave.vgb.real()/def_wave.vga.real();
		da = -coef*db;

		cur_wave.a = def_wave.a + da;
		cur_wave.b = def_wave.b + db;
		cur_wave.w = def_wave.w;
		adjustLocal(cur_wave, W_MODE);

		if (cur_wave.w.imag()<def_wave.w.imag()){
			db = -db;
			da = -coef*db;
			cur_wave.a = def_wave.a + da;
			cur_wave.b = def_wave.b + db;
			cur_wave.w = def_wave.w;
			adjustLocal(cur_wave, W_MODE);		
		}	
	} while(def_wave.w.imag()<cur_wave.w.imag());
	std::cout<<"STAT:A="<<def_wave.a.real()<<"; B="<<def_wave.b.real()<<"\n----:W="
		<<def_wave.w<<std::endl;
	return def_wave;
};

t_WaveChars t_StabSolver::_getMaxInstabTime(const t_WaveChars &init_guess){
	// empiric half percent per iteration
	double dar = 0.005*init_guess.a.real();
	t_WaveChars base_wave = _getStationaryMaxInstabTime(init_guess);
	_calcGroupVelocity(base_wave);
	t_WaveChars next_wave = base_wave;
	next_wave.a+=dar;
	next_wave.w+=base_wave.vga*dar;
	next_wave = _getStationaryMaxInstabTime(next_wave);
	// do we move in proper direction ?
	if (next_wave.w.imag()<base_wave.w.imag()){
		dar = -dar;
	};
	
	do{
		base_wave = next_wave;
		_calcGroupVelocity(base_wave);
		next_wave.a+=dar;
		next_wave.w+=base_wave.vga*dar;
		next_wave = _getStationaryMaxInstabTime(next_wave);
	} while (base_wave.w.imag()>next_wave.w.imag());
	return base_wave;
};
// input: mode defines the value to be adjusted
// see stabsolver.h for mode enum definition
void t_StabSolver::adjustLocal(t_WaveChars &a_wave_chars, t_StabSolver::t_MODE a_mode){
	// parameters
	const t_Complex d_arg = 1.0e-5;
	const int max_iters = 50;
	const double tol = 1.0e-4;
	// ensure that there is a base residual
	solve(a_wave_chars);
	t_WaveChars backup = a_wave_chars;
	t_Complex *choose_arg;
	if (a_mode==A_MODE){
		choose_arg = &(a_wave_chars.a);
	}else{
		if (a_mode==B_MODE){
			choose_arg = &(a_wave_chars.b);
		}else{
			choose_arg = &(a_wave_chars.w);
		};
	};
	t_Complex& arg = *choose_arg;
	t_Complex arg_base = arg;
	t_Complex res_base = a_wave_chars.resid;
	t_Complex deriv;
	for (int i=0; i<max_iters; i++){
		// check if we are converged
		if (std::norm(a_wave_chars.resid)<tol){
			return;
		};
		// newton iteration
		arg+=d_arg;
		solve(a_wave_chars);
		deriv = d_arg/(a_wave_chars.resid - res_base);
		arg = arg_base - res_base*deriv;
		solve(a_wave_chars);
		arg_base = arg;
		res_base = a_wave_chars.resid;
		std::cout<<"Cur SI"<<res_base<<std::endl;
	};
	std::cerr<<"Local Search Error: no convergence"<<std::endl;
	a_wave_chars = backup;
};

void t_StabSolver::setInitWaves(const std::vector<t_WaveChars>& a_inits){
	this->_initWaves = a_inits;
};

t_WaveChars t_StabSolver::getMaxWave(const int& i_ind, const int& k_ind, 
									  const int& a_nnodesStab, const std::vector<t_WaveChars>& a_inits){
	set3DContext(i_ind, k_ind, a_nnodesStab);
	setInitWaves(a_inits);
	t_WaveChars res_wave;
	//ensure
	res_wave.w = 0.0;
	t_WaveChars cur_wave;
	for (int i=0; i<_initWaves.size(); i++){
		cur_wave = _getMaxInstabTime(_initWaves[i]);
		if (cur_wave.w.imag()>res_wave.w.imag()){
			res_wave = cur_wave;
		}
	}
	return res_wave;
}
