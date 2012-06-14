#include "MeanFlow.h"
#include "common_data.h"

t_MFHSFLOW2D::t_MFHSFLOW2D():
t_MeanFlow(), t_Component(wxEmptyString, MF_HSFLOW2D_NAME){};

t_MFHSFLOW2D::t_MFHSFLOW2D(const wxString& configfile):
t_MeanFlow(), t_Component(configfile, MF_HSFLOW2D_NAME){
	_init(configfile);
};

void t_MFHSFLOW2D::initialize( const wxString& configfile ){
	_init(configfile);
};

void t_MFHSFLOW2D::_init(const wxString& configfile){
	_init_params_grps();
	_paramsFileName = configfile;
	_params.load_via_params(configfile);
	const int& Nx = _params.Nx;
	const int& Ny = _params.Ny;
	const int& Nz = _params.Nz;
	_allocate(_params.Nx, _params.Ny, _params.Nz);
	FILE* fld_file = fopen(_params.mf_bin_path.ToAscii(),"rb");
	double gmama = _params.Gamma*_params.Mach*_params.Mach;
	// in 2D fields the packing is by column
	for (int j=0; j<Ny; j++){
		for (int i=0; i<Nx; i++){
			t_Rec& rRec = _fld[i][j][0];
			fread(&rRec.x,sizeof(double),1,fld_file);
			fread(&rRec.y,sizeof(double),1,fld_file);
		}
	}
	for (int j=0; j<Ny; j++){
		for (int i=0; i<Nx; i++){
			t_Rec& rBaseRec = _fld[i][j][0];
			fread(&rBaseRec.u,sizeof(double),1,fld_file);
			fread(&rBaseRec.v,sizeof(double),1,fld_file);
			fread(&rBaseRec.p,sizeof(double),1,fld_file);
			fread(&rBaseRec.t,sizeof(double),1,fld_file);
			rBaseRec.r=gmama*rBaseRec.p/rBaseRec.t;
			for (int k=0; k<Nz; k++){
				t_Rec& rRec = _fld[i][j][k];
				if (_params.MFSym==t_MFParamsHS2D::t_AxeSym::AxeSym){
					double psi = M_PI/double(Nz-1)*k;
					rRec.x = rBaseRec.x;
					rRec.y = rBaseRec.y*cos(psi);
					rRec.z = rBaseRec.y*sin(psi);
					rRec.u = rBaseRec.u;
					rRec.v = rBaseRec.v*cos(psi);
					rRec.w = rBaseRec.v*sin(psi);
					rRec.p = rBaseRec.p;
					rRec.t = rBaseRec.t;
					rRec.r = rBaseRec.r;
				}else{
					// ???
					double z_span = 1.0;
					rRec = rBaseRec;
					rRec.z = (-0.5+double(k)/double(Nz-1))*_params.ZSpan; // z from -0.5*ZSpan to 0.5*ZSpan
					rRec.w = 0.0;
				};
			}
		}
	}
	fclose(fld_file);
};

void t_MFHSFLOW2D::_init_params_grps(){
	_mapParamsGrps.clear();
	_add_params_group(_T("default"), _params);
}

const t_MFParams& t_MFHSFLOW2D::base_params() const{
	return _params;
}
