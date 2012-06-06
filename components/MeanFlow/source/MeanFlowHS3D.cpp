#include "MeanFlow.h"
t_MFHSFLOW3D::t_MFHSFLOW3D(wxString configfile):
_params(configfile), t_MeanFlow(), t_Component(configfile){
	_allocate(_params.Nx, _params.Ny, _params.Nz);
	_init();
	// init params map etc - to bind to GUI
}

void t_MFHSFLOW3D::_init(){
	const int& Nx = _params.Nx;
	const int& Ny = _params.Ny;
	const int& Nz = _params.Nz;
	FILE* fld_file = fopen(_params.mf_bin_path.ToAscii(),"rb");
	//reverse order!! k,j,i
	for(int k=0; k<Nz; k++)	
		for (int j=0; j<Ny; j++)
			for(int i=0; i<Nx; i++){
				t_Rec& ptr = _fld[i][j][k];
				fread(&ptr.x,sizeof(double),1,fld_file);
				fread(&ptr.y,sizeof(double),1,fld_file);
				fread(&ptr.z,sizeof(double),1,fld_file);
			}
			for(int k=0; k<Nz; k++)	
				for (int j=0; j<Ny; j++)
					for(int i=0; i<Nx; i++){
						t_Rec& ptr = _fld[i][j][k];
						fread(&ptr.u,sizeof(double),1,fld_file);
						fread(&ptr.v,sizeof(double),1,fld_file);
						fread(&ptr.w,sizeof(double),1,fld_file);
						fread(&ptr.p,sizeof(double),1,fld_file);
						fread(&ptr.t,sizeof(double),1,fld_file);
						double gmama = _params.Gamma*_params.Mach*_params.Mach;
						ptr.r=gmama*ptr.p/ptr.t;
					};

			fclose(fld_file);
};

void t_MFHSFLOW3D::_init_params_grps(){
	_add_params_group(_T("default"), _params);
}

const t_MFParams& t_MFHSFLOW3D::base_params() const{
	return _params;
};
