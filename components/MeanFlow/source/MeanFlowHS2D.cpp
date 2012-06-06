#include "MeanFlow.h"
t_MFHSFLOW2D::t_MFHSFLOW2D(wxString configfile):
_params(configfile), t_MeanFlow(){
	_allocate(_params.Nx, _params.Ny, _params.Nz);
	_init();
};

void t_MFHSFLOW2D::_init(){
	const int& Nx = _params.Nx;
	const int& Ny = _params.Ny;
	const int& Nz = _params.Nz;
	FILE* fld_file = fopen(_params.mf_bin_path.ToAscii(),"rb");
	double gmama = _params.Gamma*_params.Mach*_params.Mach;
	// in 2D fields the packing is by column
	t_Rec base_rec2D;
	for (int i=0; i<Nx; i++){
		for (int j=0; j<Ny; j++){
			fread(&base_rec2D.x,sizeof(double),1,fld_file);
			fread(&base_rec2D.y,sizeof(double),1,fld_file);
			fread(&base_rec2D.u,sizeof(double),1,fld_file);
			fread(&base_rec2D.v,sizeof(double),1,fld_file);
			fread(&base_rec2D.p,sizeof(double),1,fld_file);
			fread(&base_rec2D.t,sizeof(double),1,fld_file);
			base_rec2D.r=gmama*base_rec2D.p/base_rec2D.t;
			for (int k=0; k<Nz; k++){
				t_Rec& rRec = _fld[i][j][k];
				if (_params.MFSym==t_MFParamsHS2D::t_AxeSym::AxeSym){
					double psi = 2*M_PI/double(Nz-1)*k;
					rRec.x = base_rec2D.x;
					rRec.y = base_rec2D.y*cos(psi);
					rRec.z = base_rec2D.y*sin(psi);
					rRec.u = base_rec2D.u;
					rRec.v = base_rec2D.v*cos(psi);
					rRec.w = base_rec2D.v*sin(psi);
					rRec.p = base_rec2D.p;
					rRec.t = base_rec2D.t;
					rRec.r = base_rec2D.r;
				}else{
					// ???
					double z_span = 1.0;
					rRec = base_rec2D;
					rRec.z = k-0.5; // z from -0.5 to 0.5
					rRec.w = 0.0;
				};
			}
		}
	}
	fclose(fld_file);
};

const t_MFParams& t_MFHSFLOW2D::base_params() const{
	return _params;
}
