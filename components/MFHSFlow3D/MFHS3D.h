#ifndef __MFHS3D
#define __MFHS3D

#include "PluginBase.h"
#include "MFHSDomain.h"

#include "MFHS3D_params.h"

namespace mfhs{

	class t_Block3D: public mfhs::t_Block{
	public:

		t_Block3D(const t_Domain& domain);
		void init(int nx, int ny, int nz, wxString mf_bin_path);

	};

//-----------------------------------------------------------------------------

	class t_Domain3D: public mfhs::t_Domain, public hsstab::TPlugPhysPart{

	private:

		t_Block3D _blk;

		wxString _mf_bin_path;
		mf::t_FldParams _base_params;

		void _init();

	public:

		// implement t_Domain
		t_Domain3D();
		const mf::t_FldParams& get_mf_params() const;

		// implement TPlugPhysPart
		void init(const hsstab::TPlugin& g_plug);
		const t_Block& get_blk() const;
	}; 


	}		// ~namespace mf
#endif //MFHS3D
