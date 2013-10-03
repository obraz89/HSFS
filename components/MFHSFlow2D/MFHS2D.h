#ifndef __MFHS2D
#define __MFHS2D

#include "PluginBase.h"
#include "MFHSDomain.h"

#include "MFHS2D_params.h"

namespace mfhs{

	class t_Block2D: public mfhs::t_Block{
	public:

		t_Block2D(const t_Domain& domain);
		void init(int nx, int ny, int nz, wxString mf_bin_path, const t_HSFlowParams2D& params);

	};

//-----------------------------------------------------------------------------

	class t_Domain2D: public mfhs::t_Domain, public hsstab::TPlugPhysPart{

	private:

		t_Block2D _blk;

		wxString _mf_bin_path;
		t_HSFlowParams2D _base_params;

		void _init();

	public:

		// implement t_Domain
		t_Domain2D();
		const mf::t_FldParams& get_mf_params() const;

		// implement TPlugPhysPart
		void init(const hsstab::TPlugin& g_plug);
		const t_Block& get_blk() const;
		const t_HSFlowParams2D& get_params() const;
	}; 

}			//  ~namespace mf
#endif		//  __MFHS2D
