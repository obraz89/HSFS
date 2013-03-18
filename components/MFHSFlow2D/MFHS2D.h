#ifndef __MFHS2D
#define __MFHS2D

#include "PluginBase.h"
#include "MFBlockBase.h"

#include "MFHS2D_params.h"


#include "io_helpers.h"


namespace mf{

	class t_MFHSFLOW2D: public mf::t_Block{
	private:

		wxString _mf_bin_path;
		t_HSFlowParams2D _base_params;

		void _init();

	public:

		void init(const hsstab::TPlugin& g_plug);

		const mf::t_FldParams& get_mf_params() const;
		const mf::t_HSFlowParams2D& get_params() const;
	}; 

}			//  ~namespace mf
#endif		//  __MFHS2D
