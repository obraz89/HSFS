#ifndef __MFCGNS2D
#define __MFCGNS2D

#include "PluginBase.h"
#include "MFBlockBase.h"

#include "MF_CGNS2D_params.h"


#include "io_helpers.h"


namespace mf{

	class t_MFCGNS2D: public mf::cg::TDomain{
	private:

		wxString _grd_bin_path, _fld_bin_path;
		t_CGNS2DParams _base_params;

		void _init();

		bool _parseGhostData2DfromCGNS( TcgnsContext& ctx );

		bool _parseBCData2DfromCGNS( TcgnsContext& ctx );

		bool _doLoadGrid2D_cgns( const wxString& gridFN );


	public:

		void init(const hsstab::TPlugin& g_plug);
		bool loadGrid(const wxString& gridFN );

		const mf::t_FldParams& get_mf_params() const;
		const mf::t_CGNS2DParams& get_params() const;
}; 



}			//  ~namespace mf
#endif		//  __MFCGNS2D
