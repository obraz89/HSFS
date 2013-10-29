#ifndef __MFCGNS2D
#define __MFCGNS2D

#include "PluginBase.h"
#include "MFDomainBase.h"

#include "cgns_structs.h"
#include "MF_CGNS2D_params.h"


#include "io_helpers.h"


namespace mf{

	class t_MFCGNS2D: public mf::cg::TDomain{
	protected:

		wxString _grd_bin_path, _fld_bin_path;
		t_CGNS2DParams _base_params;

		void _init();

		bool _parseGhostData2DfromCGNS( cg::TcgnsContext& ctx );

		bool _parseBCData2DfromCGNS( cg::TcgnsContext& ctx );

		bool _doLoadGrid2D_cgns( const wxString& gridFN );

		// geom staff

		void get_k_range(int iZone, int& ks, int& ke) const;


	public:

// t_DomainBase interface realization

		void init(const hsstab::TPlugin& g_plug);

		const mf::t_FldParams& get_mf_params() const;

		virtual mf::t_Rec interpolate_to_point(const t_GeomPoint& point) const;
		// TODO: implement for optimization
		//virtual void interpolate_to_point(const t_GeomPoint& point, t_Rec& rec) const=0;

		void get_rec(const mf::cg::TZone& blk, int i, int j, int k, mf::t_Rec& rec) const;
		t_Rec get_rec(const t_GeomPoint& xyz) const;

		// are we still inside domain?
		// TODO; move to CHNS-shared when implemented
		bool is_point_inside(const t_GeomPoint& xyz) const;

		// from cg::TDomain
		bool loadGrid(const wxString& gridFN );

		const mf::t_CGNS2DParams& get_params() const;
}; 



}			//  ~namespace mf
#endif		//  __MFCGNS2D
