#include "stdafx.h"

#include "MF_CGNS2D.h"

#include <cgnslib.h>
#include "cgns_structs.h"

using namespace mf;
using namespace mf::cg;

void t_MFCGNS2D::get_k_range(int iZone, int& ks, int& ke) const{
	ks = 1;
	ke = _base_params.Nz;
}

bool t_MFCGNS2D::loadGrid( const wxString& gridFN ){
	bool ok = false;

	wxLogMessage( _("* Grid: Loading from '%s' & processing..."), gridFN.c_str() );

	ok = _doLoadGrid2D_cgns( gridFN );
	if( ! ok ) return false;

	//
	// Distribute blocks through CPUs
	// [HSFlow Legacy]
	/*
	{
		if( G_Domain.nZones < G_State.mpiNProcs )
		{
			wxLogError(_("Number of domain zones is less than MPI procs"));
			return false;
		}

		const int nAvg   = G_Domain.nZones / G_State.mpiNProcs;
		const int nResdl = G_Domain.nZones % G_State.mpiNProcs;

		for( int r = 0; r < G_State.mpiNProcs; ++r )
		{
			const int bs = r * nAvg + ((r<nResdl) ?r :nResdl);
			const int be = bs + nAvg + ((r<nResdl) ?1 :0) - 1;

			wxLogMessage( _("* MPI rank %d owns:"), r );
			for( int b = bs; b <= be; ++b )
			{
				TZone& zne = G_Domain.Zones[b];

				wxLogMessage( _("    block %d: '%s' %dx%dx%d (with ghosts)"),
					b+1, wxString::FromAscii(zne.szName).c_str(),
					zne.nx, zne.ny, zne.nz
					);
			}

			if( r == G_State.mpiRank )
			{
				G_Domain.bs = bs;
				G_Domain.be = be;
			}
		}
	}
	*/

	// One core distribution
	this->bs = 0;
	this->be = nZones - 1;

	// [HSFlow legacy]
	//G_Plugins.get_mesh_caps().calcHalfnodes();

	return true;
}
/**
**/  
//-----------------------------------------------------------------------------


bool t_MFCGNS2D::_doLoadGrid2D_cgns( const wxString& gridFN )
{
	int res = CG_OK;
	char szName[33];  // names in CGNS file

	TcgnsContext& ctx = this->cgCtx;

	res = cg_open( gridFN.ToAscii(),CG_MODE_READ, &ctx.fileID );
	if( res != CG_OK )
	{
		wxLogError( _("Can't open grid file for reading") );
		return false;
	}
	ctx.iBase = 1; // assume only one base

	//
	// Space dimensions
	// 
	int dimCell = 0, dimPhys = 0;
	cg_base_read(ctx.fileID,ctx.iBase,  szName,&dimCell,&dimPhys);
	if( dimCell != 2 )
	{
		wxLogError( _("The grid is not for 2D problems") );
		return false;
	}

	//
	// Number of zones (aka blocks)
	//
	nZones = 0;  cg_nzones(ctx.fileID, ctx.iBase, &nZones);
	if( nZones < 1 )
	{
		wxLogError( _("Domain blocks are not found") );
		return false;
	}

	//
	// Allocate memory for whole computational domain
	//
	Zones = new TZone[ nZones ];
	ctx.cgZones = new TcgnsZone[ nZones ];
	if( ! Zones )  wxLogFatalError(_("Not enough memory for domain"));


	//
	// Get zones (blocks) sizes and names
	//
	for( int iZone = 1;  iZone <= nZones;  ++iZone )
	{
		CG_ZoneType_t type;  cg_zone_type(ctx.fileID,ctx.iBase,iZone, &type);
		if( type != CG_Structured )
		{
			wxLogError( _("Only structured grids are supported") );
			return false;
		}

		TZone& blk = Zones[iZone -1];

		// Get zone (block) name & size
		cgsize_t isize[6]; // NVertexI, NVertexJ, NCellI, NCellJ, NBoundVertexI, NBoundVertexJ
		res = cg_zone_read(ctx.fileID,ctx.iBase,iZone,  blk.szName,isize);
		if( res != CG_OK )
		{
			wxLogError( _("Can't read zone #%d info"), iZone );
			return false;
		}

		blk.nx = isize[0];  blk.ny = isize[1];   blk.nz = 1;

		// Indices of real (not ghost) nodes
		blk.is = 1;  blk.ie = blk.nx;
		blk.js = 1;  blk.je = blk.ny;
		blk.ks = 1;  blk.ke = 1;

		ctx.map_zoneName_id[blk.szName] = iZone;

	}


	//
	// Connectivity info
	// NB: blk.{nx, ny, nz} will be updated
	if( ! _parseGhostData2DfromCGNS(ctx) )  return false;

	//
	// Boundary conditions info
	//
	if( ! _parseBCData2DfromCGNS(ctx) )  return false;


	//
	// Read grid coordinates for all zones
	// (to simplify mesh data transfer to ghosts)
	//
	for( int b = 0;  b < nZones;  ++b )
	{
		int iZone = b + 1;
		TZone& blk = Zones[b];

		// Original block size without ghosts
		const int
			nx0 = blk.ie - blk.is + 1,
			ny0 = blk.je - blk.js + 1;

		// Zone index ranges
		cgsize_t irmin[2] = {1, 1};
		cgsize_t irmax[2] = {nx0, ny0};

		double* x = new double[nx0*ny0];
		res = cg_coord_read(ctx.fileID,ctx.iBase,iZone,"CoordinateX",CG_RealDouble,irmin,irmax, x);

		double* y = new double[nx0*ny0];
		res |= cg_coord_read(ctx.fileID,ctx.iBase,iZone,"CoordinateY",CG_RealDouble,irmin,irmax, y);

		if( res != CG_OK )  return false;

		//
		// Grid step in computational (curvilinear) coordinates
		//
		blk.grd.dx = 1./(blk.nx - 1);
		blk.grd.dy = 1./(blk.ny - 1);
		blk.grd.dz = 0.;

		//
		// Packed grid data storage
		//
		blk.grd.c2d = new TgridCell2D[blk.nx*blk.ny];
		if( ! blk.grd.c2d )  wxLogFatalError(_("Not enough memory for mesh"));

		// Pack grid coordinates
		for( int j = blk.js; j <= blk.je; ++j )
		for( int i = blk.is; i <= blk.ie; ++i )
		{
			Tmtr2D& node = blk.cell(i,j).ij;

			int n = (i - blk.is) + (j - blk.js)*nx0;
			node.x = x[n];
			node.y = y[n];
			//node.Rw = HUGE_VAL;
		}

		delete[] x, y;
	} // loop through zones


	//
	// Transfer grid data to ghost nodes
	// 
	for( int b = 0; b < nZones; ++b )
	{
		TZone& zne = Zones[b];

		for( int j = 1; j <= zne.ny; ++j )
		for( int i = 1; i <= zne.nx; ++i )
		{
			if( i >= zne.is && i <= zne.ie &&
				j >= zne.js && j <= zne.je    )
				continue; // myself

			const int locInd = zne.absIdx(i,j) - 1;
			const int globInd = zne.globIndices[ locInd ];
			if( globInd == -1 ) // unmapped
				continue;

			// Detect donor zone
			int nDnrZne = -1;
			for( nDnrZne = nZones - 1; nDnrZne >= 0; --nDnrZne )
			{
				if( globInd >= Zones[ nDnrZne ].nGlobStart )
					break;
			}

			TZone& dnrZne = Zones[ nDnrZne ];
			int dnrLocInd = dnrZne.absIdxFromGlobInd( globInd ) - 1; // 0-based

			zne.grd.c2d[ locInd ] = dnrZne.grd.c2d[ dnrLocInd ];
		}
	}

	cg_close( ctx.fileID );

	return true;
}
//-----------------------------------------------------------------------------


/**
 *  Read connectivity data from the opened CGNS file
 *  and initialize ghost indices 
 *  
 *  @param[in] ctx - context of the opened CGNS file
 *  
 *  @return true if succeeded, false if failed
**/  
bool t_MFCGNS2D::_parseGhostData2DfromCGNS( TcgnsContext& ctx )
{
	// Number of ghost nodes (half of stencil size)
	int ghostI = 2, ghostJ = 2;
	// TODO:
	//G_Plugins.get_discret_caps().getParI(11/*nsx*/, ghostI);  ghostI /= 2;
	//G_Plugins.get_discret_caps().getParI(24/*nsy*/, ghostJ);  ghostJ /= 2;

//
// (1) Read connectivity data from CGNS file
//
for( int iZone = 1;  iZone <= nZones;  ++iZone )
{
	TZone& zne = Zones[iZone-1];
	TcgnsZone& cgZne = ctx.cgZones[iZone-1];

	int n1to1 = 0;  cg_n1to1(ctx.fileID,ctx.iBase,iZone, &n1to1);
	if( n1to1 < 1 && nZones>1 )
	{
		wxLogError(_("Missing inter-block 1-to-1 connectivity info"));
		return false;
	}

	// Original block size without ghosts
	const int
		nx0 = zne.ie - zne.is + 1,
		ny0 = zne.je - zne.js + 1;

	// Loop through connectivity patches
	for( int iConn = 1; iConn <= n1to1; ++iConn )
	{
		cgsize_t idxRange[4];   // Imin, Jmin, Imax, Jmax

		cgsize_t idxRangeDonor[4];

		char szName[33], szDonor[33];
		int iTsh[2]; // short-hand notation of transform matrix
		bool res = cg_1to1_read( ctx.fileID,ctx.iBase,iZone,iConn,
			szName, szDonor, idxRange, idxRangeDonor, iTsh
		);
		if( res != CG_OK )
		{
			wxLogError( _("Can't read connection patch '%s'(#%d) of zone '%s'(#%d)"),
				wxString::FromAscii(szName).c_str(), iConn,
				wxString::FromAscii(zne.szName).c_str(), iZone
			);
			return false;
		}

		const long
			&is = idxRange[0],  &js = idxRange[1],
			&ie = idxRange[2],  &je = idxRange[3];

		const long
			&ids = idxRangeDonor[0],  &jds = idxRangeDonor[1],
			&ide = idxRangeDonor[2],  &jde = idxRangeDonor[3];

		// Check that whole faces connected
		bool ok = true;
		ok &= (is==1 && ie==1) || (is==nx0 && ie==nx0) || (is==1 && ie==nx0);
		ok &= (js==1 && je==1) || (js==ny0 && je==ny0) || (js==1 && je==ny0);
		if( ! ok )
		{
			wxLogError(
				_("Connection patch '%s'(#%d) of zone '%s'(#%d) don't cover whole block face"),
				wxString::FromAscii(szName).c_str(), iConn,
				wxString::FromAscii(zne.szName).c_str(), iZone
			);
			return false;
		}

		//
		// Detect connection face of the current block
		// and 
		// increase the block size by ghosts
		// 
		TZoneFacePos posFace = faceNone;  // block face

		if( is==ie ) // Xmin or Xmax block face
		{
			zne.nx += ghostI;
			if( is==1 )
			{
				posFace = faceXmin;
				zne.is += ghostI;
				zne.ie += ghostI;
			}
			else
				posFace = faceXmax;
		}
		else if( js==je )  // Ymin or Ymax block face
		{
			zne.ny += ghostJ;
			if( js==1 )
			{
				posFace = faceYmin;
				zne.js += ghostJ;
				zne.je += ghostJ;
			}
			else
				posFace = faceYmax;
		}

		zne.Faces[posFace].bcType = bcAbutted;

		TcgnsZone::TFace& face = cgZne.Faces[ posFace ];
		face.nDnrZne = ctx.map_zoneName_id[szDonor] - 1;  // 0-based

		face.dnr_is0 = ids;
		face.dnr_js0 = jds;

		// Transform matrix between node indexes
		// of the current and donor blocks
		// 
		// CGNS SIDS documentation,
		// section 8.2 1-to-1 Interface Connectivity Structure Definition
		// 
		// iTsh = [a, b, c] =>
		//	   | sgn(a)del(a,1) sgn(b)del(b,1) sgn(c)del(c,1) |
		// T = | sgn(a)del(a,2) sgn(b)del(b,2) sgn(c)del(c,2) |,
		//	   | sgn(a)del(a,3) sgn(b)del(b,3) sgn(c)del(c,3) |
		struct F {
			static int sgn(int x){ return (x<0)?-1:+1; }
			static int del(int x, int y){ return (abs(x)==abs(y))?+1:0; }
		};

		face.matTrans[0] = F::sgn(iTsh[0]) * F::del(iTsh[0],1);
		face.matTrans[1] = F::sgn(iTsh[1]) * F::del(iTsh[1],1);

		face.matTrans[2] = F::sgn(iTsh[0]) * F::del(iTsh[0],2);
		face.matTrans[3] = F::sgn(iTsh[1]) * F::del(iTsh[1],2);
	} // loop through connection patches
} // loop through zones


//
// (2) Initial global indices in all zones
// 
int gIdx = 0;
for( int z = 0; z < nZones; ++z )
{
	TZone& zne = Zones[z];

	zne.nGlobStart = gIdx;

	zne.globIndices = new int[ zne.nx * zne.ny ];
	memset( zne.globIndices, 0xFF, zne.nx * zne.ny * sizeof(int) );  // -1

	// Fill-in global indices of real nodes
	for( int j=zne.js; j<=zne.je; ++j )
	for( int i=zne.is; i<=zne.ie; ++i )
	{
		const int locInd = zne.absIdx(i,j) - 1;  // 0-based
		zne.globIndices[locInd] = zne.globRealInd(i,j);
	}

	// Real nodes count (without ghosts)
	const int
		nx0 = zne.ie - zne.is + 1,
		ny0 = zne.je - zne.js + 1;

	gIdx += nx0*ny0;
}


//
// (3) Parse connectivity data and fill-in ghost nodes indices
//
for( int b = 0; b < nZones; ++b )
{
	TZone& zne = Zones[b];

	// Loop through faces
	for( int f = 0; f < 4; ++f )
	{
		if( zne.Faces[f].bcType != bcAbutted ) continue;

		// Start & End indexes of ghost layer:
		int igs, ige, jgs, jge;

		// Start indexes of the zone face:
		int ifs, jfs;

		switch( f )
		{
		case faceXmin:
			ifs = zne.is;      jfs = zne.js;
			igs = 1;           jgs = zne.js;
			ige = zne.is - 1;  jge = zne.je;
			break;
		case faceXmax:
			ifs = zne.ie;      jfs = zne.js;
			igs = zne.ie + 1;  jgs = zne.js;
			ige = zne.nx;      jge = zne.je;
			break;
		case faceYmin:
			ifs = zne.is;      jfs = zne.js;
			igs = zne.is;      jgs = 1;      
			ige = zne.ie;      jge = zne.js - 1;
			break;
		case faceYmax:
			ifs = zne.is;      jfs = zne.je;
			igs = zne.is;      jgs = zne.je + 1;
			ige = zne.ie;	   jge = zne.ny;
			break;
		}

		TcgnsZone::TFace& face = ctx.cgZones[b].Faces[f];
		TZone& zneDonor = Zones[ face.nDnrZne ];

		// Donor block face starting index with ghosts
		int ids = face.dnr_is0 + zneDonor.is - 1;
		int jds = face.dnr_js0 + zneDonor.js - 1;

		for( int j=jgs; j<=jge; ++j )
		for( int i=igs; i<=ige; ++i )
		{
			// Donor block indexes
			int id = (i-ifs)*face.matTrans[0] + (j-jfs)*face.matTrans[1] + ids;
			int jd = (i-ifs)*face.matTrans[2] + (j-jfs)*face.matTrans[3] + jds;
			
			if( id > zneDonor.ie || id < zneDonor.is ||
				jd > zneDonor.je || jd < zneDonor.js
			)
			{
				wxLogError(
					_("Ghost node (%d,%d) of zone '%s' is mapped \
to nonexistent node (%d,%d) of zone '%s'. Check indices orientation!"),
					i,  j, wxString::FromAscii(zne.szName).c_str(),
					id, jd, wxString::FromAscii(zneDonor.szName).c_str()
				);
				return false;
			}

			const int locInd = zne.absIdx(i,j) - 1;  // 0-based
			zne.globIndices[locInd] = zneDonor.globRealInd(id,jd);
		}
	}  // loop through faces

} // loop through blocks


//
// (4) Fill-in remaining ghost node indices around corners
// 
for( int b = 0; b < nZones; ++b )
{
	TZone& zne = Zones[b];

	// Loop through faces
	for( int f = 0; f < 4; ++f )
	{
		if( zne.Faces[f].bcType != bcAbutted ) continue;

		// Start & End indexes of ghost corner:
		int igs, ige, jgs, jge;

		// Start indexes of the zone face:
		int ifs, jfs;

		switch( f )
		{
		case faceXmin:
			ifs = zne.is;      jfs = zne.js;
			igs = 1;           jgs = 1;
			ige = zne.is - 1;  jge = zne.ny;
			break;
		case faceXmax:
			ifs = zne.ie;      jfs = zne.js;
			igs = zne.ie + 1;  jgs = 1;
			ige = zne.nx;      jge = zne.ny;
			break;
		case faceYmin:
			ifs = zne.is;      jfs = zne.js;
			igs = 1;           jgs = 1;      
			ige = zne.nx;      jge = zne.js - 1;
			break;
		case faceYmax:
			ifs = zne.is;      jfs = zne.je;
			igs = 1;           jgs = zne.je + 1;
			ige = zne.nx;	   jge = zne.ny;
			break;
		}

		TcgnsZone::TFace& face = ctx.cgZones[b].Faces[f];

		TZone& zneDonor = Zones[ face.nDnrZne ];

		// Donor block face starting index with ghosts
		int ids = face.dnr_is0 + zneDonor.is - 1;
		int jds = face.dnr_js0 + zneDonor.js - 1;

		for( int j=jgs; j<=jge; ++j )
		for( int i=igs; i<=ige; ++i )
		{
			const int locInd = zne.absIdx(i,j) - 1;  // 0-based
			if( zne.globIndices[locInd] != -1 ) continue;

			// Donor block indexes
			int id = (i-ifs)*face.matTrans[0] + (j-jfs)*face.matTrans[1] + ids;
			int jd = (i-ifs)*face.matTrans[2] + (j-jfs)*face.matTrans[3] + jds;

			if( id > zneDonor.nx || id < 1 ||
				jd > zneDonor.ny || jd < 1    )
			{
				// ghost data unavailable
				wxLogError(_("Provided multi-block layout is not supported yet. Sorry"));
				return false;
				//continue;
			}

			int dnrLocInd = zneDonor.absIdx(id,jd) - 1;  // 0-based
			zne.globIndices[locInd] = zneDonor.globIndices[ dnrLocInd ];
		}

	}  // loop through faces
} // loop through zones

	return true;
}
//-----------------------------------------------------------------------------

// copy-paste from connectivity part of grid loading
t_ZoneNode t_MFCGNS2D::get_abutted_znode(
	const t_ZoneNode& a_znode, const int di, const int dj, const int dk) const{

		if (dk!=0){
			wxLogError(_T("Get abutted znode 2D: k-shift is not allowed (2D configuration)"));
		}

		int b = a_znode.iZone - 1;

		TZone& zne = Zones[b];

		if (a_znode.iFacePos==faceNone) 
			wxLogError(_T("In get_abutted_znode: ifacepos not specified for input znode param"));

		TZoneFacePos f = a_znode.iFacePos;

		TcgnsZone::TFace& face = cgCtx.cgZones[b].Faces[f];

		TZone& zneDonor = Zones[ face.nDnrZne ];

		// Start indexes of the zone face:
		int ifs, jfs;

		switch( f )
		{
		case faceXmin:
			ifs = zne.is;      jfs = zne.js;
			break;
		case faceXmax:
			ifs = zne.ie;      jfs = zne.js;
			break;
		case faceYmin:
			ifs = zne.is;      jfs = zne.js;
			break;
		case faceYmax:
			ifs = zne.is;      jfs = zne.je;
			break;
		}

		const int i = a_znode.iNode.i + di;
		const int j = a_znode.iNode.j + dj;
		const int k = a_znode.iNode.k + dk;

		// Donor block face starting index with ghosts
		int ids = face.dnr_is0 + zneDonor.is - 1;
		int jds = face.dnr_js0 + zneDonor.js - 1;

		// Donor block indexes
		int id = (i-ifs)*face.matTrans[0] + (j-jfs)*face.matTrans[1] + ids;
		int jd = (i-ifs)*face.matTrans[2] + (j-jfs)*face.matTrans[3] + jds;

		if( id > zneDonor.nx || id < 1 ||
			jd > zneDonor.ny || jd < 1    )
		{
			wxLogError(_("Get abutted znode 2D : Failed to get corrrect abutted znode index"));
		}

		t_ZoneNode ret;

		// 1-based zone id
		ret.iZone = face.nDnrZne + 1;

		// k is unchanged
		ret.iNode = t_BlkInd(id, jd, a_znode.iNode.k);

		// face position unknown after shift
		ret.iFacePos = faceNone;

		return ret;
};


bool t_MFCGNS2D::_parseBCData2DfromCGNS( TcgnsContext& ctx )
{

for( int iZone = 1;  iZone <= nZones;  ++iZone )
{
	TZone& blk = Zones[iZone-1];

	// Original block size without ghosts
	const int
		nx0 = blk.ie - blk.is + 1,
		ny0 = blk.je - blk.js + 1;

	// Number of BCs in the Zone
	int nBCs = 0;  cg_nbocos(ctx.fileID,ctx.iBase,iZone, &nBCs);

	for( int iBC = 1; iBC <= nBCs; ++iBC )
	{
		CG_BCType_t iBCtype;

		CG_PointSetType_t pntSet;
		cgsize_t nPnts = -1; // number of points defining the BC region

		// Normals to the BC patch
		int iNorm[2]; // normal as index vector (computational coords)
		cgsize_t normListSize;  CG_DataType_t normDataType; // normals in phys coords

		int nDatasets = 0; // number of datasets with additional info for the BC

		char szName[33];

		cg_boco_info( ctx.fileID, ctx.iBase, iZone, iBC,
			szName, &iBCtype,
			&pntSet, &nPnts,
			iNorm, &normListSize, &normDataType,
			&nDatasets
		);
		if( pntSet != CG_PointRange && nPnts != 2 )
		{
			wxLogError(
				_("Boundary condition patch '%s'(#%d) of zone '%s'(#%d) isn't defined as point range"),
				wxString::FromAscii(szName).c_str(), iBC,
				wxString::FromAscii(blk.szName).c_str(), iZone
			);
			return false;
		}

		// Family name
		cg_goto(ctx.fileID,ctx.iBase, "Zone_t",iZone, "ZoneBC",0, "BC_t",iBC, NULL);
		if( cg_famname_read(szName) == CG_OK )
		{
			// Read BC-family "Fam_Descr_Name" generated by Pointwise 16.03
			if( cg_goto(ctx.fileID,ctx.iBase, szName,0, NULL) == CG_OK )
			{
				char szDescrName[33] = "";
				char* szDescrText = NULL;
				if( cg_descriptor_read(1, szDescrName, &szDescrText) == CG_OK )
				{
					if( strcmp(szDescrName, "Fam_Descr_Name") == 0 )
					{
						strcpy(szName, szDescrText);
						cg_free(szDescrText);
					}
				}
			}
		}


		// Read BC patch point range
		cgsize_t idxRange[4];   // Imin, Jmin, Imax, Jmax

		cg_boco_read(ctx.fileID,ctx.iBase,iZone, iBC, idxRange, NULL);

		const long
			&is = idxRange[0],  &js = idxRange[1],
			&ie = idxRange[2],  &je = idxRange[3];

		// Check that BC patch cover whole face
		bool ok = true;
		ok &= (is==1 && ie==1) || (is==nx0 && ie==nx0) || (is==1 && ie==nx0);
		ok &= (js==1 && je==1) || (js==ny0 && je==ny0) || (js==1 && je==ny0);
		if( ! ok )
		{
			wxLogError(
				_("Boundary condition patch '%s'(#%d) of zone '%s'(#%d) don't cover whole face"),
				wxString::FromAscii(szName).c_str(), iBC,
				wxString::FromAscii(blk.szName).c_str(), iZone
			);
			return false;
		}

		// Detect BC face
		TZoneFacePos posFace = faceNone;
		if( is==ie ) // Xmin or Xmax block face
			posFace = ( is==1 ) ? faceXmin : faceXmax;
		else if( js==je )  // Ymin or Ymax block face
			posFace = ( js==1 ) ? faceYmin : faceYmax;

		strcpy(blk.Faces[posFace].szBCFamName, szName);
	} // for iBC
} // for iZone

	return true;
}


//-----------------------------------------------------------------------------