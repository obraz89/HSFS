#include "stdafx.h"

#include "MF_CGNS3D.h"

#include <cgnslib.h>
#include "cgns_structs.h"

using namespace mf;
using namespace mf::cg;

void t_MFCGNS3D::get_k_range(int iZone, int& ks, int& ke) const{

	const TZone& blk = Zones[iZone-1];

	ks = blk.ks;
	ke = blk.ke;
}

bool t_MFCGNS3D::loadGrid( const wxString& gridFN ){
	bool ok = false;

	wxLogMessage( _("* Grid: Loading from '%s' & processing..."), gridFN.c_str() );

	ok = _doLoadGrid3D_cgns( gridFN );
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
 *  Loads grid data from the file
 *  
 *  @param[in] gridFN
 *  
 *  @return true if succeeded, false if failed
**/  

bool t_MFCGNS3D::_doLoadGrid3D_cgns( const wxString& gridFN )
{
	int res = CG_OK;
	char szName[33];  // names in CGNS file

	TcgnsContext ctx;

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
	if( dimCell != 3 )
	{
		wxLogError( _("The grid is not for 3D problems") );
		return false;
	}

	//
	// Number of zones (aka blocks)
	// 
	nZones = 0;  cg_nzones(ctx.fileID,ctx.iBase, &nZones);
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
		ZoneType_t type;  cg_zone_type(ctx.fileID,ctx.iBase,iZone, &type);
		if( type != Structured )
		{
			wxLogError( _("Only structured grids are supported") );
			return false;
		}

		TZone& blk = Zones[iZone -1];

		// Get zone (block) name & size
		cgsize_t isize[9]; // NVertexI, NVertexJ, NVertexK,
		              // NCellI, NCellJ, NCellK,
		              // NBoundVertexI, NBoundVertexJ, NBoundVertexK
		res = cg_zone_read(ctx.fileID,ctx.iBase,iZone,  blk.szName,isize);
		if( res != CG_OK )
		{
			wxLogError( _("Can't read zone #%d info"), iZone );
			return false;
		}

		blk.nx = isize[0];  blk.ny = isize[1];   blk.nz = isize[2];

		// Indices of real (not ghost) nodes
		blk.is = 1;  blk.ie = blk.nx;
		blk.js = 1;  blk.je = blk.ny;
		blk.ks = 1;  blk.ke = blk.nz;

		ctx.map_zoneName_id[blk.szName] = iZone;
	}


	//
	// Connectivity info
	// NB: blk.{nx, ny, nz} will be updated
	//
	if( ! _parseGhostData3DfromCGNS(ctx) )  return false;

	//
	// Boundary conditions info
	//
	if( ! _parseBCData3DfromCGNS(ctx) )  return false;


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
			ny0 = blk.je - blk.js + 1,
			nz0 = blk.ke - blk.ks + 1;

		// Zone index ranges
		cgsize_t irmin[3] = {1, 1, 1};
		cgsize_t irmax[3] = {nx0, ny0, nz0};

		double* x = new double[nx0*ny0*nz0];
		res = cg_coord_read(ctx.fileID,ctx.iBase,iZone,"CoordinateX",RealDouble,irmin,irmax, x);

		if( res != CG_OK )
			wxLogError(_T("cg_ccord_read error:%s"), wxString::FromAscii(cg_get_error()).c_str());

		double* y = new double[nx0*ny0*nz0];
		res |= cg_coord_read(ctx.fileID,ctx.iBase,iZone,"CoordinateY",RealDouble,irmin,irmax, y);

		double* z = new double[nx0*ny0*nz0];
		res |= cg_coord_read(ctx.fileID,ctx.iBase,iZone,"CoordinateZ",RealDouble,irmin,irmax, z);

		if( res != CG_OK )  return false;

		//
		// Grid step in computational (curvilinear) coordinates
		// 
		blk.grd.dx = 1./(blk.nx - 1);
		blk.grd.dy = 1./(blk.ny - 1);
		blk.grd.dz = 1./(blk.nz - 1);

		//
		// Packed grid data storage
		// 
		blk.grd.c3d = new TgridCell3D[blk.nx*blk.ny*blk.nz];
		if( ! blk.grd.c3d )  wxLogFatalError(_("Not enough memory for mesh"));

		// Pack grid coordinates
		for( int k=blk.ks; k <= blk.ke; ++k )
		for( int j=blk.js; j <= blk.je; ++j )
		for( int i=blk.is; i <= blk.ie; ++i )
		{
			Tmtr3D& node = blk.cell(i,j,k).ijk;

			int n = (i - blk.is) + (j - blk.js)*nx0 + (k - blk.ks)*nx0*ny0;
			node.x = x[n];
			node.y = y[n];
			node.z = z[n];
			//node.Rw = HUGE_VAL;
		}

			delete[] x, y, z;
	} // loop through zones


	//
	// Transfer grid data to ghost nodes
	// 
	for( int b = 0; b < nZones; ++b )
	{
		TZone& zne = Zones[b];

		for( int k = 1; k <= zne.nz; ++k )
		for( int j = 1; j <= zne.ny; ++j )
		for( int i = 1; i <= zne.nx; ++i )
		{
			if( i >= zne.is && i <= zne.ie &&
				j >= zne.js && j <= zne.je &&
				k >= zne.ks && k <= zne.ke    )
				continue; // myself

			const int locInd = zne.absIdx(i,j) - 1;
			const int globInd = zne.globIndices[ locInd ];
			if( globInd == -1 )  // unmapped
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

			zne.grd.c3d[ locInd ] = dnrZne.grd.c3d[ dnrLocInd ];
		}
	}

	cg_close( ctx.fileID );
	return true;
}
//-----------------------------------------------------------------------------

// ~_doLoadGrid_cgns from HSFlow(20140508)
// version of _doLoadGrid_cgns from _HSFlow(20140508)


/**
 *  Read connectivity data from the opened CGNS file
 *  and initialize ghost indices 
 *  
 *  @param[in] ctx - context of the opened CGNS file
 *  
 *  @return true if succeeded, false if failed
**/  
bool t_MFCGNS3D::_parseGhostData3DfromCGNS( TcgnsContext& ctx )
{
	// Number of ghost nodes (half of stencil size)
	int ghostI = 2, ghostJ = 2, ghostK = 2;
	// TODO:
	//G_Plugins.get_discret_caps().getParI(11/*nsx*/, ghostI);  ghostI /= 2;
	//G_Plugins.get_discret_caps().getParI(24/*nsy*/, ghostJ);  ghostJ /= 2;

//
// (1) Read connectivity data from CGNS file
// 
for( int iZone = 1;  iZone <= nZones;  ++iZone )
{
	TZone& blk = Zones[iZone-1];
	TcgnsZone& cgZne = ctx.cgZones[iZone-1];

	int n1to1 = 0;  cg_n1to1(ctx.fileID,ctx.iBase,iZone, &n1to1);
	if( n1to1 < 1 && nZones > 1 )
	{
		wxLogError(_("Missing inter-block 1-to-1 connectivity info"));
		return false;
	}

	// Original block size without ghosts
	const int
		nx0 = blk.ie - blk.is + 1,
		ny0 = blk.je - blk.js + 1,
		nz0 = blk.ke - blk.ks + 1;


	// Loop through connectivity patches
	for( int iConn = 1; iConn <= n1to1; ++iConn )
	{
		cgsize_t idxRange[6];   // Imin, Jmin, Kmin, Imax, Jmax, Kmin

		cgsize_t idxRangeDonor[6];

		char szName[33], szDonor[33];
		int iTsh[3]; // short-hand notation of transform matrix
		bool res = cg_1to1_read( ctx.fileID,ctx.iBase,iZone,iConn,
			szName, szDonor, idxRange, idxRangeDonor, iTsh
		);
		if( res != CG_OK )
		{
			wxLogError( _("Can't read connection patch '%s'(#%d) of zone '%s'(#%d)"),
				wxString::FromAscii(szName).c_str(), iConn,
				wxString::FromAscii(blk.szName).c_str(), iZone
			);
			return false;
		}

		const long
			&is = idxRange[0],  &js = idxRange[1],  &ks = idxRange[2],
			&ie = idxRange[3],  &je = idxRange[4],  &ke = idxRange[5];

		const long
			&ids = idxRangeDonor[0],  &jds = idxRangeDonor[1],  &kds = idxRangeDonor[2],
			&ide = idxRangeDonor[3],  &jde = idxRangeDonor[4],  &kde = idxRangeDonor[5];

		// Check that whole faces connected
		bool ok = true;
		ok &= (is==1 && ie==1) || (is==nx0 && ie==nx0) || (is==1 && ie==nx0);
		ok &= (js==1 && je==1) || (js==ny0 && je==ny0) || (js==1 && je==ny0);
		ok &= (ks==1 && ke==1) || (ks==nz0 && ke==nz0) || (ks==1 && ke==nz0);
		if( ! ok )
		{
			wxLogError(
				_("Connection patch '%s'(#%d) of zone '%s'(#%d) don't cover whole block face"),
				wxString::FromAscii(szName).c_str(), iConn,
				wxString::FromAscii(blk.szName).c_str(), iZone
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
			blk.nx += ghostI;
			if( is==1 )
			{
				posFace = faceXmin;
				blk.is += ghostI;
				blk.ie += ghostI;
			}
			else
				posFace = faceXmax;

		}
		else if( js==je )  // Ymin or Ymax block face
		{
			blk.ny += ghostJ;
			if( js==1 )
			{
				posFace = faceYmin;
				blk.js += ghostJ;
				blk.je += ghostJ;
			}
			else
				posFace = faceYmax;
		}
		else if( ks==ke )  // Zmin or Zmax block face
		{
			blk.nz += ghostK;
			if( ks==1 )
			{
				posFace = faceZmin;
				blk.ks += ghostK;
				blk.ke += ghostK;
			}
			else
				posFace = faceZmax;
		}

		blk.Faces[posFace].bcType = bcAbutted;

		TcgnsZone::TFace& face = cgZne.Faces[posFace];
		face.nDnrZne = ctx.map_zoneName_id[szDonor] - 1; // 0-based

		face.dnr_is0 = ids;
		face.dnr_js0 = jds;
		face.dnr_ks0 = kds;

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
		face.matTrans[2] = F::sgn(iTsh[2]) * F::del(iTsh[2],1);

		face.matTrans[3] = F::sgn(iTsh[0]) * F::del(iTsh[0],2);
		face.matTrans[4] = F::sgn(iTsh[1]) * F::del(iTsh[1],2);
		face.matTrans[5] = F::sgn(iTsh[2]) * F::del(iTsh[2],2);

		face.matTrans[6] = F::sgn(iTsh[0]) * F::del(iTsh[0],3);
		face.matTrans[7] = F::sgn(iTsh[1]) * F::del(iTsh[1],3);
		face.matTrans[8] = F::sgn(iTsh[2]) * F::del(iTsh[2],3);
	} // loop through connection patches
} // loop through zones



//
// (2) Generate global indices of real nodes
// 
int gIdx = 0;
for( int b = 0; b < nZones; ++b )
{
	TZone& zne = Zones[b];

	zne.nGlobStart = gIdx;

	zne.globIndices = new int[ zne.nx * zne.ny * zne.nz ];
	memset( zne.globIndices, 0xFF, zne.nx* zne.ny* zne.nz * sizeof(int) );  // -1

	// Fill-in global indices of real nodes
	for( int k=zne.ks; k<=zne.ke; ++k )
	for( int j=zne.js; j<=zne.je; ++j )
	for( int i=zne.is; i<=zne.ie; ++i )
	{
		const int locInd = zne.absIdx(i,j,k) - 1;  // 0-based
		zne.globIndices[locInd] = zne.globRealInd(i,j,k);
	}

	// Real nodes count (without ghosts)
	const int
		nx0 = zne.ie - zne.is + 1,
		ny0 = zne.je - zne.js + 1,
		nz0 = zne.ke - zne.ks + 1;

	gIdx += nx0*ny0*nz0;
}


//
// (3) Parse connectivity data and fill in ghost nodes indices
//
for( int b = 0; b < nZones; ++b )
{
	TZone& zne = Zones[b];

	// Loop through faces
	for( int f = 0; f < 6; ++f )
	{
		if( zne.Faces[f].bcType != bcAbutted ) continue;

		// Start & End indexes of ghost layer
		int igs, ige,   jgs, jge,   kgs, kge;

		// Start indexes of the block face 'bF'
		int ifs, jfs, kfs;

		switch( f )
		{
		case faceXmin:
			ifs = zne.is;      jfs = zne.js;      kfs = zne.ks;
			igs = 1;           jgs = zne.js;      kgs = zne.ks;
			ige = zne.is - 1;  jge = zne.je;      kge = zne.ke;
			break;
		case faceXmax:
			ifs = zne.ie;      jfs = zne.js;      kfs = zne.ks;
			igs = zne.ie + 1;  jgs = zne.js;      kgs = zne.ks;
			ige = zne.nx;      jge = zne.je;      kge = zne.ke;
			break;
		case faceYmin:
			ifs = zne.is;      jfs = zne.js;      kfs = zne.ks;
			igs = zne.is;      jgs = 1;           kgs = zne.ks;
			ige = zne.ie;      jge = zne.js - 1;  kge = zne.ke;
			break;
		case faceYmax:
			ifs = zne.is;      jfs = zne.je;      kfs = zne.ks;
			igs = zne.is;      jgs = zne.je + 1;  kgs = zne.ks;
			ige = zne.ie;      jge = zne.ny;      kge = zne.ke;
			break;
		case faceZmin:
			ifs = zne.is;      jfs = zne.js;      kfs = zne.ks;
			igs = zne.is;      jgs = zne.js;      kgs = 1;
			ige = zne.ie;      jge = zne.je;      kge = zne.ks - 1;
			break;
		case faceZmax:
			ifs = zne.is;      jfs = zne.js;      kfs = zne.ke;
			igs = zne.is;      jgs = zne.js;      kgs = zne.ke + 1;
			ige = zne.ie;      jge = zne.je;      kge = zne.nz;
			break;
		}

		TcgnsZone::TFace& face = ctx.cgZones[b].Faces[f];
		TZone& zneDonor = Zones[ face.nDnrZne ];

		// Donor block face starting index with ghosts
		int ids = face.dnr_is0 + zneDonor.is - 1;
		int jds = face.dnr_js0 + zneDonor.js - 1;
		int kds = face.dnr_ks0 + zneDonor.ks - 1;

		for( int k=kgs; k<=kge; ++k )
		for( int j=jgs; j<=jge; ++j )
		for( int i=igs; i<=ige; ++i )
		{
			// Donor block indexes
			int id = (i-ifs)*face.matTrans[0] + (j-jfs)*face.matTrans[1] + (k-kfs)*face.matTrans[2] + ids;
			int jd = (i-ifs)*face.matTrans[3] + (j-jfs)*face.matTrans[4] + (k-kfs)*face.matTrans[5] + jds;
			int kd = (i-ifs)*face.matTrans[6] + (j-jfs)*face.matTrans[7] + (k-kfs)*face.matTrans[8] + kds;

			if( id > zneDonor.ie || id < zneDonor.is ||
				jd > zneDonor.je || jd < zneDonor.js ||
				kd > zneDonor.ke || kd < zneDonor.ks
				)
			{
				wxLogError(
					_("Ghost node (%d,%d,%d) in zone '%s' tried to be mapped to \
					  nonexistent node (%d,%d,%d) in zone '%s'. Check indices orientation!"),
					i,  j,  k,  wxString::FromAscii(zne.szName).c_str(),
					id, jd, kd, wxString::FromAscii(zneDonor.szName).c_str()
				);
				return false;
			}

			const int locInd = zne.absIdx(i,j,k) - 1;
			zne.globIndices[ locInd ] = zneDonor.globRealInd(id,jd,kd);
		}

	} 	// loop through faces

} // loop through zones


//
// (4) Fill-in remaining ghost node indices around corners
//
for( int b = 0; b < nZones; ++b )
{
	TZone& zne = Zones[b];

	// Loop through faces
	for( int f = 0; f < 6; ++f )
	{
		if( zne.Faces[f].bcType != bcAbutted ) continue;

		// Start & End indexes of ghost layer
		int igs, ige,   jgs, jge,   kgs, kge;

		// Start indexes of the block face 'bF'
		int ifs, jfs, kfs;

		switch( f )
		{
		case faceXmin:
			ifs = zne.is;      jfs = zne.js;      kfs = zne.ks;
			igs = 1;           jgs = 1;           kgs = 1;
			ige = zne.is - 1;  jge = zne.ny;      kge = zne.nz;
			break;
		case faceXmax:
			ifs = zne.ie;      jfs = zne.js;      kfs = zne.ks;
			igs = zne.ie + 1;  jgs = 1;           kgs = 1;
			ige = zne.nx;      jge = zne.ny;      kge = zne.nz;
			break;
		case faceYmin:
			ifs = zne.is;      jfs = zne.js;      kfs = zne.ks;
			igs = 1;           jgs = 1;           kgs = 1;
			ige = zne.nx;      jge = zne.js - 1;  kge = zne.nz;
			break;
		case faceYmax:
			ifs = zne.is;      jfs = zne.je;      kfs = zne.ks;
			igs = 1;           jgs = zne.je + 1;  kgs = 1;
			ige = zne.nx;      jge = zne.ny;      kge = zne.nz;
			break;
		case faceZmin:
			ifs = zne.is;      jfs = zne.je;      kfs = zne.ks;
			igs = 1;           jgs = 1;           kgs = 1;
			ige = zne.nx;      jge = zne.ny;      kge = zne.ks - 1;
			break;
		case faceZmax:
			ifs = zne.ie;      jfs = zne.js;      kfs = zne.ks;
			igs = 1;           jgs = 1;           kgs = zne.ke + 1;
			ige = zne.nx;      jge = zne.ny;      kge = zne.nz;
			break;
		}

		TcgnsZone::TFace& face = ctx.cgZones[b].Faces[f];
		TZone& zneDonor = Zones[ face.nDnrZne ];

		// Donor block face starting index with ghosts
		int ids = face.dnr_is0 + zneDonor.is - 1;
		int jds = face.dnr_js0 + zneDonor.js - 1;
		int kds = face.dnr_ks0 + zneDonor.ks - 1;

		for( int k=kgs; k<=kge; ++k )
		for( int j=jgs; j<=jge; ++j )
		for( int i=igs; i<=ige; ++i )
		{
			const int locInd = zne.absIdx(i,j,k) - 1;
			if( zne.globIndices[ locInd ] != -1 ) continue;

			// Donor block indexes
			int id = (i-ifs)*face.matTrans[0] + (j-jfs)*face.matTrans[1] + (k-kfs)*face.matTrans[2] + ids;
			int jd = (i-ifs)*face.matTrans[3] + (j-jfs)*face.matTrans[4] + (k-kfs)*face.matTrans[5] + jds;
			int kd = (i-ifs)*face.matTrans[6] + (j-jfs)*face.matTrans[7] + (k-kfs)*face.matTrans[8] + kds;
			if( id > zneDonor.nx || id < 1 ||
				jd > zneDonor.ny || jd < 1 ||
				kd > zneDonor.nz || kd < 1    )
			{
				// ghost data unavailable
				wxLogError(_("Provided multi-block layout is not supported yet. Sorry"));
				return false;
				//continue;
			}

			int dnrLocInd = zneDonor.absIdx(id,jd,kd) - 1;  // 0-based

			zne.globIndices[locInd] = zneDonor.globIndices[ dnrLocInd ];
		}

	} 	// loop through faces

} // loop through zones


	return true;
}
//-----------------------------------------------------------------------------


bool t_MFCGNS3D::_parseBCData3DfromCGNS( TcgnsContext& ctx )
{

for( int iZone = 1;  iZone <= nZones;  ++iZone )
{
	TZone& blk = Zones[iZone-1];

	// Original block size without ghosts
	const int
		nx0 = blk.ie - blk.is + 1,
		ny0 = blk.je - blk.js + 1,
		nz0 = blk.ke - blk.ks + 1;

	// Number of BCs in the Zone
	int nBCs = 0;  cg_nbocos(ctx.fileID,ctx.iBase,iZone, &nBCs);

	for( int iBC = 1; iBC <= nBCs; ++iBC )
	{
		BCType_t iBCtype;

		PointSetType_t pntSet;
		cgsize_t nPnts = -1; // number of points defining the BC region

		// Normals to the BC patch
		int iNorm[3]; // normal as index vector (computational coords)
		cgsize_t normListSize;  DataType_t normDataType; // normals in phys coords

		int nDatasets = 0; // number of datasets with additional info for the BC

		char szName[33];

		cg_boco_info( ctx.fileID, ctx.iBase, iZone, iBC,
			szName, &iBCtype,
			&pntSet, &nPnts,
			iNorm, &normListSize, &normDataType,
			&nDatasets
			);
		if( pntSet != PointRange && nPnts != 2 )
		{
			wxLogError(
				_("Boundary condition patch '%s'(#%d) of zone '%s'(#%d) isn't defined as point range"),
				wxString::FromAscii(szName).c_str(), iBC,
				wxString::FromAscii(blk.szName).c_str(), iZone
			);
			return false;
		}
		/*
		if( iBCtype == FamilySpecified )
		{
			// Support for Pointwise version >= 16.04
			// It writes meaningless BC names, but correct BC FamilyName
			cg_goto(ctx.fileID,ctx.iBase, "Zone_t",iZone, "ZoneBC",0, "BC_t",iBC, NULL);
			cg_famname_read(szName);
		}*/

		// Family name
		cg_goto(ctx.fileID,ctx.iBase, "Zone_t",iZone, "ZoneBC",0, "BC_t",iBC, NULL);
		if( cg_famname_read(szName) == CG_OK )
		{
			// Read "Fam_Descr_Name" generated by Pointwise 16.03
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
		cgsize_t idxRange[6];   // Imin, Jmin, Kmin, Imax, Jmax, Kmax

		cg_boco_read(ctx.fileID,ctx.iBase,iZone, iBC, idxRange, NULL);

		// Check that BC patch cover whole face

		const long
			&is = idxRange[0],  &js = idxRange[1],  &ks = idxRange[2],
			&ie = idxRange[3],  &je = idxRange[4],  &ke = idxRange[5];

		bool ok = true;
		ok &= (is==1 && ie==1) || (is==nx0 && ie==nx0) || (is==1 && ie==nx0);
		ok &= (js==1 && je==1) || (js==ny0 && je==ny0) || (js==1 && je==ny0);
		ok &= (ks==1 && ke==1) || (ks==nz0 && ke==nz0) || (ks==1 && ke==nz0);
		if( ! ok )
		{
			wxLogError(
				_("Boundary condition patch '%s'(#%d) of zone '%s'(#%d) don't cover whole block face"),
				wxString::FromAscii(szName).c_str(), iBC,
				wxString::FromAscii(blk.szName).c_str(), iZone
			);
			return false;
		}

		// Detect BC face
		TZoneFace* bF = NULL;
		if( is==ie ) // Xmin or Xmax block face
			bF = ( is==1 ) ? &blk.Faces[faceXmin] : &blk.Faces[faceXmax];
		else if( js==je )  // Ymin or Ymax block face
			bF = ( js==1 ) ? &blk.Faces[faceYmin] : &blk.Faces[faceYmax];
		else if( ks==ke )  // Zmin or Zmax block face
			bF = ( ks==1 ) ? &blk.Faces[faceZmin] : &blk.Faces[faceZmax];

		strcpy(bF->szBCFamName, szName);
	} // for iBC
} // for iZone

	return true;
}
//-----------------------------------------------------------------------------
