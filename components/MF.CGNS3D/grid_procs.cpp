#include "stdafx.h"

#include "MF_CGNS3D.h"

#include <cgnslib.h>
#include "cgns_structs.h"

using namespace cgns_my;
using namespace mf;

/**
 *  Loads grid data from the file
 *  
 *  @param[in] gridFN
 *  
 *  @return true if succeeded, false if failed
**/  

bool t_MFCGNS3D::_load_grid_data_3D(const wxString& gridFN){

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
	G_Domain.nZones = 0;  cg_nzones(ctx.fileID,ctx.iBase, &G_Domain.nZones);
	if( G_Domain.nZones < 1 )
	{
		wxLogError( _("Domain blocks are not found") );
		return false;
	}

	//
	// Allocate memory for whole computational domain
	// 
	G_Domain.Zones = new TZone[ G_Domain.nZones ];
	ctx.cgZones = new TcgnsZone[ G_Domain.nZones ];
	if( ! G_Domain.Zones )  wxLogFatalError(_("Not enough memory for domain"));


	//
	// Get zones (blocks) sizes and names
	// 
	for( int iZone = 1;  iZone <= G_Domain.nZones;  ++iZone )
	{
		ZoneType_t type;  cg_zone_type(ctx.fileID,ctx.iBase,iZone, &type);
		if( type != Structured )
		{
			wxLogError( _("Only structured grids are supported") );
			return false;
		}

		TZone& blk = G_Domain.Zones[iZone -1];

		// Get zone (block) name & size

		int isize[9]; // NVertexI, NVertexJ, NVertexK,
		// NCellI, NCellJ, NCellK,
		// NBoundVertexI, NBoundVertexJ, NBoundVertexK
		res = cg_zone_read(ctx.fileID,ctx.iBase,iZone,  blk.szName,isize);
		if( res != CG_OK )
		{
			wxLogError( _("Can't read zone #%d info"), iZone );
			return false;
		}

		blk.nx = isize[0];  blk.ny = isize[1];   blk.nz = isize[2];

		// TODO: proper place
		Nx = blk.nx;
		Ny = blk.ny;
		Nz = blk.nz;

		// Indices of real (not ghost) nodes
		blk.is = 1;  blk.ie = blk.nx;
		blk.js = 1;  blk.je = blk.ny;
		blk.ks = 1;  blk.ke = blk.nz;

		ctx.map_zoneName_id[blk.szName] = iZone;

		// Initial face names, will be redefined below
		sprintf(blk.Faces[faceXmin].szName, "blk%d-Xmin", iZone);
		sprintf(blk.Faces[faceXmax].szName, "blk%d-Xmax", iZone);
		sprintf(blk.Faces[faceYmin].szName, "blk%d-Ymin", iZone);
		sprintf(blk.Faces[faceYmax].szName, "blk%d-Ymax", iZone);
		sprintf(blk.Faces[faceZmin].szName, "blk%d-Zmin", iZone);
		sprintf(blk.Faces[faceZmax].szName, "blk%d-Zmax", iZone);
	}


	//
	// Connectivity info
	// NB: blk.{nx, ny, nz} will be updated
	//
	if( ! parseGhostData3DfromCGNS(ctx) )  return false;

	//
	// Boundary conditions info
	//
	if( ! parseBCData3DfromCGNS(ctx) )  return false;


	//
	// Read grid coordinates for all zones
	// (to simplify mesh data transfer to ghosts)
	//
	for( int b = 0;  b < G_Domain.nZones;  ++b )
	{
		int iZone = b + 1;
		TZone& blk = G_Domain.Zones[b];

		// Original block size without ghosts
		const int
			nx0 = blk.ie - blk.is + 1,
			ny0 = blk.je - blk.js + 1,
			nz0 = blk.ke - blk.ks + 1;

		// Zone index ranges
		int irmin[3] = {1, 1, 1};
		int irmax[3] = {nx0, ny0, nz0};

		double* x = new double[nx0*ny0*nz0];
		res = cg_coord_read(ctx.fileID,ctx.iBase,iZone,"CoordinateX",RealDouble,irmin,irmax, x);

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
					node.Rw = HUGE_VAL;
				}

				delete[] x, y, z;
	} // loop through zones


	//
	// Transfer grid data to ghost nodes
	// 
	for( int b = 0; b < G_Domain.nZones; ++b )
	{
		TZone& zne = G_Domain.Zones[b];

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
					for( nDnrZne = G_Domain.nZones - 1; nDnrZne >= 0; --nDnrZne )
					{
						if( globInd >= G_Domain.Zones[ nDnrZne ].nGlobStart )
							break;
					}

					TZone& dnrZne = G_Domain.Zones[ nDnrZne ];
					int dnrLocInd = dnrZne.absIdxFromGlobInd( globInd ) - 1; // 0-based

					zne.grd.c3d[ locInd ] = dnrZne.grd.c3d[ dnrLocInd ];
				}
	}

	cg_close( ctx.fileID );
	return true;

}
