#include "stdafx.h"

#include "MF_CGNS2D.h"

#include <cgnslib.h>
#include "cgns_structs.h"

using namespace mf;
using namespace mf::cg;

#if CG_BUILD_SCOPE
  #define CG_MY_RealDouble CG_RealDouble
  #define CG_MY_ZoneType_t CG_ZoneType_t
  #define CG_MY_Structured CG_Structured
  #define CG_MY_DataType_t CG_DataType_t
  #define CG_MY_BCType_t CG_BCType_t
  #define CG_MY_PointSetType_t CG_PointSetType_t
  #define CG_MY_PointRange CG_PointRange
#else
  #define CG_MY_RealDouble RealDouble
  #define CG_MY_ZoneType_t ZoneType_t
  #define CG_MY_Structured Structured
  #define CG_MY_DataType_t DataType_t
  #define CG_MY_BCType_t BCType_t
  #define CG_MY_PointSetType_t PointSetType_t
  #define CG_MY_PointRange PointRange
#endif

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

	TcgnsContext& ctx = cgCtx;

	res = cg_open(gridFN.ToAscii(), CG_MODE_READ, &ctx.fileID);
	if (res != CG_OK)
	{
		wxLogError(
			_("Can't open grid file '%s' for reading (%s)"),
			gridFN.c_str(), wxString::FromAscii(cg_get_error()).c_str()
		);
		return false;
	}
	ctx.iBase = 1; // assume only one base

				   //
				   // Space dimensions
				   // 
	int dimCell = 0, dimPhys = 0;
	cg_base_read(ctx.fileID, ctx.iBase, szName, &dimCell, &dimPhys);
	if (dimCell != 2)
	{
		wxLogError(_("The grid is not for 2D problems"));
		return false;
	}

	//
	// Number of zones
	//
	nZones = 0;  cg_nzones(ctx.fileID, ctx.iBase, &nZones);
	if (nZones < 1)
	{
		wxLogError(_("Domain zones are not found"));
		return false;
	}

	//
	// Allocate memory for whole computational domain
	//
	try {
		Zones = new TZone[nZones];
	}
	catch (std::bad_alloc& ba)
	{
		wxLogError(_("Can't allocate memory for domain zones container. Asked for %.3f MB."),
			(float)sizeof(TZone) * nZones / (1024.*1024.)
		);
		return false;
	}


	//
	// Get zones sizes and names
	//
	for (int iZone = 1; iZone <= nZones; ++iZone)
	{
		CG_ZoneType_t type;  cg_zone_type(ctx.fileID, ctx.iBase, iZone, &type);
		if (type != CG_Structured)
		{
			wxLogError(_("Only structured grids are supported"));
			return false;
		}

		TZone& zne = Zones[iZone - 1];

		// Get zone name & size
		cgsize_t isize[6]; // NVertexI, NVertexJ, NCellI, NCellJ, NBoundVertexI, NBoundVertexJ
		res = cg_zone_read(ctx.fileID, ctx.iBase, iZone, zne.szName, isize);
		if (res != CG_OK)
		{
			wxLogError(_("Can't read zone #%d info"), iZone);
			return false;
		}

		zne.nx = isize[0];  zne.ny = isize[1];   zne.nz = 1;

		// Indices of real (not ghost) nodes
		zne.is = 1;  zne.ie = zne.nx;
		zne.js = 1;  zne.je = zne.ny;
		zne.ks = 1;  zne.ke = 1;

		ctx.map_zoneName_id[zne.szName] = iZone;
	}

	//

	//
	// Inter-zones connectivity info
	// Updates {zne,cgZne}.{is,ie,js,je} & zne.{nx,ny}
	if (!_parseGhostData2DfromCGNS(ctx))  return false;


	// Boundary conditions info
	if (!_parseBCData2DfromCGNS(ctx))  return false;

	// bs and be kept from Nova's code for compatibility
	// One core distribution
	this->bs = 0;
	this->be = nZones - 1;

	//
	// Read grid coordinates
	//
	for (int b = bs; b <= be; ++b)
	{
		int iZone = b + 1;
		TZone& zne = Zones[b];
		TcgnsZone& cgZne = ctx.cgZones[b];

		// Original zone size (not accounting added ghosts & skipped layers)
		const int
			nx0 = cgZne.ie1 - cgZne.is1 + 1,
			ny0 = cgZne.je1 - cgZne.js1 + 1;

		// Zone index ranges
		cgsize_t irmin[2] = { 1, 1 };
		cgsize_t irmax[2] = { nx0, ny0 };

		double* x = new double[nx0*ny0];
		res = cg_coord_read(ctx.fileID, ctx.iBase, iZone, "CoordinateX", CG_RealDouble, irmin, irmax, x);

		double* y = new double[nx0*ny0];
		res |= cg_coord_read(ctx.fileID, ctx.iBase, iZone, "CoordinateY", CG_RealDouble, irmin, irmax, y);

		if (res != CG_OK)  return false;

		{
			zne.grd.dx = 1. / (zne.nx - 1);
			zne.grd.dy = 1. / (zne.ny - 1);
		}

		//
		// Packed grid data storage
		//
		try {
			zne.grd.c2d = new TgridCell2D[zne.nx * zne.ny];
		}
		catch (std::bad_alloc& ba)
		{
			wxLogError(_("Can't allocate memory for grid data in zone '%s'#%d: asked for %.3f MB."),
				wxString::FromAscii(zne.szName).c_str(), iZone,
				(float)sizeof(TgridCell2D) * zne.nx*zne.ny / (1024.*1024.)
			);
			return false;
		}

		// Pack grid coordinates
		for (int j = 1; j <= zne.ny; ++j)
			for (int i = 1; i <= zne.nx; ++i)
			{
				Tmtr2D& node = zne.cell(i, j).ij;

				if (j >= zne.js && j <= zne.je
					&&
					i >= zne.is && i <= zne.ie
					)
				{
					int n = (i - cgZne.is1) + (j - cgZne.js1)*nx0;
					node.x = x[n];
					node.y = y[n];
				}
				else
				{
					static const double NaN = fmod(1., 0.);
					node.x = NaN;
					node.y = NaN;
				}
			}

		delete[] x;   delete[] y;
	} // loop zones

	cg_close(ctx.fileID);

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

	ctx.cgZones = new TcgnsZone[nZones];

	//
	// (1) Read connectivity data from CGNS file
	//
	for (int iZone = 1; iZone <= nZones; ++iZone)
	{
		TZone& zne = Zones[iZone - 1];
		TcgnsZone& cgZne = ctx.cgZones[iZone - 1];

		cg_n1to1(ctx.fileID, ctx.iBase, iZone, &cgZne.nPatches);
		if (cgZne.nPatches < 1 && nZones > 1)
		{
			wxLogError(_("Missing inter-zone 1-to-1 connectivity info"));
			return false;
		}

		cgZne.Patches = new TcgnsZone::TFacePatch[cgZne.nPatches];

		// Loop through connectivity patches
		for (int iPatch = 1; iPatch <= cgZne.nPatches; ++iPatch)
		{
			TcgnsZone::TFacePatch& cgFacePatch = cgZne.Patches[iPatch - 1];

			char szName[33], szDonor[33];
			cgsize_t idxRange[4], idxRangeDonor[4];
			int iTsh[2]; // short-hand notation of transform matrix

			if (cg_1to1_read(ctx.fileID, ctx.iBase, iZone, iPatch,
				szName, szDonor, idxRange, idxRangeDonor, iTsh) != CG_OK)
			{
				wxLogError(_("Can't read connection patch '%s'(#%d) of zone '%s'(#%d)"),
					wxString::FromAscii(szName).c_str(), iPatch,
					wxString::FromAscii(zne.szName).c_str(), iZone
				);
				return false;
			}

			const cgsize_t
				&is = idxRange[0], &js = idxRange[1],
				&ie = idxRange[2], &je = idxRange[3];
			const cgsize_t
				&ids = idxRangeDonor[0], &jds = idxRangeDonor[1],
				&ide = idxRangeDonor[2], &jde = idxRangeDonor[3];

			//
			// Detect face, the current connection patch belongs to
			// 
			if (is == ie)  cgFacePatch.posFace = (is == 1) ? faceXmin : faceXmax;
			else if (js == je)  cgFacePatch.posFace = (js == 1) ? faceYmin : faceYmax;

			cgFacePatch.is0 = is;   cgFacePatch.ie0 = ie;
			cgFacePatch.js0 = js;   cgFacePatch.je0 = je;

			//
			// Increase the zone size by ghosts count
			//
			if (zne.Faces[cgFacePatch.posFace].bcType != bcAbutted) // not processed yet
			{
				zne.Faces[cgFacePatch.posFace].bcType = bcAbutted;

				if (zne.isFrozen) continue;

				switch (cgFacePatch.posFace)
				{
				case faceXmin:  zne.nx += ghostI;  zne.is += ghostI;  zne.ie += ghostI;  break;
				case faceXmax:  zne.nx += ghostI;                                        break;

				case faceYmin:  zne.ny += ghostJ;  zne.js += ghostJ;  zne.je += ghostJ;  break;
				case faceYmax:  zne.ny += ghostJ;                                        break;
				}
			}

			//
			// Detect donor zone
			//
			if (ctx.map_zoneName_id.find(szDonor) == ctx.map_zoneName_id.end())
			{
				wxLogError(
					_("Connection patch '%s'#%d of zone '%s'#%d abuts to nonexistent zone '%s'"),
					wxString::FromAscii(szName).c_str(), iPatch,
					wxString::FromAscii(zne.szName).c_str(), iZone,
					wxString::FromAscii(szDonor).c_str()
				);
				return false;
			}
			cgFacePatch.nDnrZne = ctx.map_zoneName_id[szDonor] - 1;  // 0-based

																	 //
																	 // Detect face, the donor connection patch belongs to.
																	 // Check if connectivity info is correct
																	 // 
			const TZone& zneDnr = Zones[cgFacePatch.nDnrZne];
			if (ids == ide)
			{
				if (ids == 1)
					cgFacePatch.dnrPosFace = faceXmin;
				else if (ids == (zneDnr.ie - zneDnr.is + 1))
					cgFacePatch.dnrPosFace = faceXmax;
			}
			else if (jds == jde)
			{
				if (jds == 1)
					cgFacePatch.dnrPosFace = faceYmin;
				if (jds == (zneDnr.je - zneDnr.js + 1))
					cgFacePatch.dnrPosFace = faceYmax;
			}

			if (cgFacePatch.dnrPosFace == faceNone)
			{
				wxLogError(
					_("Connection patch '%s'#%d of zone '%s'#%d doesn't abut to any face"),
					wxString::FromAscii(szName).c_str(), iPatch,
					wxString::FromAscii(zne.szName).c_str(), iZone
				);
				return false;
			}

			cgFacePatch.dnr_is0 = ids;
			cgFacePatch.dnr_js0 = jds;

			// Transform matrix between node indexes
			// of the current and donor zones
			// 
			// CGNS SIDS documentation,
			// section 8.2 1-to-1 Interface Connectivity Structure Definition
			// 
			// iTsh = [a, b, c] =>
			//	   | sgn(a)del(a,1) sgn(b)del(b,1) sgn(c)del(c,1) |
			// T = | sgn(a)del(a,2) sgn(b)del(b,2) sgn(c)del(c,2) |,
			//	   | sgn(a)del(a,3) sgn(b)del(b,3) sgn(c)del(c,3) |
			struct F {
				static int sgn(int x) { return (x<0) ? -1 : +1; }
				static int del(int x, int y) { return (abs(x) == abs(y)) ? +1 : 0; }
			};

			int* mT = cgFacePatch.matTrans;
			mT[0] = F::sgn(iTsh[0]) * F::del(iTsh[0], 1);
			mT[1] = F::sgn(iTsh[1]) * F::del(iTsh[1], 1);

			mT[2] = F::sgn(iTsh[0]) * F::del(iTsh[0], 2);
			mT[3] = F::sgn(iTsh[1]) * F::del(iTsh[1], 2);
		} // loop connection patches
	} // loop zones


	//
	// (2) Assign owners for abutted faces
	// 

	//
	// (2.1) Mark faces that should be skipped and processed by abutted zone
	//
	// NOTE: skipped faces disabled
	for (int b = 0; b < nZones; ++b)
	{
		TZone& zne = Zones[b];
		TcgnsZone& cgZne = ctx.cgZones[b];

		for (int p = 0; p < cgZne.nPatches; ++p) // face patches
		{
			TcgnsZone::TFacePatch& cgFacePatch = cgZne.Patches[p];
			TZone& zneDnr = Zones[cgFacePatch.nDnrZne];

			if (!zne.isFrozen && zneDnr.isFrozen)
			{
				zne.Faces[cgFacePatch.posFace].isSkipped = false; // true;
				continue;
			}

			// Lengths of the abutted faces with ghosts
			int lenMy = 0;
			switch (cgFacePatch.posFace)
			{
			case faceXmin:   case faceXmax:
				lenMy = zne.ny;
				break;
			case faceYmin:   case faceYmax:
				lenMy = zne.nx;
				break;
			}

			int lenDnr = 0;
			switch (cgFacePatch.dnrPosFace)
			{
			case faceXmin:   case faceXmax:
				lenDnr = zneDnr.ny;
				break;
			case faceYmin:   case faceYmax:
				lenDnr = zneDnr.nx;
				break;
			}

			// Skip face for processing by the abutted zone,
			// if donor face is shorter, i.e. has less ghost nodes, has boundary nearby
			if (lenDnr < lenMy)
				zne.Faces[cgFacePatch.posFace].isSkipped = false; // true;
		}
	}
	


	// FIXME: Check triple-points
	//
	// (2.2) Resolve conflicts if faces are owned by both abutted zones
	// 
	// NOTE: skipped faces disabled
	for (int b = 0; b < nZones; ++b)
	{
		TcgnsZone& cgZne = ctx.cgZones[b];
		for (int p = 0; p < cgZne.nPatches; ++p) // abutted face patches
		{
			TcgnsZone::TFacePatch& cgPatch = cgZne.Patches[p];

			TZoneFace& face = Zones[b].Faces[cgPatch.posFace];
			TZoneFace& faceDnr = Zones[cgPatch.nDnrZne].Faces[cgPatch.dnrPosFace];

			if (!face.isSkipped && !faceDnr.isSkipped)
				faceDnr.isSkipped = false; // true;
		}
	} // loop zones

	


	//
	// (3) Update domain zones dimensions and indices using abutted faces info
	//
	// NOTE: skipped faces disabled
	for (int b = 0; b < nZones; ++b)
	{
		TZone& zne = Zones[b];
		TcgnsZone& cgZne = ctx.cgZones[b];

		cgZne.is1 = zne.is;   cgZne.ie1 = zne.ie;
		cgZne.js1 = zne.js;   cgZne.je1 = zne.je;

		for (int f = 0; f < 4; ++f) // faces
		{
			TZoneFace& face = zne.Faces[f];
			if (face.bcType != bcAbutted)  continue;

			if (face.isSkipped)
			{
				// Face is owned by donor zone, skip it here -> cut nodes layer
				switch (f)
				{
				case faceXmin:
					zne.nx -= 1;
					cgZne.ie1 = (zne.ie -= 1);
					cgZne.is1 -= 1;
					break;
				case faceXmax:
					zne.nx -= 1;
					zne.ie -= 1;
					break;
				case faceYmin:
					zne.ny -= 1;
					cgZne.je1 = (zne.je -= 1);
					cgZne.js1 -= 1;
					break;
				case faceYmax:
					zne.ny -= 1;
					zne.je -= 1;
					break;
				}
			}
		}
	} // loop through zones

	



	  //
	  // (4) Generate global indices in all zones
	  // 
	TZoneGlobIndices** znesGlobIndices = new TZoneGlobIndices*[nZones];

	int gIdx = 0;
	for (int b = 0; b < nZones; ++b)
	{
		TZone& zne = Zones[b];
		znesGlobIndices[b] = new TZoneGlobIndices(zne);

		zne.nGlobStart = gIdx;

		// Real nodes count (without ghosts)
		const int
			nx0 = zne.ie - zne.is + 1,
			ny0 = zne.je - zne.js + 1;

		gIdx += nx0*ny0;
	}


	//
	// (3) Parse connectivity data and fill-in ghost nodes indices
	//
	for (int stage = 0; stage < 2; ++stage) {
		for (int b = 0; b < nZones; ++b)
		{
			TZone& zne = Zones[b];
			TcgnsZone& cgZne = ctx.cgZones[b];
			TZoneGlobIndices& zneGlobInd = *znesGlobIndices[b];

			// Loop through face patches
			for (int p = 0; p < cgZne.nPatches; ++p)
			{
				const TcgnsZone::TFacePatch& cgPatch = cgZne.Patches[p];

				// Start indices of the original patch (assuming no layers were skipped) in working numbering (with ghosts)
				const int ips = cgPatch.is0 + cgZne.is1 - 1;
				const int jps = cgPatch.js0 + cgZne.js1 - 1;

				//
				// Start & End indices of ghost layer
				int igs, ige, jgs, jge;

				// Indices range of the patch in working numbering assuming no layers were skipped
				igs = cgPatch.is0 + cgZne.is1 - 1;
				ige = cgPatch.ie0 + cgZne.is1 - 1;
				jgs = cgPatch.js0 + cgZne.js1 - 1;
				jge = cgPatch.je0 + cgZne.js1 - 1;

				if (stage > 0)
				{
					// Extend patch nodes by ghosts if patch edge coinside with the face edge
					if (cgPatch.is0 == 1)                            igs = 1;
					if (cgPatch.ie0 == (cgZne.ie1 - cgZne.is1 + 1))  ige = zne.nx;

					if (cgPatch.js0 == 1)                            jgs = 1;
					if (cgPatch.je0 == (cgZne.je1 - cgZne.js1 + 1))  jge = zne.ny;
				}

				// Nodes of the ghost layer extruded from the patch
				switch (cgPatch.posFace)
				{
				case faceXmin:  igs = 1;          ige = zne.is - 1;   break;
				case faceXmax:  igs = zne.ie + 1;   ige = zne.nx;     break;
				case faceYmin:  jgs = 1;          jge = zne.js - 1;   break;
				case faceYmax:  jgs = zne.je + 1;   jge = zne.ny;     break;
				}

				const TZone& zneDnr = Zones[cgPatch.nDnrZne];
				const TcgnsZone& cgZneDnr = ctx.cgZones[cgPatch.nDnrZne];
				TZoneGlobIndices& zneDnrGlobInd = *znesGlobIndices[cgPatch.nDnrZne];

				// Start indices of the original donor patch (assuming no layers were skipped) in working numbering
				const int ids = cgPatch.dnr_is0 + cgZneDnr.is1 - 1;
				const int jds = cgPatch.dnr_js0 + cgZneDnr.js1 - 1;

				for (int j = jgs; j <= jge; ++j) {
					for (int i = igs; i <= ige; ++i)
					{
						// Donor zone indices
						int id = (i - ips)*cgPatch.matTrans[0] + (j - jps)*cgPatch.matTrans[1] + ids;
						int jd = (i - ips)*cgPatch.matTrans[2] + (j - jps)*cgPatch.matTrans[3] + jds;

						if (0 == stage)
						{
							if (id < cgZneDnr.is1 || id > cgZneDnr.ie1 ||
								jd < cgZneDnr.js1 || jd > cgZneDnr.je1)
							{
								wxLogError(
									_("Ghost node (%d,%d) of zone '%s' is mapped \
to nonexistent node (%d,%d) of zone '%s'. Check indices orientation!"),
i, j, wxString::FromAscii(zne.szName).c_str(),
id, jd, wxString::FromAscii(zneDnr.szName).c_str()
);
								return false;
							}

							if (id < zneDnr.is || id > zneDnr.ie ||
								jd < zneDnr.js || jd > zneDnr.je)
							{
								continue;  // process on next stage
							}

							zneGlobInd(i, j) = zneDnr.globRealInd(id, jd);
						}
						else if (1 == stage)  // resolve cross references (nodes around corners)
						{
							if (zneGlobInd(i, j) != -1)  continue;

							if (id < 1 || id > zneDnr.nx ||
								jd < 1 || jd > zneDnr.ny)
							{
								// ghost data unavailable
								continue;
							}

							zneGlobInd(i, j) = zneDnrGlobInd(id, jd);
						}
					}
				} // for i,j
			}  // loop faces
		}
	} // loop zones, stages


	  //
	  // Fill-in globIndices[] of zones in current MPI-rank
	  //
	for (int b = bs; b <= be; ++b)
	{
		TZone& zne = Zones[b];
		TZoneGlobIndices& zneGlobInd = *znesGlobIndices[b];

		try {
			zne.globIndices = new int[zne.nx * zne.ny];
		}
		catch (std::bad_alloc& e)
		{
			printf(
				"[%d] %s(): Memory allocation for TZone::globIndices[] in zone #%d/%d FAILED. Asked for %.3f MB.\n",
				G_State.mpiRank, __FUNCTION__,
				b + 1, nZones,
				zne.nx*zne.ny * (float)sizeof(int) / (1024.*1024.)
			);
			return false;
		}

		// Fill-in global indices of real nodes
		for (int j = 1; j <= zne.ny; ++j) {
			for (int i = 1; i <= zne.nx; ++i)
			{
				zne.globInd(i, j) = zneGlobInd(i, j);
			}
		}
	}

	for (int b = 0; b < nZones; ++b)  delete znesGlobIndices[b];
	delete[] znesGlobIndices;

	return true;
}
//-----------------------------------------------------------------------------

const TcgnsZone::TFacePatch& t_MFCGNS2D::get_face_patch(const t_ZoneNode& a_znode) const {

	int b = a_znode.iZone - 1;

	const int nd_i = a_znode.iNode.i;
	const int nd_j = a_znode.iNode.j;
	const TcgnsZone& cgZne = cgCtx.cgZones[b];

	// Start & End indices of ghost layer
	int ips, ipe, jps, jpe;

	for (int ip = 0; ip < cgZne.nPatches; ip++) {

		TcgnsZone::TFacePatch& cgPatch = cgZne.Patches[ip];

		// important: when on edge, skip edges of patches from adjacent faces
		if (cgPatch.posFace != a_znode.iFacePos) continue;

		bool ok = true;

		// Indices range of the patch in working numbering assuming no layers were skipped
		ips = cgPatch.is0 + cgZne.is1 - 1;
		ipe = cgPatch.ie0 + cgZne.is1 - 1;
		jps = cgPatch.js0 + cgZne.js1 - 1;
		jpe = cgPatch.je0 + cgZne.js1 - 1;

		ok = ok && (nd_i >= ips) && (nd_i <= ipe);
		ok = ok && (nd_j >= jps) && (nd_j <= jpe);

		if (ok) {

			return cgPatch;

		}

	}

	wxLogMessage(_T("t_MFCGNS2D: get_face_patch error: provided znode is not on zone face patches..."));

	return TcgnsZone::TFacePatch();

}

// get donor zone-node index of the specified zone-node index
// we have a_znode, shifting with di, dj, dk (not exceeding size of ghost layer!)
// then we get znode index for the specified znode
// (Zone, i, j, k) <=> (ZoneDonor, i_donor, j_donor, k_donor)
t_ZoneNode t_MFCGNS2D::get_abutted_znode(
	const t_ZoneNode& a_znode, const int di, const int dj, const int dk) const{

		if (dk!=0){
			wxLogError(_T("MFCGNS2D: get_abutted_znode error: k-shift is not allowed (2D configuration)"));
		}

		int b = a_znode.iZone - 1;

		TZone& zne = Zones[b];
		TcgnsZone& cgZne = cgCtx.cgZones[b];

		if (a_znode.iFacePos == faceNone)
			wxLogError(_T("In get_abutted_znode: ifacepos not specified for input znode param"));

		const TcgnsZone::TFacePatch& cgPatch = get_face_patch(a_znode);

		TZone& zneDonor = Zones[cgPatch.nDnrZne];
		const TcgnsZone& cgZneDnr = cgCtx.cgZones[cgPatch.nDnrZne];

		// Start indices of the original patch (assuming no layers were skipped) 
		// in working numbering (with ghosts)
		const int ips = cgPatch.is0 + cgZne.is1 - 1;
		const int jps = cgPatch.js0 + cgZne.js1 - 1;

		const int i = a_znode.iNode.i + di;
		const int j = a_znode.iNode.j + dj;
		const int k = a_znode.iNode.k + dk;

		const int ids = cgPatch.dnr_is0 + cgZneDnr.is1 - 1;
		const int jds = cgPatch.dnr_js0 + cgZneDnr.js1 - 1;

		// Donor zone indices
		int id = (i - ips)*cgPatch.matTrans[0] + (j - jps)*cgPatch.matTrans[1] + ids;
		int jd = (i - ips)*cgPatch.matTrans[2] + (j - jps)*cgPatch.matTrans[3] + jds;

		if( id > zneDonor.nx || id < 1 ||
			jd > zneDonor.ny || jd < 1    )
		{
			wxLogError(_("Get abutted znode 2D : Failed to get corrrect abutted znode index"));
		}

		t_ZoneNode ret;

		// 1-based zone id
		ret.iZone = cgPatch.nDnrZne + 1;

		// k is unchanged
		ret.iNode = t_BlkInd(id, jd, a_znode.iNode.k);

		// face position unknown after shift
		ret.iFacePos = faceNone;

		return ret;
};


bool t_MFCGNS2D::_parseBCData2DfromCGNS( TcgnsContext& ctx )
{

	assert(ctx.cgZones);

	for (int iZone = 1; iZone <= nZones; ++iZone)
	{
		TZone& zne = Zones[iZone - 1];
		TcgnsZone& cgZne = ctx.cgZones[iZone - 1];

		// Original zone size (not accounting added ghosts & skipped layers)
		const int
			nx0 = cgZne.ie1 - cgZne.is1 + 1,
			ny0 = cgZne.je1 - cgZne.js1 + 1;

		// Number of BCs in the Zone
		int nBCs = 0;  cg_nbocos(ctx.fileID, ctx.iBase, iZone, &nBCs);

		for (int iBC = 1; iBC <= nBCs; ++iBC)
		{
			CG_BCType_t iBCtype;

			CG_PointSetType_t pntSet;
			cgsize_t nPnts = 0; // number of points defining the BC region

								// Normals to the BC patch
			int iNorm[2]; // normal as index vector (computational coords)
			cgsize_t normListSize;  CG_DataType_t normDataType; // normals in phys coords

			int nDatasets = 0; // number of datasets with additional info for the BC

			char szName[33];

			cg_boco_info(ctx.fileID, ctx.iBase, iZone, iBC,
				szName, &iBCtype,
				&pntSet, &nPnts,
				iNorm, &normListSize, &normDataType,
				&nDatasets
			);
			if (pntSet != CG_PointRange && nPnts != 2)
			{
				wxLogError(
					_("Boundary condition patch '%s'#%d of zone '%s'#%d isn't defined as point range"),
					wxString::FromAscii(szName).c_str(), iBC,
					wxString::FromAscii(zne.szName).c_str(), iZone
				);
				return false;
			}

			// Family name
			cg_goto(ctx.fileID, ctx.iBase, "Zone_t", iZone, "ZoneBC", 0, "BC_t", iBC, NULL);
			if (cg_famname_read(szName) == CG_OK)
			{
				// Read BC-family "Fam_Descr_Name" generated by Pointwise 16.03
				if (cg_goto(ctx.fileID, ctx.iBase, szName, 0, NULL) == CG_OK)
				{
					char szDescrName[33] = "";
					char* szDescrText = NULL;
					if (cg_descriptor_read(1, szDescrName, &szDescrText) == CG_OK)
					{
						if (strcmp(szDescrName, "Fam_Descr_Name") == 0)
						{
							strcpy(szName, szDescrText);
							cg_free(szDescrText);
						}
					}
				}
			}

			// Read BC patch point range
			cgsize_t idxRange[4];   // Imin, Jmin, Imax, Jmax
			cg_boco_read(ctx.fileID, ctx.iBase, iZone, iBC, idxRange, NULL);

			const long
				&is = idxRange[0], &js = idxRange[1],
				&ie = idxRange[2], &je = idxRange[3];

			// Check that BC patch cover whole face
			bool ok = true;
			ok &= (is == 1 && ie == 1) || (is == nx0 && ie == nx0) || (is == 1 && ie == nx0);
			ok &= (js == 1 && je == 1) || (js == ny0 && je == ny0) || (js == 1 && je == ny0);
			if (!ok)
			{
				wxLogError(
					_("Boundary condition patch '%s'#%d of zone '%s'#%d doesn't cover whole face"),
					wxString::FromAscii(szName).c_str(), iBC,
					wxString::FromAscii(zne.szName).c_str(), iZone
				);
				return false;
			}

			// Detect BC face
			TZoneFacePos posFace = faceNone;
			if (is == ie)
				posFace = (is == 1) ? faceXmin : faceXmax;
			else if (js == je)
				posFace = (js == 1) ? faceYmin : faceYmax;

			strcpy(zne.Faces[posFace].szBCFamName, szName);
		} // for iBC

		  //
		  // Fill-in missed BCs
		  //
		for (int f = 0; f < 4; ++f)
		{
			TZoneFace& bF = zne.Faces[f];
			if (bF.bcType != bcTypeH || bF.szBCFamName[0])  continue;

			const char* names[] = { "Imin", "Imax", "Jmin", "Jmax" };
			sprintf(bF.szBCFamName, "%s-%s", zne.szName, names[f]);

			wxLogWarning(
				_("Zone '%s' has face w/o BC; name '%s' was assigned"),
				wxString::FromAscii(zne.szName).c_str(),
				wxString::FromAscii(bF.szBCFamName).c_str()
			);
		}

	} // for iZone

	return true;
}


//-----------------------------------------------------------------------------