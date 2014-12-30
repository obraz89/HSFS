#include "stdafx.h"
#include "cgns_structs.h"

#include <cgnslib.h>

using namespace mf::cg;

// Computational domain composed of blocks with own grid and field

static int g_time = 0.0;

bool mf::cg::TDomain::initField(const wxString& a_fldFName)
{
	g_time = 0.0;
	//const int &nu = G_Domain.nu;

	//
	// Memory allocation for field data
	//

	// Space to store field data of zones local to current MPI rank
	for( int b = bs; b <= be; ++b )
	{
		TZone& zne = Zones[b];
		const int &nx = zne.nx,  &ny = zne.ny,  &nz = zne.nz;

		// Current (n+1) time layer
		zne.U = new double[nu*nx*ny*nz];

		// Previous (n) time layer
		// TODO: can i set this to NULL as well as Unm1
		zne.Un = new double[nu*nx*ny*nz];

		// Pre-previous (n-1) time layer
		if( /*G_Params.timeOrder==2*/ false )
			zne.Unm1 = new double[nu*nx*ny*nz];
		else
			zne.Unm1 = NULL;

		if( (! zne.U) || (! zne.Un) )
		{
			wxLogError(_("Can't allocate memory for solution"));
			return false;
		}
	}

	// Additional space to store whole solution on the root MPI rank
	// TODO: all already allocated when we have 1 core ?
	/*
	if( G_State.mpiRank == 0 )
	{
		for( int z = 0; z < G_Domain.nZones; ++z )
		{
			TZone& zne = G_Domain.Zones[z];
			const int &nx = zne.nx,  &ny = zne.ny,  &nz = zne.nz;

			if( zne.U ) continue;  // already allocated

			// Current (n+1) time layer
			zne.U = new double[nu*nx*ny*nz];
			if( ! zne.U )
			{
				wxLogError(_("Can't allocate memory for solution"));
				return false;
			}
		}
	}
	*/


	//
	// Load initial estimation from file
	//
	{
		wxLogMessage(
			_("* Reading initial field from file '%s'..."),
			a_fldFName.c_str()
		);

		if( ! loadField(a_fldFName, false) )
			return false;
	}

	// Set initial time in time-stepping
	//G_Plugins.get_discret_caps().setParR(5/*tau*/, g_time);
	
	return true;
}
//-----------------------------------------------------------------------------

/**
 *  Loads field data from the file
 *  
 *  @param[in] fileName
 *  @param[in] bIsPrev - the field is for previous time step
 *  
 *  @return true if succeeded, false if failed
**/  
bool mf::cg::TDomain::loadField( const wxString& fileName, bool bIsPrev )
{
	//const int &nu = G_Domain.nu;

	if( fileName.EndsWith(_T(".cgns")) )
	{
		// CGNS file format
		//
		int f = -1;
		if( cg_open( fileName.mb_str(),CG_MODE_READ, &f ) != CG_OK )
		{
			wxLogError(
				_("Can't open field file '%s' for reading. Possibly \
unaccessible file, broken internal links (e.g. to grid) \
or unsupported CGNS version"),
				fileName.c_str()
			);
			return false;
		}

		// Assume only one base in the file
		const int iBase = 1;

		// Space dimensions
		int dimCell = 0, dimPhys = 0;
		char szName[33];
		cg_base_read(f,iBase,  szName,&dimCell,&dimPhys);

		if( dimCell != nDim )
		{
			wxLogError(_("CGNS: Inconsistent space dimensions"));
			return false;
		}

		// Number of zones (aka blocks)
		int nBlk = 0;  cg_nzones(f,iBase, &nBlk);
		if( nBlk != nZones )
		{
			wxLogError( _("CGNS: Inconsistent number of blocks") );
			return false;
		}

		// Loop through blocks
		for( int b = bs; b <= be; ++b )
		{
			int iZone = b + 1;
			TZone& blk = Zones[b];
			double* dstU = (bIsPrev) ? blk.Unm1 : blk.U;

			// Block size without ghosts
			const int  
				nx = blk.ie - blk.is + 1,
				ny = blk.je - blk.js + 1,
				nz = blk.ke - blk.ks + 1;

			// Data of the block excluding ghosts!!!
			TDims dims;  double *newGrid = NULL, *newField = NULL;

			if( ! readBlockFromCGNS(f,iZone,  dims, &newGrid, &newField) )
				return false;

			if( dims.numX==nx && dims.numY==ny && dims.numZ==nz
				&& /*g_genOpts.remapField==false*/ true       )
			{
				// Loaded and current grids are the same -> copy field
				const double* srcU = newField;
				for( int k = blk.ks;  k <= blk.ke;  ++k )
				for( int j = blk.js;  j <= blk.je;  ++j )
				for( int i = blk.is;  i <= blk.ie;  ++i )
				{
					double* U = dstU + nu*( blk.absIdx(i,j,k) - 1 );

					memcpy( U, srcU, nu*sizeof(double) );
					srcU += nu;
				}
			}
			else
			{
				wxLogError(_T("Loaded and current grids do not match, wrong grid file?"));
				ssuTHROW(t_GenException, _T("in LoadField: Loaded and current grids do not match"));
				// Loaded and current grids are different -> remap field
				//remapField(b, bIsPrev, dims, newGrid, newField);
			}

			delete[] newGrid, newField;
		} // Loop through blocks

		
		// Get time
		/*
		if( g_genOpts.initTimeFromFile )
		do
		{
			int nTmSteps = 0;  cg_biter_read(f,iBase, szName, &nTmSteps);
			if( nTmSteps < 1 ) break;

			if( cg_goto(f,iBase,"BaseIterativeData_t",1, NULL) != CG_OK )
				break;

			int nArrs = 0;  cg_narrays(&nArrs);
			if( nArrs < 1 )  break;

			int arrIdx;
			for( arrIdx = 1; arrIdx <= nArrs; ++arrIdx )
			{
				DataType_t type;  int dim;  int len[3];
				cg_array_info(arrIdx, szName, &type, &dim, len);
				if( strcmp( szName, "TimeValues" ) == 0 )
					break;
			}
			if( arrIdx <= nArrs )
				cg_array_read_as(arrIdx, RealDouble, &g_time);
		}
		while(false);
		*/

		cg_close(f);
	}
	else
	{
		// Legacy file format (.ttl, .hsx)
		// 
		wxLogError(_T("Load Fld Error: Not a CGNS extension[.cgns]"));
		ssuTHROW(t_GenException, _T("In LoadField: Not a CGNS format"));
	}

	return true;
}

/**
 *  Read block data from the opened CGNS file
 *  
 *  @param[in] fileID - ID of the opened CGNS file
 *  @param[in] iZone - 1-based zone (block) number in the file
 *  @param[out] dims - dimensions of the block, without ghosts
 *  @param[out] grid - node coordinates of the loaded block's grid
 *  @param[out] field - loaded field data
 *  
 *  @return true if succeeded, false if failed
**/  
bool mf::cg::TDomain::readBlockFromCGNS(
	int fileID, int iZone,
	mf::cg::TDims& dims, double** grid, double** field)
{
	int r = CG_OK;
	const int iBase = 1;  // NB: Assume only one base

	// Get zone size and name
	char szZone[33];
	cgsize_t isize[9];  // NVertexI, NVertexJ, NVertexK,
	               // NCellI, NCellJ, NCellK,
	               // NBoundVertexI, NBoundVertexJ, NBoundVertexK
	if( cg_zone_read(fileID,iBase,iZone,  szZone,isize) != CG_OK )
	{
		wxLogError( _("Can't read CGNS zone #%d"), iZone );
		return false;
	}

	// Block size without ghosts
	unsigned int &nx = dims.numX;  nx = isize[0];
	unsigned int &ny = dims.numY;  ny = isize[1];
	unsigned int &nz = dims.numZ;  nz = (nDim==2) ? 1 : isize[2];

	// Indexes ranges
	cgsize_t irmin[3] = {1, 1, 1};
	cgsize_t irmax[3] = {nx, ny, nz};

	//
	// Read grid coordinates
	// 
	*grid = new double[nDim * nx*ny*nz];
	{
		double* x = new double[nx*ny*nz];
		r |= cg_coord_read(fileID,iBase,iZone,"CoordinateX",RealDouble,irmin,irmax, x);
		if( r != CG_OK ){

			wxLogError(_T("cg_ccord_read error"));
			wxLogError(_T("%s"), wxString::FromAscii(cg_get_error()).c_str());

		}  

		double* y = new double[nx*ny*nz];
		r |= cg_coord_read(fileID,iBase,iZone,"CoordinateY",RealDouble,irmin,irmax, y);

		double* z = NULL;
		if( nDim == 3 )
		{
			z = new double[nx*ny*nz];
			r |= cg_coord_read(fileID,iBase,iZone,"CoordinateZ",RealDouble,irmin,irmax, z);
		}

		if( r != CG_OK )
		{
			wxLogError( _("Can't read coordinates from zone '%s'(#%d)"),
				wxString::FromAscii(szZone).c_str(), iZone
			);
			//return false;
		}

		for( int k=0; k<nz; ++k )
		for( int j=0; j<ny; ++j )
		for( int i=0; i<nx; ++i )
		{
			int n = i + j*nx + k*nx*ny;
			double* G = *grid + nDim * n;

			G[0] = x[n];
			G[1] = y[n];
			if( nDim==3 )  G[2] = z[n];
		}

		delete[] x, y, z;
	}


	//
	// Get solution info
	// FIXME: flow assumed existing
	// 
	int iFlow = 1;  char szFlow[33];  GridLocation_t loc;
	r = cg_sol_info(fileID,iBase,iZone,iFlow,  szFlow,&loc);
	if( r != CG_OK )
	{
		wxLogError( _("Can't read flow info from zone '%s'(#%d)"),
			wxString::FromAscii(szZone).c_str(), iZone
		);
		return false;
	}

	if( loc != Vertex )
	{
		wxLogError(_("CGNS: GridLocation must be Vertex"));
		return false;
	}


	//
	// Load needed functions
	//
	*field = new double[nu * nx*ny*nz];
	{
		dims.numFunc = 0;
		double* FF = new double[nx*ny*nz];

		for( int iFun = 0; iFun < nu; ++iFun )
		{
			const char* name = G_vecCGNSFuncNames[iFun].c_str();

			r = cg_field_read(fileID,iBase,iZone,iFlow,(char*)name,RealDouble,irmin,irmax, FF);
			if( r != CG_OK )
			{
				wxLogError(_T("%s"), wxString::FromAscii(cg_get_error()).c_str());
				wxLogError(
					_("CGNS: Can't read '%s' data from zone '%s'(#%d)"),
					wxString::FromAscii(name).c_str(),
					wxString::FromAscii(szZone).c_str(), iZone
				);
				return false;
			}

			for( int k=0; k<nz; ++k )
			for( int j=0; j<ny; ++j )
			for( int i=0; i<nx; ++i )
			{
				int n = i + j*nx + k*nx*ny;
				int pos = nu * n;

				(*field)[pos + iFun] = FF[n];
			}

			dims.numFunc++;
		}
		delete[] FF;
	}

	return true;
}

// free fld and grd mem - shared part

mf::cg::TDomain::~TDomain(){
	delete[] Zones;
}

bool mf::cg::TDomain::_is_face_of_bcwall_type(const char* faceBCFamName) const{

	std::vector<std::string>::const_iterator iter = _vecBCWallNames.begin();

	std::string facename_str(faceBCFamName);

	for (;iter<_vecBCWallNames.end(); iter++){

		if (*iter==facename_str) 
			return true;

	}

	return false;

};

void TDomain::get_rec(const TZone& zone, const t_BlkInd& ind, mf::t_Rec& rec) const{
	get_rec(zone, ind.i, ind.j, ind.k, rec);
}

void TDomain::get_rec(const t_ZoneNode& znode, mf::t_Rec& rec) const{
	const TZone& blk = Zones[znode.iZone-1];
	get_rec(blk, znode.iNode.i, znode.iNode.j, znode.iNode.k, rec);
}
//IMPORTANT TODO: check and verify {written at midnight!!!} 
void TDomain::extract_profile_data(const mf::t_GeomPoint& xyz, 
	const mf::t_ProfDataCfg& prdata_cfg, std::vector<mf::t_Rec>& data) const{

		double bl_thick;
		t_ZoneNode surf_znode, outer_znode;

		_calc_bl_thick(xyz, bl_thick, surf_znode, outer_znode);

		double total_thick = bl_thick * prdata_cfg.ThickCoef;

		int total_nodes = 0;
//////

		double eta = 0.0;
		t_GeomPoint cur_xyz, surf_xyz;
		t_Vec3Dbl rvec;

		t_ZoneNode cur_znode = surf_znode;
		mf::t_Rec cur_rec;
		get_rec(surf_znode, cur_rec);
		surf_xyz.set(cur_rec);

		int*pd = NULL;
		int dd = 0;
		switch (surf_znode.iFacePos)
		{
			case (faceXmin):
				pd = &(cur_znode.iNode.i);
				dd = 1;
				break;
			case (faceXmax):
				pd = &(cur_znode.iNode.i);
				dd = -1;
				break;
			case (faceYmin):
				pd = &(cur_znode.iNode.j);
				dd = 1;
				break;
			case (faceYmax):
				pd = &(cur_znode.iNode.j);
				dd = -1;
				break;
			case (faceZmin):
				pd = &(cur_znode.iNode.k);
				dd = 1;
				break;
			case (faceZmax):
				pd = &(cur_znode.iNode.k);
				dd = -1;
				break;

		}
		do 
		{
			get_rec(cur_znode, cur_rec);
			cur_xyz.set(cur_rec);
			matrix::base::minus<double, double>(cur_xyz, surf_xyz, rvec);
			eta = rvec.norm();
			*pd+=dd;total_nodes++;

		} while (eta<total_thick);

		data.resize(total_nodes);

		for (int d=0; d<total_nodes; d++){
			cur_znode = surf_znode;
			*pd+=d*dd;
			get_rec(cur_znode, data[d]);
		}
}
