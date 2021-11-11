#include "stdafx.h"
#include "cgns_structs.h"

#include <cgnslib.h>

#include "io_helpers.h"

#include <fstream>

using namespace mf::cg;

#if CG_BUILD_SCOPE
  #define CG_MY_RealDouble CG_RealDouble
  #define CG_MY_GridLocation_t CG_GridLocation_t
  #define CG_MY_Vertex CG_Vertex
#else
  #define CG_MY_RealDouble RealDouble
  #define CG_MY_GridLocation_t GridLocation_t
  #define CG_MY_Vertex Vertex
#endif

#define OPT_FACE_POS_XMIN _T("Xmin")
#define OPT_FACE_POS_XMAX _T("Xmax")
#define OPT_FACE_POS_YMIN _T("Ymin")
#define OPT_FACE_POS_YMAX _T("Ymax")
#define OPT_FACE_POS_ZMIN _T("Zmin")
#define OPT_FACE_POS_ZMAX _T("Zmax")

static const char* g_cgCoordNames[] = { "CoordinateX", "CoordinateY", "CoordinateZ" };

mf::t_DomainCGNSParams::t_DomainCGNSParams() :t_FldParams() {

	FACE_POS_StART_STR.clear();

	FACE_POS_StART_STR.insert(std::make_pair(OPT_FACE_POS_XMIN, mf::cg::TZoneFacePos::faceXmin));
	FACE_POS_StART_STR.insert(std::make_pair(OPT_FACE_POS_XMAX, mf::cg::TZoneFacePos::faceXmax));
	FACE_POS_StART_STR.insert(std::make_pair(OPT_FACE_POS_YMIN, mf::cg::TZoneFacePos::faceYmin));
	FACE_POS_StART_STR.insert(std::make_pair(OPT_FACE_POS_YMAX, mf::cg::TZoneFacePos::faceYmax));
	FACE_POS_StART_STR.insert(std::make_pair(OPT_FACE_POS_ZMIN, mf::cg::TZoneFacePos::faceZmin));
	FACE_POS_StART_STR.insert(std::make_pair(OPT_FACE_POS_ZMAX, mf::cg::TZoneFacePos::faceZmax));
}

// full or truncated (in case surf_znode is internal but carries face_pos of starting face) 
// grid line from a starting znode to its opposite face position;
// end node is always on real zone surface
t_ZoneGrdLine::t_ZoneGrdLine(const TDomain& a_dom, const t_ZoneNode& surf_znode){

	iZone = surf_znode.iZone;

	zne = &a_dom.getZone(iZone);

	assert(surf_znode.iNode.i >= zne->is && surf_znode.iNode.i <=zne->ie);
	assert(surf_znode.iNode.j >= zne->js && surf_znode.iNode.j <=zne->je);
	// no assert for k because of 2D fields (we do 2D to 3D pseudo conversion)

	ind_s = surf_znode.iNode;

	ind_e = ind_s;

	di = dj = dk = 0;

	switch(surf_znode.iFacePos){
		case(faceXmin):
			di = 1;
			ind_s.i = surf_znode.iNode.i;
			ind_e.i = zne->ie;
			idx_change = I_LINE;
			break;
		case(faceXmax):
			di = -1;
			ind_s.i = surf_znode.iNode.i;
			ind_e.i = zne->is;
			idx_change = I_LINE;
			break;
		case(faceYmin):
			dj = 1;
			ind_s.j = surf_znode.iNode.j;
			ind_e.j = zne->je;
			idx_change = J_LINE;
			break;
		case(faceYmax):
			dj = -1;
			ind_s.j = surf_znode.iNode.j;
			ind_e.j = zne->js;
			idx_change = J_LINE;
			break;
		case(faceZmin):
			dk = 1;
			a_dom.get_k_range(surf_znode.iZone, ind_s.k, ind_e.k);
			ind_s.k = surf_znode.iNode.k;
			idx_change = K_LINE;
			break;
		case(faceZmax):
			dk=-1;
			a_dom.get_k_range(surf_znode.iZone, ind_e.k, ind_s.k);
			ind_s.k = surf_znode.iNode.k;
			idx_change = K_LINE;
			break;
		default:
			wxLogError(_T("ZoneGrdLine Error: start point is not on zone surface"));
			break;
	}

	switch (idx_change)
	{
	case I_LINE:
		nnodes = di*(ind_e.i - ind_s.i) + 1;
		break;
	case J_LINE:
		nnodes = dj*(ind_e.j - ind_s.j) + 1;
		break;
	case K_LINE:
		nnodes = dk*(ind_e.k - ind_s.k) + 1;
		break;
	}
}

void t_DomainGrdLine::_add(const t_ZoneGrdLine& grd_line ){

	for (int p=0; p<grd_line.nnodes; p++){
		znodes.push_back(
			t_ZoneNode(grd_line.iZone, 
						t_BlkInd(grd_line.ind_s, 
							p*grd_line.di, p*grd_line.dj, p*grd_line.dk)));
	}

}

void t_DomainGrdLine::init(const t_ZoneNode& a_face_znode){

	t_ZoneNode start_znode = a_face_znode;

	TZoneFacePos cur_face_pos;

	do 
	{

		t_ZoneGrdLine grd_line(dom, start_znode);

		_add(grd_line);

		// debug, dump at every iter to get last line before crash

		dump("output/dom_grdline.dat");

		// move on to the abutted zone if final node is on zone interface

		cur_face_pos = TZoneFace::get_brick_opposite_face(start_znode.iFacePos);

		const TZone& zne = dom.getZone(start_znode.iZone);
		//const TcgnsZone& cgZne = dom.getCgZone(start_znode.iZone);

		if (zne.Faces[cur_face_pos].bcType!=TBlockBCType::bcAbutted)
			break;

		t_ZoneNode znode_face(grd_line.iZone, grd_line.ind_e, cur_face_pos);
		
		// donor znode on (f)ace
		t_ZoneNode znode_dnr_f = dom.get_abutted_znode(znode_face, 0, 0, 0);

		// debug, checking that the nodes coincide
		/*
		const TZone& zneDnr = dom.getZone(znode_dnr_f.iZone);
		const TcgnsZone& cgZneDnr = dom.getCgZone(znode_dnr_f.iZone);

		mf::t_Rec rec1, rec2;

		std::ofstream ofstr("output/znode_face_check.dat", std::ios::app);
		dom.get_rec(znode_face, rec1);
		ofstr << "znode    :" << znode_face.iZone <<"\t"
			<<znode_face.iNode.i<<"\t"<<znode_face.iNode.j<<"\t"<<znode_face.iNode.k<<"\t"
			<< rec1.x << "\t" << rec1.y << "\t" << rec1.z << "\n";
		dom.get_rec(znode_dnr_f, rec2);
		ofstr << "znode_dnr:" << znode_dnr_f.iZone << "\t"
			<< znode_dnr_f.iNode.i << "\t" << znode_dnr_f.iNode.j << "\t" << znode_dnr_f.iNode.k << "\t" 
			<< rec2.x << "\t" << rec2.y << "\t" << rec2.z << "\n";
		ofstr.flush(); ofstr.close();
		*/
		
		// first (i)nterior donor znode
		t_ZoneNode znode_dnr_i = 
			dom.get_abutted_znode(znode_face, grd_line.di, grd_line.dj, grd_line.dk);

		// now we have new di, dj, dk for a donor zone

		int dnr_di = znode_dnr_i.iNode.i - znode_dnr_f.iNode.i;
		int dnr_dj = znode_dnr_i.iNode.j - znode_dnr_f.iNode.j;
		int dnr_dk = znode_dnr_i.iNode.k - znode_dnr_f.iNode.k;

		// check that only one of di, dj, dk is nonzero and equals +1 or -1

		int sum = dnr_di*dnr_di + dnr_dj*dnr_dj + dnr_dk*dnr_dk;

		if (sum!=1){

			wxLogError(_T("Domain Grid Line: failed to get donor zne starting indexes"));
			break;

		};


		// get direction of dnr zone gridline:

		TZoneFacePos new_face_pos = faceNone;

		if (dnr_di==1) new_face_pos = faceXmin;
		if (dnr_di==-1) new_face_pos = faceXmax;

		if (dnr_dj==1) new_face_pos = faceYmin;
		if (dnr_dj==-1) new_face_pos = faceYmax;

		if (dnr_dk==1) new_face_pos = faceZmin;
		if (dnr_dk==-1) new_face_pos = faceZmax;

		start_znode.iFacePos = new_face_pos;
		start_znode.iZone = znode_dnr_i.iZone;

		// to avoid duplication of face nodes start from first interior node;
		// this node does NOT belong to zone surface but carries facePos information
		// to initialize zone grd line
		start_znode.iNode = znode_dnr_i.iNode;

	} while (true);

}

void t_DomainGrdLine::dump(std::string fname) const{

	std::ofstream ofstr(fname);
	ofstr.width(12);
	ofstr.precision(6);

	mf::t_Rec rec;

	ofstr << "iZone\ti\tj\tk\tx\ty\tz\tu\tv\tw\n";

	for (int i= 0; i < znodes.size(); i++) {

		const t_ZoneNode& z = znodes[i];
		dom.get_rec(z, rec);
		ofstr << z.iZone << "\t" << z.iNode.i << "\t" << z.iNode.j << "\t" << z.iNode.k << "\t";
		ofstr << rec.x << "\t" << rec.y << "\t" << rec.z << "\t"<< rec.u << "\t" << rec.v << "\t" << rec.w << "\n";

	}
	ofstr.flush();

}

void t_DomainGrdLine::dump_as_ppoints(std::string fname) const {

	std::ofstream ofstr(fname);

	ofstr.setf(std::ios_base::left, std::ios_base::adjustfield);
	ofstr.setf(std::ios_base::scientific, std::ios_base::floatfield);

	ofstr.width(12);
	ofstr.precision(6);

	mf::t_Rec rec;

	ofstr << znodes.size() << "\n";

	for (int i = 0; i < znodes.size(); i++) {

		const t_ZoneNode& z = znodes[i];
		dom.get_rec(z, rec);
		ofstr <<rec.x << "\t" << rec.y << "\t" << rec.z << "\n";

	}
	ofstr.flush();

}

// Computational domain composed of blocks with own grid and field

static int g_time = 0.0;

mf::cg::TDomain::TDomain() {
	_pGrdLine = new t_GrdLineFullCashed();
	_pProfRaw = new t_ProfileRawCashed();
}


void mf::cg::TDomain::_read_parse_bc_wall_names(const wxString& strBCWallFamNames){

	io_hlp::read_parse_str_array(strBCWallFamNames, _T(','), _vecBCWallNames);

}

void mf::cg::TDomain::_read_parse_func_names(const wxString& str){

	io_hlp::read_parse_str_array(str, _T(','), G_vecCGNSFuncNames);

}



void mf::cg::TDomain::_read_parse_bl_thick_calc_type(const wxString& a_str){

	wxString str = a_str;

	str.Trim(true);str.Trim(false);

	if (str.CmpNoCase(_T("BY_VELO_DERIV"))==0){

		_profile_cfg.BLThickCalcType = mf::t_BLThickCalcType::BLTHICK_BY_VDERIV;
		return;

	}

	if (str.CmpNoCase(_T("BY_ENTHALPY"))==0){

		_profile_cfg.BLThickCalcType = mf::t_BLThickCalcType::BLTHICK_BY_ENTHALPY;
		return;

	}

	if (str.CmpNoCase(_T("FULL_GRIDLINE"))==0){

		_profile_cfg.BLThickCalcType = mf::t_BLThickCalcType::BLTHICK_FULL_GRIDLINE;
		return;

	}

	wxString msg(_T("Unknown BLCalcType value, supported options are: BY_VELO_DERIV, BY_ENTHALPY"));
	wxLogError(msg); ssuGENTHROW(msg);

}

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

	if (fileName.EndsWith(_T(".cgns")))
	{
		// CGNS file format
		//
		int f = -1;
		if (cg_open(fileName.mb_str(), CG_MODE_READ, &f) != CG_OK)
		{
			wxLogError(
				_("Can't open field file '%s' for reading: %s"),
				fileName.c_str(), wxString::FromAscii(cg_get_error()).c_str()
			);
			return false;
		}

		// Assume only one base in the file
		const int iBase = 1;

		// Space dimensions
		int dimCell = 0, dimPhys = 0;
		char szName[33];
		cg_base_read(f, iBase, szName, &dimCell, &dimPhys);

		if (dimCell != nDim)
		{
			wxLogError(_("CGNS: Inconsistent space dimensions"));
			return false;
		}

		// Number of zones (aka blocks)
		int nBlk = 0;  cg_nzones(f, iBase, &nBlk);
		if (nBlk != nZones)
		{
			wxLogError(_("CGNS: Inconsistent number of zones"));
			return false;
		}

		// Loop through zones
		for (int b = bs; b <= be; ++b)
		{
			int iZone = b + 1;
			TZone& zne = Zones[b];
			if (zne.isFrozen && bIsPrev) continue;

			double* dstU = (bIsPrev) ? zne.Unm1 : zne.U;

			// Original zone size
			// (not accounting added ghosts & skipped layers)
			const int is0 = zne.is - (zne.Faces[faceXmin].isSkipped ? 1 : 0);
			const int ie0 = zne.ie +(zne.Faces[faceXmax].isSkipped ? 1 : 0);
			const int js0 = zne.js -(zne.Faces[faceYmin].isSkipped ? 1 : 0);
			const int je0 = zne.je +(zne.Faces[faceYmax].isSkipped ? 1 : 0);
			const int ks0 = zne.ks -(zne.Faces[faceZmin].isSkipped ? 1 : 0);
			const int ke0 = zne.ke +(zne.Faces[faceZmax].isSkipped ? 1 : 0);

			const unsigned int
				nx0 = ie0 - is0 + 1,
				ny0 = je0 - js0 + 1,
				nz0 = ke0 - ks0 + 1;

			// Data of the zone excluding ghosts!!!
			TDims dims;  double *newGrid = NULL, *newField = NULL;
			
			if (!readCGNSZone(f, iBase, iZone, dims, &newGrid, &newField)) {
				wxLogMessage(_T("Error: Failed to read cgns zone iZone=%d"), iZone);
				return false;
			}

			if (dims.numX != nx0 || dims.numY != ny0 || dims.numZ != nz0)
			{
				wxLogError(
					_("Zone #%d has different dimensions in the loaded field (%dx%dx%d) and in the current grid (%dx%dx%d)."),
					dims.numX, dims.numY, dims.numZ,
					nx0, ny0, nz0
				);
				return false;
			}

			// Copy field to internal structs
			// NB: skipped grid layers are also filled to support immediate saving (e.g. slices) without scatternig ghosts
			for (int k = ks0; k <= ke0; ++k) {
				for (int j = js0; j <= je0; ++j) {
					for (int i = is0; i <= ie0; ++i)
					{
						double* U = dstU + nu*(zne.absIdx(i, j, k) - 1);
						double* srcU = newField + nu*((i - is0) + (j - js0)*nx0 + (k - ks0)*nx0*ny0);

						memcpy(U, srcU, nu * sizeof(double));
					}
				}
			}

			delete[] newGrid;   delete[] newField;
		} // Loop through blocks

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
bool mf::cg::TDomain::readCGNSZone(
	const int fileID, const int iBase, const int iZone,
	TDims& dims, double** grid, double** field)
{
	char cgName[33];  // temp string to store CGNS node names

					  // Get zone size and name
	char szZone[33];
	cgsize_t isize[9];  // NVertexI, NVertexJ, NVertexK,
						// NCellI, NCellJ, NCellK,
						// NBoundVertexI, NBoundVertexJ, NBoundVertexK
	if (cg_zone_read(fileID, iBase, iZone, szZone, isize) != CG_OK)
	{
		wxLogError(_("Can't read CGNS zone #%d ( %s )"),
			iZone, wxString::FromAscii(cg_get_error()).c_str());
		return false;
	}

	if (strcmp(szZone, Zones[iZone - 1].szName) != 0)
	{
		wxLogWarning(
			_("Inconsistent zone #%d names: '%s' -> '%s'"), iZone,
			wxString::FromAscii(szZone).c_str(),
			wxString::FromAscii(Zones[iZone - 1].szName).c_str()
		);
	}


	// Block size without ghosts
	unsigned int &nx = dims.numX;  nx = isize[0];
	unsigned int &ny = dims.numY;  ny = isize[1];
	unsigned int &nz = dims.numZ;  nz = (nDim == 2) ? 1 : isize[2];

	// Indexes faces
	cgsize_t irmin[3] = { 1, 1, 1 };
	cgsize_t irmax[3] = { nx, ny, nz };

	// Array for reading CGNS values into
	double* VV = new double[nx*ny*nz];

	//
	// Read grid coordinates
	// 
	*grid = new double[nDim * nx*ny*nz];
	for (int iCoord = 0; iCoord < nDim; ++iCoord)
	{
		const char* name = g_cgCoordNames[iCoord];

		if (cg_coord_read(fileID, iBase, iZone, name, CG_RealDouble, irmin, irmax, VV) != CG_OK)
		{
			wxLogError(
				_("Can't read %s from zone '%s'#%d ( %s )"),
				wxString::FromAscii(name).c_str(),
				wxString::FromAscii(szZone).c_str(), iZone,
				wxString::FromAscii(cg_get_error()).c_str()
			);
			return false;
		}

		// Pack to output array
		for (unsigned int k = 0; k<nz; ++k) {
			for (unsigned int j = 0; j<ny; ++j) {
				for (unsigned int i = 0; i<nx; ++i)
				{
					int n = i + j*nx + k*nx*ny;
					int pos = nDim * n + iCoord;

					(*grid)[pos] = VV[n];
				}
			}
		}
	} // for iCoord


	  //
	  // Get solution info
	  // FIXME: flow assumed existing
	  // 
	int iSol = 1;
	{
		CG_GridLocation_t loc;
		if (cg_sol_info(fileID, iBase, iZone, iSol, cgName, &loc) != CG_OK)
		{
			wxLogError(_("Can't read flow info from zone '%s'#%d ( %s )"),
				wxString::FromAscii(szZone).c_str(), iZone,
				wxString::FromAscii(cg_get_error()).c_str()
			);
			return false;
		}

		if (loc != CG_Vertex)
		{
			wxLogError(_("CGNS: GridLocation must be Vertex"));
			return false;
		}
	}


	//
	// Read functions
	//
	*field = new double[nu * nx*ny*nz];

	dims.numFunc = 0;
	for (int iFun = 0; iFun < nu; ++iFun)
	{
		const char* name = G_vecCGNSFuncNames[iFun].c_str();

		int r = cg_field_read(fileID, iBase, iZone, iSol, name,
			CG_RealDouble, irmin, irmax, VV);

		if (r != CG_OK && r != CG_NODE_NOT_FOUND)
		{
			wxLogError(
				_("Can't read '%s' from zone '%s'#%d ( %s )"),
				wxString::FromAscii(name).c_str(),
				wxString::FromAscii(szZone).c_str(), iZone,
				wxString::FromAscii(cg_get_error()).c_str()
			);
			return false;
		}

		if (r == CG_NODE_NOT_FOUND)
		{
			wxLogWarning(
				_("%s not found in zone '%s'#%d, using free stream values"),
				wxString::FromAscii(name).c_str(),
				wxString::FromAscii(szZone).c_str(), iZone
			);

			for (unsigned int k = 0; k<nz; ++k) {
				for (unsigned int j = 0; j<ny; ++j) {
					for (unsigned int i = 0; i<nx; ++i)
					{
						int n = i + j*nx + k*nx*ny;
						int pos = nu * n + iFun;

						(*field)[pos] = infVals[iFun];
					}
				}
			}

			continue;
		}

		bool checkPositive =
			(strcmp(name, "Pressure") == 0) ||
			(strcmp(name, "Temperature") == 0) ||
			(strstr(name, "TurbulentEnergy") == name); // begins with "Turb..."

													   // Pack to output array
		for (unsigned int k = 0; k<nz; ++k) {
			for (unsigned int j = 0; j<ny; ++j) {
				for (unsigned int i = 0; i<nx; ++i)
				{
					const int n = i + j*nx + k*nx*ny;
					const int pos = nu * n + iFun;

					if (checkPositive && VV[n] < 0)
					{
						wxLogWarning(
							_("Negative value %s=%g @ '%s'(%d,%d,%d). Using absolute value."),
							wxString::FromAscii(name).c_str(), VV[n],
							wxString::FromAscii(szZone).c_str(), i + 1, j + 1, k + 1
						);
						VV[n] *= -1;
					}

					(*field)[pos] = VV[n];
				}
			}
		}

		dims.numFunc++;
	}

	delete[] VV;

	return true;
}

// free fld and grd mem - shared part

mf::cg::TDomain::~TDomain(){
	delete[] Zones;

	delete _pGrdLine;
	delete _pProfRaw;
}

// TODO: make less restrictive comparison (case insensitive maybe, trimming)
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

// calculate derivative of some flow field variable 
// in mathematical reference frame (ksi, eta, dzeta)
// flow var is get by name:(x,y,z,u,v,w,p,t,r)
void TDomain::_calc_scalar_ked_deriv(const t_ZoneNode& znode, char fun_name, t_Vec3Dbl& df_dked) const {

	const int i = znode.iNode.i;
	const int j = znode.iNode.j;
	const int k = znode.iNode.k;

	const t_BlkInd ijk = znode.iNode;

	const TZone& zne = Zones[znode.iZone -1];

	//IMPORTANT TODO: ask Nova about globInd check (see HSFlow grid-data-3d.cpp)
	const bool noBakI = (i == zne.is);//|| (globInd(i - 1, j, k) == -1);
	const bool noBakJ = (j == zne.js);//|| (globInd(i, j - 1, k) == -1);
	const bool noBakK = (k == zne.ks);//|| (globInd(i, j, k - 1) == -1);

	// Can't go forward
	const bool noFwdI = (i == zne.ie);// || (globInd(i + 1, j, k) == -1);
	const bool noFwdJ = (j == zne.je);// || (globInd(i, j + 1, k) == -1);
	const bool noFwdK = (k == zne.ke);// || (globInd(i, j, k + 1) == -1);

	mf::t_Rec r0, r1, r2;
	double df_dksi, df_deta, df_ddze;
	//
	// d?_dksi
	//
	if (noBakI)
	{
		get_rec(zne, i + 0, j, k, r0);
		get_rec(zne, i + 1, j, k, r1);
		get_rec(zne, i + 2, j, k, r2);

		double v0 = r0.get_val(fun_name);
		double v1 = r1.get_val(fun_name);
		double v2 = r2.get_val(fun_name);

		df_dksi = 0.5*(-3.*v0 + 4.*v1 - v2);
	}
	else if (noFwdI)
	{
		get_rec(zne, i - 0, j, k, r0);
		get_rec(zne, i - 1, j, k, r1);
		get_rec(zne, i - 2, j, k, r2);

		double v0 = r0.get_val(fun_name);
		double v1 = r1.get_val(fun_name);
		double v2 = r2.get_val(fun_name);

		df_dksi = 0.5*(3.*v0 - 4.*v1 + v2);
	}
	else
	{

		get_rec(zne, i - 1, j, k, r1);
		get_rec(zne, i + 1, j, k, r2);

		double v1 = r1.get_val(fun_name);
		double v2 = r2.get_val(fun_name);

		df_dksi = 0.5*(v2 - v1);

		// Guard a point of tripple zones intersection
		//assert(fabs(D.x_dksi) + fabs(D.y_dksi) + fabs(D.z_dksi) > 1e-10);
	}
	//
	// d?_deta
	//
	if (noBakJ)
	{
		get_rec(zne, i, j + 0, k, r0);
		get_rec(zne, i, j + 1, k, r1);
		get_rec(zne, i, j + 2, k, r2);

		double v0 = r0.get_val(fun_name);
		double v1 = r1.get_val(fun_name);
		double v2 = r2.get_val(fun_name);

		df_deta = 0.5*(-3.*v0 + 4.*v1 - v2);
	}
	else if (noFwdJ)
	{
		get_rec(zne, i, j - 0, k, r0);
		get_rec(zne, i, j - 1, k, r1);
		get_rec(zne, i, j - 2, k, r2);

		double v0 = r0.get_val(fun_name);
		double v1 = r1.get_val(fun_name);
		double v2 = r2.get_val(fun_name);

		df_deta = 0.5*(3.*v0 - 4.*v1 + v2);
	}
	else
	{
		get_rec(zne, i, j - 1, k, r1);
		get_rec(zne, i, j + 1, k, r2);

		double v1 = r1.get_val(fun_name);
		double v2 = r2.get_val(fun_name);

		df_deta = 0.5*(v2 - v1);

		// Guard a point of tripple zones intersection
		// assert(fabs(D.x_deta) + fabs(D.y_deta) + fabs(D.z_deta) > 1e-10);
	}
	//
	// d?_dzet
	//
	if (noBakK)
	{
		get_rec(zne, i, j, k + 0, r0);
		get_rec(zne, i, j, k + 1, r1);
		get_rec(zne, i, j, k + 2, r2);

		double v0 = r0.get_val(fun_name);
		double v1 = r1.get_val(fun_name);
		double v2 = r2.get_val(fun_name);

		df_ddze = 0.5*(-3.*v0 + 4.*v1 - v2);
	}
	else if (noFwdK)
	{
		get_rec(zne, i, j, k - 0, r0);
		get_rec(zne, i, j, k - 1, r1);
		get_rec(zne, i, j, k - 2, r2);

		double v0 = r0.get_val(fun_name);
		double v1 = r1.get_val(fun_name);
		double v2 = r2.get_val(fun_name);

		df_ddze = 0.5*(3.*v0 - 4.*v1 + v2);
	}
	else
	{
		get_rec(zne, i, j, k - 1, r1);
		get_rec(zne, i, j, k + 1, r2);

		double v1 = r1.get_val(fun_name);
		double v2 = r2.get_val(fun_name);

		df_ddze = 0.5*(v2 - v1);

		// Guard a point of tripple zones intersection
		//assert(fabs(D.x_dzet) + fabs(D.y_dzet) + fabs(D.z_dzet) > 1e-10);
	}
	//===================~calc deriv
	if (isnan(df_dksi) || isnan(df_deta) || isnan(df_ddze)) {
		wxLogError(_T("Error:TDomain::_calc_scalar_ked_deriv:Izone=%d, node=(%d, %d, %d)"),
			znode.iZone, i, j, k);
	}
	df_dked.set(df_dksi, df_deta, df_ddze);

}

void TDomain::calc_rec_grad(const t_ZoneNode& znode, mf::t_RecGrad& rec_grad) const {

	Tmtr3D mtr;

	t_Vec3Dbl v;

	_calc_scalar_ked_deriv(znode, 'x', v);
	
	mtr.x_dksi = v[0];
	mtr.x_deta = v[1];
	mtr.x_dzet = v[2];

	_calc_scalar_ked_deriv(znode, 'y', v);

	mtr.y_dksi = v[0];
	mtr.y_deta = v[1];
	mtr.y_dzet = v[2];

	_calc_scalar_ked_deriv(znode, 'z', v);

	mtr.z_dksi = v[0];
	mtr.z_deta = v[1];
	mtr.z_dzet = v[2];

	mtr.jac = jacobian(mtr);

	Tmtr3D::inv mtr_inv = mtr.getInverseMetric();

	_calc_scalar_ked_deriv(znode, 'u', v);
	mtr_inv.calc_xyz_deriv(v, rec_grad.ug);

	_calc_scalar_ked_deriv(znode, 'v', v);
	mtr_inv.calc_xyz_deriv(v, rec_grad.vg);

	_calc_scalar_ked_deriv(znode, 'w', v);
	mtr_inv.calc_xyz_deriv(v, rec_grad.wg);

	_calc_scalar_ked_deriv(znode, 'p', v);
	mtr_inv.calc_xyz_deriv(v, rec_grad.pg);

	_calc_scalar_ked_deriv(znode, 't', v);
	mtr_inv.calc_xyz_deriv(v, rec_grad.tg);

}

void TDomain::extract_profile_data(const mf::t_GeomPoint &xyz, const mf::t_ProfDataCfg& init_cfg, 
	std::vector<mf::t_Rec> &data, std::vector<mf::t_RecGrad> &data_derivs) const{

	wxLogMessage(_T("Extracting profile data"));

	t_ZoneNode surf_znode;

	_extract_profile_data_grdline(xyz, surf_znode);

	switch (_profile_cfg.BLThickCalcType)
	{
	case t_BLThickCalcType::BLTHICK_BY_VDERIV:
	case t_BLThickCalcType::BLTHICK_BY_ENTHALPY:

		_extract_profile_data_blbound(surf_znode, init_cfg, data, data_derivs);
		break;

	case t_BLThickCalcType::BLTHICK_FULL_GRIDLINE:

		_get_profile_data_grdline(data, data_derivs);
		break;

	default:
		ssuGENTHROW(_T("MF.CGNS error: profile exctraction not implemented for BLCalcType option specified"));
	};

}

// IMPORTANT TODO: this only works for ortho grids (see eta calculations)
void TDomain::_extract_profile_data_blbound(const t_ZoneNode& surf_znode, const mf::t_ProfDataCfg& init_cfg, 
	std::vector<mf::t_Rec>& data, std::vector<mf::t_RecGrad>& data_derivs) const{

		t_ProfScales bl_scales;
		t_ZoneNode outer_znode;

		std::vector<t_ZoneNode> raw_profile;

		_calc_bl_thick(surf_znode, bl_scales, raw_profile);


		//const double total_thick = 
		//	(init_cfg.ThickFixed >0.0) ? init_cfg.ThickFixed : bl_scales.thick_scale * init_cfg.ThickCoef;

		double total_thick = bl_scales.thick_scale * init_cfg.ThickCoef;
		if (init_cfg.ThickFixed > 0.0) total_thick = init_cfg.ThickFixed;

		// debug
		wxLogMessage(_T("Thick fixed = %lf"), init_cfg.ThickFixed);
		wxLogMessage(_T("Total thick = %lf"), total_thick);

		t_GeomPoint cur_xyz, surf_xyz;
		t_Rec cur_rec;
		t_Vec3Dbl rvec;

		get_rec(raw_profile[0], cur_rec);
		surf_xyz.set(cur_rec);

		int total_nodes = 0;
		double eta = 0.0;

		bool thick_ok = false;

		for (int m = 0; m < raw_profile.size(); m++) {

			get_rec(raw_profile[m], cur_rec);
			cur_xyz.set(cur_rec);
			matrix::base::minus<double, double>(cur_xyz, surf_xyz, rvec);
			eta = rvec.norm();
			total_nodes = m + 1;

			if (eta > total_thick) {

				thick_ok = true;

				break;
			}

		}

		data.resize(total_nodes);
		data_derivs.resize(total_nodes);

		if (!thick_ok) {
			wxLogError(_T("Error in _extract_profile_data_blbound: Requested thickness \
						  is greater than extracted from cgns domain"));

			for (int p = 0; p < total_nodes; p++) {

				get_rec(raw_profile[p], data[p]);
				calc_rec_grad(raw_profile[p], data_derivs[p]);
			} 
		}
		else {

			for (int p = 0; p < total_nodes - 1; p++) {

				get_rec(raw_profile[p], data[p]);
				calc_rec_grad(raw_profile[p], data_derivs[p]);

			} 

			// do some interpolation for the last point
			// we want thickness to be exactly total_thick
			t_GeomPoint xyz1, xyz2;
			mf::t_Rec r1, r2;
			mf::t_RecGrad r1_d, r2_d;

			get_rec(raw_profile[total_nodes - 2], r1);
			xyz1.set(r1);
			calc_rec_grad(raw_profile[total_nodes - 2], r1_d);

			matrix::base::minus<double, double>(xyz1, surf_xyz, rvec);

			double d1 = rvec.norm();

			get_rec(raw_profile[total_nodes - 1], r2);
			xyz2.set(r2);
			calc_rec_grad(raw_profile[total_nodes - 1], r2_d);

			matrix::base::minus<double, double>(xyz2, surf_xyz, rvec);

			double d2 = rvec.norm();

			double coef = (total_thick - d1) / (d2 - d1);

			if (coef<0.0 || coef>1.0)
				wxLogMessage(_T("Error in _extract_profile_data_blbound: interpolation coef is %lf, should be between 0 and 1"), coef);

			mf::t_Rec rec_intp = mf::t_Rec::lin_comb(1.0 - coef, r1, coef, r2);
			mf::t_RecGrad rec_derivs_intp = mf::t_RecGrad::lin_comb(1.0 - coef, r1_d, coef, r2_d);

			data[total_nodes - 1] = rec_intp;
			data_derivs[total_nodes - 1] = rec_derivs_intp;

		}

}

// Extract full domain gridline as a profile
void TDomain::_get_profile_data_grdline(std::vector<mf::t_Rec>& data, std::vector<mf::t_RecGrad> & data_derivs) const{

	const int nnodes = _pGrdLine->znodes.size();

	data.resize(nnodes);
	data_derivs.resize(nnodes);

	for (int p = 0; p < nnodes; p++) {

		get_rec(_pGrdLine->znodes[p], data[p]);
		calc_rec_grad(_pGrdLine->znodes[p], data_derivs[p]);
	}

}

void TDomain::_extract_profile_data_grdline(const t_ZoneNode& surf_znode) const{

	if (surf_znode != _pGrdLine->surf_znode) {

		wxLogMessage(_T("Extract full gridline"));

		t_DomainGrdLine grd_line(*this);

		grd_line.init(surf_znode);

		_pGrdLine->set(surf_znode, grd_line);
	}
	else {
		wxLogMessage(_T("Gridline already extracted, using cashed data"));
	}

}

// Extract full domain gridline as a set of indexes
void TDomain::_extract_profile_data_grdline(const mf::t_GeomPoint& xyz, t_ZoneNode& surf_znode) const{

	surf_znode = _get_nrst_node_surf(xyz);

	_extract_profile_data_grdline(surf_znode);
}

void TDomain::get_wall_gridline(const mf::t_GeomPoint& xyz) {

	t_DomainGrdLine grdLine(*this);

	mf::t_Rec rec;

	t_ZoneNode znode = _get_nrst_node_surf(xyz);

	wxLogMessage(_T("GetWallGridLine: starting from zone #%d, [i,j,k]=[%d, %d, %d]"),
		znode.iZone, znode.iNode.i, znode.iNode.j, znode.iNode.k);

	wxLogMessage(_T("GetWallGridLine: assuming wall grid line is along Xaxis (FacePos=Xmin)!!!"));
	wxLogMessage(_T("GetWallGridLine: marching from starting point in i+ direction ('downstream')!!!"));

	znode.iFacePos = get_cgns_params().FacePosStarting;

	grdLine.init(znode);

	wxLogMessage(_T("Writing wall gridline to output/ppoints_wall_grdline.dat"));

	grdLine.dump_as_ppoints("output/ppoints_wall_grdline.dat");

};

void TDomain::dump_rec_derivs(const mf::t_GeomPoint& xyz) const {

	t_ZoneNode znode = _get_nrst_node_raw(xyz);

	mf::t_Rec rec;

	get_rec(znode, rec);

	mf::t_GeomPoint xyz_node = rec.get_xyz();

	double dst = (xyz - xyz_node).norm();

	wxLogMessage(_T("distance to node:%lf"), dst);

	t_RecGrad rec_grad;
	calc_rec_grad(znode, rec_grad);

	wxLogMessage(rec_grad.to_wxstr());

};

void t_GrdLineFullCashed::set(const t_ZoneNode& a_surf_znode, const t_DomainGrdLine& grd_line) {

	surf_znode = a_surf_znode;

	const int nnodes_new = grd_line.znodes.size();

	if (nnodes_new != znodes.size()) znodes.resize(nnodes_new);

	for (int p = 0; p<nnodes_new; p++) znodes[p] = grd_line.znodes[p];
};
