#include "stdafx.h"

#include "tasks.h"

#include "solvers_glob.h"

#include "mpi.h"

#include "hdf5.h"

#include "msg_pack.h"

#include "WPTrackBase.h"

using namespace stab;
using namespace hsstab;

static const int MAX_FNAME_LEN = 100;

/************************************************************************/
/* serializer for a gs record
** used to send info from worker to master
** then it is unpacked and written to stability db
*/
/************************************************************************/

/*class t_MsgWPLine : public t_MsgBase {

public:

	int SerSize;

	int getSerSize() const { return SerSize; };

	t_MsgWPLine(const t_WPTrackBase& wp_line);
	void write2file(std::wofstream& ofstr) const;

};

t_MsgGSRec::t_MsgGSRec(double *cont) :t_MsgBase(cont) {}

t_MsgGSRec::t_MsgGSRec(const t_MsgBase& bb) : t_MsgBase(bb) {}

const int t_MsgGSRec::SerSize = 17;

void t_MsgGSRec::write2file(std::wofstream& ofstr) const {

	int pid = _pCont[0];
	double x = _pCont[1];
	double y = _pCont[2];
	double z = _pCont[3];
	int ok = _pCont[4];

	double w_ar = _pCont[5];
	double w_ai = _pCont[6];
	double w_br = _pCont[7];
	double w_bi = _pCont[8];
	double w_wr = _pCont[9];
	double w_wi = _pCont[10];

	double d_ar = _pCont[11];
	double d_ai = _pCont[12];
	double d_br = _pCont[13];
	double d_bi = _pCont[14];
	double d_wr = _pCont[15];
	double d_wi = _pCont[16];


	ofstr << pid << "\t" << x << "\t" << y << "\t" << z << "\t" << ok << "\t"
		<< w_ar << "\t" << w_ai << "\t" << w_br << "\t" << w_bi << "\t" << w_wr << "\t" << w_wi << "\t"
		<< d_ar << "\t" << d_ai << "\t" << d_br << "\t" << d_bi << "\t" << d_wr << "\t" << d_wi
		<< "\n";


};*/

void task::retrace_MPI(stab::t_WPRetraceMode a_mode_retrace) {

	char szFname[MAX_FNAME_LEN];
	char fout_maxnfactor_str[MAX_FNAME_LEN];
	char fout_wplines_str[MAX_FNAME_LEN];
	char fout_wplines_disp_str[MAX_FNAME_LEN];
	char fout_wpline_disp_dump[MAX_FNAME_LEN];

	sprintf(fout_wplines_str, "%s/wplines_mode%d.dat",
		hsstab::OUTPUT_DIR.ToAscii(), a_mode_retrace);

	sprintf(fout_maxnfactor_str, "%s/max_N_mode%d.dat",
		hsstab::OUTPUT_DIR.ToAscii(), a_mode_retrace);

	TCapsWPTrack& caps_wp = G_Plugins.get_caps_wp();

	stab::t_WPTrackBase* wp_line = caps_wp.create_wp_track(*g_pMFDomain);

	wp_line->init(G_Plugins.get_plugin(plgWPTrack));

	const task::TTaskParams& gtp = g_taskParams;

	// single-proc variant
	// nwp - number of wape pack line
	const int nwp_s = 0;
	const int nwp_e = gtp.N_b * gtp.N_w - 1;

	// npp - number of pave point
	// move from back to leading edge (assuming x increasing in pave points)
	const int npp_s = g_pStabDB->get_npoints() - 1;
	const int npp_e = 0;

	for (int i=0; i<gtp.N_w; i++)
		for (int j=0; j<gtp.N_b; j++){

		 int nwp = j*gtp.N_w + i;

		 double db_dim = (gtp.N_b > 1) ? (gtp.b_dim_max - gtp.b_dim_min) / double(gtp.N_b - 1) : 0.0;
		 double dw_dim = (gtp.N_w > 1) ? (gtp.w_dim_max - gtp.w_dim_min) / double(gtp.N_w - 1) : 0.0;

		 double b_ldim = gtp.b_dim_min + j*db_dim;
 		 double w_ldim = gtp.w_dim_min + j*dw_dim;


		 for (int npp = npp_s; npp >= npp_e; npp--) {

			try {

				const mf::t_GeomPoint& test_xyz = g_pStabDB->get_pave_pt(npp).xyz;

				g_pGSSolverSpat->setContext(test_xyz);

				g_pStabSolver->setContext(test_xyz);

				const t_StabScales& sc = g_pStabSolver->get_stab_scales();

				t_WCharsLocDim wch_ldim;
				wch_ldim.set_treat(stab::t_TaskTreat::SPAT);
				wch_ldim.set_scales(sc);

				wch_ldim.a = 0.0;
				wch_ldim.b = b_ldim;
				wch_ldim.w = w_ldim;

				t_WCharsLoc wchars = wch_ldim.to_nondim(sc);

				// do global search of instab
				bool gs_success = false;
				stab::t_LSCond cond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED);
				wxLogMessage(_T("start gs for WPLine = [i_w=%d, j_b=%d]\n\t pave_point=%d"), i, j, npp);
				{
					std::vector<t_WCharsLoc> init_waves_raw;
					std::vector<t_WCharsLoc> init_waves_filtered;

					init_waves_raw = g_pGSSolverSpat->getInstabModes(wchars);

					init_waves_filtered = g_pStabSolver->filter_gs_waves_spat(init_waves_raw, cond);

					if (init_waves_filtered.size()>0){

						gs_success = true;
						wchars = t_WCharsLoc::find_max_instab_spat(init_waves_filtered);

					}

				}

				// Init wave found, retrace WP
				if (gs_success) {

					int perc_complete = double(nwp) / double(nwp_e)*100.;
					wxLogMessage(_T("\tstart retrace for WPLine"));

					std::wcout << _T("Init for wpline:") << wchars;

					wp_line->retrace(test_xyz, wchars, *g_pStabSolver, *g_pGSSolverSpat, a_mode_retrace);

					wp_line->print_to_file(fout_wplines_str, std::ios::app);

					g_pStabDB->update(*wp_line);

					break;

				}
				else {
					wxLogMessage(_T("\tgs failed for WPLine"));
				}

			}
			catch (const t_GenException& ex) {
				wxLogMessage(ex.what());
			}
			catch (...) {

				wxLogMessage(_T("RetraceMPI Error, WPLine = [i_w=%d, j_b=%d]\n\t pave_point=%d"), i, j, npp);
			}

		}	// ~pave points loop

		
	}	// ~WP Lines loop

	//StabDB.to_cone_ref_frame(HALF_CONE_ANGLE);

	g_pStabDB->write_to_file(fout_maxnfactor_str);

	delete wp_line;
	return;

}

void hdf5_test() {
	hid_t       file_id;

	herr_t      status;

	file_id = H5Fcreate("output/file.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Fclose(file_id);

	/*

	// Create the dataspace for the dataset.
	dims[0] = 4;
	dims[1] = 6;
	dataspace_id = H5Screate_simple(2, dims, NULL);


	// Create the dataset.
	dataset_id = H5Dcreate (file_id, "/dset", H5T_STD_I32BE,
	dataspace_id, H5P_DEFAULT, H5P_DEFAULT,
	H5P_DEFAULT);

	// Close the dataset and dataspace
	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);

	// read - write datasets

	status = H5Dwrite (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, dset_data);

	status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, dset_data);


	*/

}

