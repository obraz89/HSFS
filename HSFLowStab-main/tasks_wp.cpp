#include "stdafx.h"

#include "tasks.h"

#include "solvers_glob.h"

#include "mpi.h"

#include "hdf5.h"

#include "msg_pack.h"

#include "WPTrackBase.h"

using namespace stab;
using namespace hsstab;

#define NMAX_FNAME_LEN 100
#define FILENAME_H5    "output/stab.h5"

void retrace_single_WP(const int wp_id, const int ppid, const bool a_do_simple_retrace, 
	const stab::t_WPRetraceMode a_mode_retrace, t_WPLine2H5Arr& a_arr);

void write_wpdata(hid_t file, hid_t group, char* dsname, const t_WPLine2H5Arr& arr);
void read_wpdata(hid_t file, char* ds_abs_name, t_WPLine2H5Arr& arr);

void task::retrace_MPI(stab::t_WPRetraceMode a_mode_retrace) {

	// rank of worker, total number of workers 
	int mpi_rank, mpi_size, len;

	char hostname[MPI_MAX_PROCESSOR_NAME];

	MPI_Get_processor_name(hostname, &len);

	const int MASTER = 0;

	// distribute wps through workers
	const task::TTaskParams& gtp = g_taskParams;

	const bool do_simple_retrace = (a_mode_retrace == stab::t_WPRetraceMode::WB_FIXED);

	const int npp = g_pStabDB->get_npoints();

	const int NWPGlob = do_simple_retrace ? gtp.N_b * gtp.N_w : gtp.N_b*gtp.N_w*npp;

	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	const int nAvg = NWPGlob / mpi_size;
	const int nResdl = NWPGlob % mpi_size;

	// interval of global pave points on current worker

	int wpid_s = mpi_rank * nAvg + ((mpi_rank<nResdl) ? mpi_rank : nResdl);
	int wpid_e = wpid_s + nAvg + ((mpi_rank<nResdl) ? 1 : 0) - 1;

	const int NWPLoc = wpid_e - wpid_s + 1;

	wxLogMessage(_("* MPI rank %d owns range of wp_ids : %d-%d"), mpi_rank, wpid_s, wpid_e);

	// table to map wpid (or "tag" in mpi msg) to mpi rank
	int* wpid_to_rank_map = new int[NWPGlob];

	for (int i = 0; i < mpi_size; i++) {

		int wpid_s = i * nAvg + ((i<nResdl) ? i : nResdl);
		int wpid_e = wpid_s + nAvg + ((i<nResdl) ? 1 : 0) - 1;

		for (int j = wpid_s; j <= wpid_e; j++) wpid_to_rank_map[j] = i;

	}

	int* wpid_to_ppid_map = new int[NWPGlob];
	
	if (!do_simple_retrace) {

		for (int k = 0; k<npp; k++)
			for (int j = 0; j<gtp.N_w; j++)
				for (int i = 0; i < gtp.N_b; i++) {
					int wpid = k*gtp.N_w*gtp.N_b + j*gtp.N_w + i;
					wpid_to_ppid_map[wpid] = k;
		}

	}
	else { for (int i = 0; i < NWPGlob; i++) wpid_to_ppid_map[i] = -1; }


	// retrace & pack

	t_WPLine2H5Arr* arr_pack = new t_WPLine2H5Arr[NWPLoc];

	for (int wpid = wpid_s; wpid <= wpid_e; wpid++) {

		wxLogMessage(_T("rank=%d, starting retrace wpid=%d, wpid_s=%d, wpid_e=%d"), mpi_rank, wpid, wpid_s, wpid_e);

		retrace_single_WP(wpid, wpid_to_ppid_map[wpid], do_simple_retrace, a_mode_retrace, arr_pack[wpid-wpid_s]);

	}

	wxLogMessage(_T("rank=%d finished retrace, waiting other workers"), mpi_rank);

	MPI_Barrier(MPI_COMM_WORLD);

	// send messages

	if (mpi_rank!=MASTER){

		double* mpi_buff = new double[NMAX_WPBUFF_DBL];

		for (int wpid = wpid_s; wpid <= wpid_e; wpid++) {

			arr_pack[wpid - wpid_s].pack_to_mpi_msg(mpi_buff);

			MPI_Send(mpi_buff, NMAX_WPBUFF_DBL, MPI_DOUBLE, 0, wpid, MPI_COMM_WORLD);			

		}
	}

	// do io
	if (mpi_rank==MASTER) {

		hid_t        file, group;
		herr_t status;

		file = H5Fcreate(FILENAME_H5, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

		group = H5Gcreate(file, "/WPData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		const hsize_t dims_attr = 1;
		hid_t ds_attr_id = H5Screate_simple(1, &dims_attr, NULL);

		/* Create a group attribute - number of wplines. */
		hid_t attr_id = H5Acreate2(group, "nwplines", H5T_NATIVE_INT, ds_attr_id,
			H5P_DEFAULT, H5P_DEFAULT);

		status = H5Awrite(attr_id, H5T_NATIVE_INT, &NWPGlob);

		status = H5Aclose(attr_id);

		// write local data

		for (int wpid = 0; wpid < NWPLoc; wpid++) {

			char dsname[32];

			sprintf(dsname, "%s%d", "Data", wpid);

			write_wpdata(file, group, dsname, arr_pack[wpid]);
		}

		double* mpi_buff = new double[NMAX_WPBUFF_DBL];

		t_WPLine2H5Arr arr;

		// receive messages for the rest of wpline ids
		for (int wpid = NWPLoc; wpid < NWPGlob; wpid++) {

			int rcv_from_rank = wpid_to_rank_map[wpid];

			MPI_Status mpi_status;

			MPI_Recv(mpi_buff, NMAX_WPBUFF_DBL, MPI_DOUBLE, rcv_from_rank, wpid, MPI_COMM_WORLD, &mpi_status);

			arr.unpack_from_mpi_msg(mpi_buff);

			char dsname[32];

			sprintf(dsname, "%s%d", "Data", wpid);

			write_wpdata(file, group, dsname, arr);
		}


		status = H5Gclose(group);
		status = H5Fclose(file);

	}

}

bool search_init_wave(const int wpid, const int a_ppoint_id, const bool a_do_simple_retrace, const bool is_update_dels,
	const stab::t_WPRetraceMode a_mode_retrace, t_WCharsLoc& a_wchars, mf::t_GeomPoint& a_start_xyz) {

	bool gs_success = false;

	// npp - number of pave point
	const int npp = g_pStabDB->get_npoints();

	const task::TTaskParams& gtp = g_taskParams;

	const int Nw = gtp.N_w;
	const int Nb = gtp.N_b;

	int i, j;

	if (a_do_simple_retrace) {

		i = wpid % Nw;
		j = wpid / Nw;

	}
	else {

		j = wpid / (npp*Nw);
		i = (wpid - j*npp*Nw) / npp;

	};

	double db_dim = (gtp.N_b > 1) ? (gtp.b_dim_max - gtp.b_dim_min) / double(gtp.N_b - 1) : 0.0;
	double dw_dim = (gtp.N_w > 1) ? (gtp.w_dim_max - gtp.w_dim_min) / double(gtp.N_w - 1) : 0.0;

	double b_ldim = gtp.b_dim_min + j*db_dim;
	double w_ldim = gtp.w_dim_min + i*dw_dim;

	// if doing simplified retrace search for instab at all pave points until instab is found
	// if doing exact retrace, do gs at a given point

	const int npp_s = a_do_simple_retrace ? 0 : a_ppoint_id;
	const int npp_e = a_do_simple_retrace ? npp - 1 : a_ppoint_id;

	for (int npp = npp_s; npp <= npp_e; npp++) {

		try {

			const mf::t_GeomPoint& test_xyz = g_pStabDB->get_pave_pt(npp).xyz;

			g_pGSSolverSpat->setContext(test_xyz);

			// TODO: write fixed dels, multiproc variant
			if (is_update_dels) {

				g_pMFDomain->set_stored_dels(g_pGSSolverSpat->get_stab_scales().Dels);
			}

			g_pStabSolver->setContext(test_xyz);

			const t_StabScales& sc = g_pStabSolver->get_stab_scales();

			t_WCharsLocDim wch_ldim;
			wch_ldim.set_treat(stab::t_TaskTreat::SPAT);
			wch_ldim.set_scales(sc);

			wch_ldim.a = 0.0;
			wch_ldim.b = b_ldim;
			wch_ldim.w = w_ldim;

			t_WCharsLoc wchars = wch_ldim.to_nondim(sc);

			gs_success = false;
			stab::t_LSCond cond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED);
			wxLogMessage(_T("start gs for WPLine = [i_w=%d, j_b=%d]\n\t pave_point=%d"), i, j, a_ppoint_id);
			{
				std::vector<t_WCharsLoc> init_waves_raw;
				std::vector<t_WCharsLoc> init_waves_filtered;

				init_waves_raw = g_pGSSolverSpat->getInstabModes(wchars);

				init_waves_filtered = g_pStabSolver->filter_gs_waves_spat(init_waves_raw, cond);

				if (init_waves_filtered.size()>0) {

					t_WCharsLoc wch_cand = t_WCharsLoc::find_max_instab_spat(init_waves_filtered);
					// TODO: option to start retrace from stable regions
					// for now, start only from unstable
					if (wch_cand.a.imag() < 0.0) {
						gs_success = true;
						a_wchars = t_WCharsLoc::find_max_instab_spat(init_waves_filtered);
						a_start_xyz = test_xyz;
						return true;
					}

				}

			}



		}
		catch (const t_GenException& ex) {
			wxLogMessage(ex.what());
			return false;
		}
		catch (...) {

			wxLogMessage(_T("RetraceMPI Error, WPLine = [i_w=%d, j_b=%d]\n\t pave_point=%d"), i, j, npp);
			return false;
		}

	}	// ~pave points loop

	return false;

}


void retrace_single_WP(const int wpid, const int ppid, const bool a_do_simple_retrace,
	stab::t_WPRetraceMode a_mode_retrace, t_WPLine2H5Arr& a_arr) {

	char szFname[NMAX_FNAME_LEN];
	char fout_maxnfactor_str[NMAX_FNAME_LEN];
	char fout_wplines_str[NMAX_FNAME_LEN];
	char fout_wplines_disp_str[NMAX_FNAME_LEN];
	char fout_wpline_disp_dump[NMAX_FNAME_LEN];

	sprintf(fout_wplines_str, "%s/wplines_mode%d.dat",
		hsstab::OUTPUT_DIR.ToAscii(), a_mode_retrace);

	sprintf(fout_maxnfactor_str, "%s/max_N_mode%d.dat",
		hsstab::OUTPUT_DIR.ToAscii(), a_mode_retrace);

	TCapsWPTrack& caps_wp = G_Plugins.get_caps_wp();

	stab::t_WPTrackBase* wp_line = caps_wp.create_wp_track(*g_pMFDomain);

	wp_line->init(G_Plugins.get_plugin(plgWPTrack));

	const task::TTaskParams& gtp = g_taskParams;

	{

	
		t_WCharsLoc wchars;
		mf::t_GeomPoint start_xyz;
			try {

				bool is_update_dels = wp_line->is_update_dels_at_rstart();

				bool gs_ok = search_init_wave(wpid, ppid, a_do_simple_retrace, is_update_dels, a_mode_retrace, wchars, start_xyz);

				// Init wave found, retrace WP
				if (gs_ok) {

					wxLogMessage(_T("\tstart retrace for WPLine"));

					std::wcout << _T("Init for wpline:") << wchars;

					wp_line->retrace(start_xyz, wchars, *g_pStabSolver, *g_pGSSolverSpat, a_mode_retrace);

					wp_line->print_to_file(fout_wplines_str, std::ios::app);

					wp_line->pack_to_arr(a_arr);

					//g_pStabDB->update(*wp_line);

					// no need to retrace same WP from different points
					// correct only for 2D configurations

				}
				else {
					wxLogMessage(_T("\tgs failed for WPLine"));
				}

			}
			catch (const t_GenException& ex) {
				wxLogMessage(ex.what());
			}
			catch (...) {

				wxLogMessage(_T("RetraceMPI Error, WPLine=%d]\n\t pave_point=%d"), wpid, ppid);
			}

		
	}	// ~WP Lines loop

	delete wp_line;
	return;

}

void task::postproc_retrace() {

	hid_t file = H5Fopen(FILENAME_H5, H5F_ACC_RDONLY, H5P_DEFAULT);
	hid_t group = H5Gopen(file, "/WPData", H5P_DEFAULT);
	herr_t status;

	t_WPLine2H5Arr arr;

	hid_t attr_nwp = H5Aopen_name(group,"nwplines");
	int nwp=0;

	H5Aread(attr_nwp, H5T_NATIVE_INT, &nwp);

	wxLogMessage(_T("Number of WpLines to read:nwplines=%d"), nwp);

	int npnts = g_pStabDB->get_npoints();

	t_WPLineRec* env_data = new t_WPLineRec[npnts];

	int* env_wpids = new int[npnts];

	for (int i = 0; i < nwp; i++) {

		char wp_dset_name[32];

		sprintf(wp_dset_name, "/WPData/Data%d", i);

		wxLogMessage(_T("\t Parsing %s"), wxStrdup(wxConvertMB2WX(wp_dset_name)));

		read_wpdata(file, wp_dset_name, arr);

		t_WPLineRec cur_env_rec;

		for (int n = 0; n < npnts; n++) {

			const stab::t_PavePoint& pnt = g_pStabDB->get_pave_pt(n);

			t_WPLineRec& env_data_rec = env_data[n];

			arr.interpolate_to_point(pnt.xyz, cur_env_rec);

			if (cur_env_rec.n_factor > env_data_rec.n_factor) {

				env_data_rec = cur_env_rec;

				env_wpids[n] = i;
			}

		}
	}

	// output envelope data
	{
		std::ofstream ofstr_env("output/N_fact_envelope.dat");

		ofstr_env <<"x\ty\tz\tN_fact\talpha_dim\tkn_dim\tbeta_dim\tfreq_dim\t";
		ofstr_env <<"dN_dw_gndim\tdN_db_gndim\td2N_dw2_gndim\td2N_dwb_gndim\td2N_db2_gndim\n";

		for (int i = 0; i < npnts; i++) {
			const t_WPLineRec& env_rec = env_data[i];
			const mf::t_GeomPoint& xyz = g_pStabDB->get_pave_pt(i).xyz;
			const t_WCharsGlobDim& wc = env_rec.wave_chars.to_dim();
			ofstr_env << xyz.x() << "\t" << xyz.y() << "\t" << xyz.z() << "\t"
				<< env_rec.n_factor << "\t"
				<< wc.a.real() << "\t" << wc.kn.real() << "\t" << wc.b.real() << "\t"
				<< wc.w.real() << "\t"
				<< env_rec.dN_dw_gndim << "\t" << env_rec.dN_db_gndim << "\t"
				<< env_rec.d2N_dw2_gndim << "\t" << env_rec.d2N_dwb_gndim << "\t" << env_rec.d2N_db2_gndim << "\n";
		}
	}

	// output dispersion data for wplines that produce envelope

	TCapsWPTrack& caps_wp = G_Plugins.get_caps_wp();

	stab::t_WPTrackBase* wp_line = caps_wp.create_wp_track(*g_pMFDomain);

	for (int i = 0; i < npnts; i++) {

		char wp_dset_name[32];

		sprintf(wp_dset_name, "/WPData/Data%d", env_wpids[i]);

		wxLogMessage(_T("\t Reading envelope wpline : %s"), wxStrdup(wxConvertMB2WX(wp_dset_name)));

		read_wpdata(file, wp_dset_name, arr);

		wp_line->unpack_from_arr(arr);

		char wp_disp_fname[32];

		sprintf(wp_disp_fname, "output/disp_data_env_%d.dat", env_wpids[i]);

		wp_line->print_dispersion_data_full(wp_disp_fname);

	}

	status = H5Aclose(attr_nwp);
	status = H5Fclose(file);

	delete[] env_data;

}
// rank of wpline data array in hdf5 file
#define RWP 2

void write_wpdata(hid_t file, hid_t group,char* dsname, const t_WPLine2H5Arr& arr){

	hsize_t dims[RWP] = { arr.nrecs , N_WPREC_H5_LEN };
	hid_t dataspace = H5Screate_simple(RWP, dims, NULL);

	hid_t dataset = H5Dcreate2(group, dsname, H5T_NATIVE_DOUBLE, dataspace,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	double* data = arr.cont;

	herr_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		H5P_DEFAULT, data);

	// attributes

	const hsize_t dims_attr = 1;
	hid_t ds_attr_id = H5Screate_simple(1, &dims_attr, NULL);

	/* Create a dataset attribute. */
	hid_t attr_id = H5Acreate2(dataset, "nrecs", H5T_STD_I32BE, ds_attr_id,
		H5P_DEFAULT, H5P_DEFAULT);

	/* Write the attribute data. */
	int attr_data = arr.nrecs;
	status = H5Awrite(attr_id, H5T_NATIVE_INT, &attr_data);

	/* Close the attribute. */
	status = H5Aclose(attr_id);

	/* Close the dataspace. */
	status = H5Sclose(ds_attr_id);


	status = H5Dclose(dataset);
	status = H5Sclose(dataspace);

	return;
}

void read_wpdata(hid_t file, char* ds_abs_name, t_WPLine2H5Arr& arr){

	herr_t       status;

	hid_t dataset = H5Dopen(file, ds_abs_name, H5P_DEFAULT);

	hid_t attr = H5Aopen_by_name(file, ds_abs_name,
		"nrecs", H5P_DEFAULT, H5P_DEFAULT);

	status = H5Aread(attr, H5T_NATIVE_INT, &arr.nrecs);

	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5P_DEFAULT, H5P_DEFAULT,
		H5P_DEFAULT, arr.cont);
		
}


