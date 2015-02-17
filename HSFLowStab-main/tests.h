namespace test{

	void king_m35();
	// TODO: after tests trnasform into universal
	// calc_eN_time_envelope
	void king_m35_eN_time_envelope();
	// TODO: after tests transform into universal
	// calc_eN_spat_fixedB
	void king_m35_eN_spat_fixedB();
	// test spatial gs vs time gs
	void king_m35_gs_spat_vs_time();

	void king_m35_generate_pave_points();

	void transhyb_base_08();
	void transhyb_base_wartman();
	void profile_compar();
	void itam_hz();
	void selfsim_M45_second_mode();
	void selfsim_M45_spectrum();
	void selfsim_M3_first_mode();
	//void selfsim_M3_spectrum(); 

	void selfsim_M2_CF();

	//math testing
	void test_smat();
	void test_conj_grad_min2D();

	//e-N closures testing
	void test_fixed_beta_calc();

	// test robustness of wchars checker by phase speed
	void test_wchars_filter_c();

	// verify Lapack GS & test against Petsc
	void gs_lapack_vs_petsc();

};