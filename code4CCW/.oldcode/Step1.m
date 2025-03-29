directories_allData_v1
%% Prepare
doall_data_pos2ang
doall_add_thetawin_v2
doall_ratemap_singlelap_v2
doall_ratemap_AllLaps_v2
doall_draw_ratemap_single_v0
doall_draw_ratemap_All_v0
%% Decoding
doall_decoding_CT_v4
drawDecoding
%% Power+Phase+EEG
doall_get_gamma_v7
%% Theta cycle
doall_get_phase_v0
doall_get_peak
doall_get_info
doall_get_power
doall_draw_cycle_v0
%% Phase precession
PhasePrecession_v5 % get pos X theta phase
SGPhasePrecession_v5
FGPhasePrecession_v5
PPshow % draw theta phase precession
%
doall_slowgamma_cycle % get List of each slow gamma cycle£¨3£©
doall_fastgamma_cycle % get List of each fast gamma cycle£¨5£©
doall_peak_firingrate % get peak firing rate of each cell
doall_theta_slowgammaP1_5 % peak firing rate>5Hz
doall_theta_slowgammaP2 % draw 2d histogram
doall_theta_fastgammaP1_5
doall_theta_fastgammaP2
doall_get_R % get R to measure the circular-linear fit
doall_get_aR_lap
%% Phase locking
doall_phaselock_fg_v2 % get the fast gamma phase-locking value of spikes
doall_phaselock_sg_v2 % get the slow gamma phase-locking value of spikes
doall_draw_R_PL % draw R and PLV of slow and fast gamma
