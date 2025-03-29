% main flow
%% 预处理
% 获取神经元和视频信息，计算ratemap
doall_data_pos2ang
doall_ratemap_singlelap_v2_wn
doall_ratemap_Allsegment_v1_WN
doall_ratemap_ALLlaps_v2_wn

%% 解码
% 全部session的ratemap解码
doall_decoding_CT_wn_v1
% 单圈的ratemap解码
doall_decoding_CT_v5
% 获取节律信息
doall_get_gamma_v7
fullAllLapScores

%% 切分theta sequence
theta_phase_cutting_traverseCW
theta_sequence_averaging_traverseCW
%% 找一个最好的相位出来（半自动）
theta_sequence_averaging_traverseCW_alllapseq
%% 计算相锁
doall_phaselocking_wn_v2
%% 观察相锁神经元的情况，计算相锁程度（2）
doall_phaselocking_alllap_obsinglelaps
doall_phaselocking_firstlap_obsinglelaps
%% 用最好的相位来切分theta相位
doall_theta_seqAveraging_cw
%% 计算gamma功率
doall_get_power_wn
%% 去掉快gamma相锁神经元解码+去除相同数量的非相锁神经元
doall_decoding_CT_exclude_cellfg
doall_decoding_CT_downsamp_cellfg
%% 用相同是时间窗口观察去除神经元后的解码效果
theta_sequence_averaging_exfg_samewin
theta_sequence_averaging_dsfg_samewin
%% 统计gamma功率（1）
doall_stat_power_cw
%% 计算sequence指标（3）
doall_plot_ProbDiff
doall_plot_WeightCorr
doall_plot_FlipCorr


