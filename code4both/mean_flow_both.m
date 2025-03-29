% mean flow for both

%% CW 和 CCW各自预处理 并解码

%% 计算最好的切分相位
cd E:\code\theta_precession_gamma\code4both
% 先遍历
theta_sequence_averaging_traverse_both
% 再通过计算指标来计算最好的相位
theta_sequence_averaging_goodphase

%% 把goodphase写在directories里

%% 计算三个条件的sequence
doall_Avg_theta_seq_both
doall_Avg_theta_seq_both_exfg
doall_Avg_theta_seq_both_dsfg
%% 计算三个条件的sequence，但是分成三个部分
doall_Avg_theta_seq_both_section
doall_Avg_theta_seq_both_exfg_section
doall_Avg_theta_seq_both_dsfg_section

%% 计算指标
doall_plot_ProbDiff_both
doall_plot_WeightCorr_both
doall_plot_WeightCorr_each_both
doall_plot_FlipCorr_both

doall_plot_ProbDiff_both_section
doall_plot_WeightCorr_both_section
doall_plot_WeightCorr_each_both_section  
doall_plot_FlipCorr_both_section



