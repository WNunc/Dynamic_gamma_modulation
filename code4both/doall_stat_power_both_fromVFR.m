%根据 王雪玲 doall_VFR_gammaZ 脚本修改
% created by 王宁 on 2023年2月2日
%% static of slow and fast power from VFR
% 从VFR中找能量数据进行统计

clear
directories_allData_v0_allgood
%频段和速度上平均
resuletFolder = 'H:\neuralynx\gamma in sequence result';
file_input = 'Data_VFR_prerunning_all_trail_allgood.mat';% VFR result
file_output_mean = 'Data_VFR_prerunning_mean_allgood.mat';

% vr1=0.65;%wxl
vr1=5.05;
vr2=35.05;
% logvr = 0:0.2:6.2;
logvr = 3.35:0.15:6.2;
bin=0.5;
vr=2.^logvr*bin;
vr_center = mean([vr(1:end-1);vr(2:end)]);
[~,v1]=min(abs(vr_center-vr1));
[~,v2]=min(abs(vr_center-vr2));

isession = size(path,1);
ndirs = 1;
mean_VFR = cell(ndirs,1 );
mean_VFR_all = cell(ndirs, 5);%跑5次
theta_low = 6;
theta_high = 12;
gamma_slow_low = 25;
gamma_slow_high = 45;
gamma_fast_low = 55;
gamma_fast_high = 100;

sample_VFR_all = []; 
test_VFR_all = [];
pre_VFR_all = [];
nsn = 0; % 实际的session数
for ns = [1,2,4:6,8,10:13,15,16]%1:isession%
    nsn = nsn + 1;
    path_ns = path{ns};
    inFolder = [resuletFolder path_ns(13:end)];
    cd(inFolder);
    
    VFR_ns = strcat(inFolder,file_input);
    load(VFR_ns);
    trackdata_ns = trackdata{ns};
%     load(trackdata_ns,'n_prerunning','n_sample_test','sign_correct_test','sign_correct_sample');
%     load('Data_angle_ontrack.mat','n_test');
    speed_lim = v1:v2;
    
    theta_num = find(freq_gamma>=theta_low & freq_gamma<=theta_high);
    gamma_slow_num = find(freq_gamma>=gamma_slow_low & freq_gamma<=gamma_slow_high);
    gamma_fast_num = find(freq_gamma>=gamma_fast_low & freq_gamma<=gamma_fast_high);
    ind_fre={theta_num,gamma_slow_num,gamma_fast_num};
%     ntrail_session=[]; aa=[];
%     for nd= 1:ndirs
%         for nl=1:length(VFR_exp{1,1})
%             aa(nd,nl) = ~isempty(VFR_exp{1,1}{nd,nl});
%         end
%         ntrail_session(nd) = sum(aa(nd,:));
%     end
    
    for nd = 1:ndirs
        if nd == 1
            for nl = 1:length(VFR_trial)%1:ntrail_session(nd)
                for fre=1:length(ind_fre) %三个频段
                    pre_VFR {nsn,1}(nl,fre) = nanmean(nanmean(VFR_trial{nl}(ind_fre{fre},speed_lim)));
                end
            end
        else
            for nl = 1:ntrail_session(nd) %8个lap
                theta = NaN;
                gamma_slow = NaN;
                gamma_fast = NaN;
                if ~isempty(VFR_exp{1,1}{nd,nl})
                    if sign_correct_sample(nl) == 1&(sign_correct_test(nl)==0||sign_correct_test(nl)==1)
                        theta = nanmean(nanmean(VFR_exp{1,1}{nd,nl}(theta_num,:)));
                        gamma_slow = nanmean(nanmean(VFR_exp{1,1}{nd,nl}(gamma_slow_num,:)));
                        gamma_fast = nanmean(nanmean(VFR_exp{1,1}{nd,nl}(gamma_fast_num,:)));
                    end
                end
                if nd==2
                    sample_VFR {ns,1}(nl,1)= theta;
                    sample_VFR {ns,1}(nl,2)= gamma_slow;
                    sample_VFR {ns,1}(nl,3)= gamma_fast;
                else
                    test_VFR {ns,1}(nl,1)= theta;
                    test_VFR {ns,1}(nl,2)= gamma_slow;
                    test_VFR {ns,1}(nl,3)= gamma_fast;
                end
            end
        end
%         theta_ns = nanmean(nanmean(VFR_nd{1,nd}(theta_num,:)));%每天实验的平均
%         gamma_slow_ns = nanmean(nanmean(VFR_nd{1,nd}(gamma_slow_num,:)));
%         gamma_fast_ns = nanmean(nanmean(VFR_nd{1,nd}(gamma_fast_num,:)));
%         mean_VFR{nd,1}(1,ns) = theta_ns;
%         mean_VFR{nd,1}(2,ns) = gamma_slow_ns;
%         mean_VFR{nd,1}(3,ns) = gamma_fast_ns;
    end
%     sample_VFR_all = [sample_VFR_all;sample_VFR{ns,1}]; %把每天每一圈的连起来，做统计
%     test_VFR_all = [test_VFR_all;test_VFR{ns,1}];
    pre_VFR_all = [pre_VFR_all;pre_VFR{nsn,1}];
    
%     %===========================先把功率谱连起来，再算VFR=======================%
%     for fre=1:length(ind_fre)
%         pre_VFR_2 (ns,fre) = nanmean(nanmean(VFR_session_all{1}(ind_fre{fre},speed_lim))); %每个session内每一圈连起来的VFR
%         sample_VFR_2 (ns,fre) = nanmean(nanmean(VFR_session_all{2}(ind_fre{fre},speed_lim))); %每个session内每一圈连起来的VFR
%         test_VFR_2 (ns,fre) = nanmean(nanmean(VFR_session_all{3}(ind_fre{fre},speed_lim))); %每个session内每一圈连起来的VFR
%     end
end

thpower = reshape(pre_VFR_all(:,1),[5,nsn]);
sgpower = reshape(pre_VFR_all(:,2),[5,nsn]);
fgpower = reshape(pre_VFR_all(:,3),[5,nsn]);

plot_sem(thpower','lap num','zscored theta power')
ylim([-0.2 0.3]);axis square
plot_sem(sgpower','lap num','zscored sgamma power')
ylim([-0.2 0.3]);axis square
plot_sem(fgpower','lap num','zscored fgamma power')
ylim([-0.2 0.3]);axis square





% %把AD和con的数据分别存起来
% if icon == 1
%     mean_VFR_Con = mean_VFR;
%     sample_all_Con = sample_VFR_all;%每一天中的每一圈，不按天分开
%     test_all_Con = test_VFR_all;
%     pre_all_Con = pre_VFR_all;
%     sample_Con = sample_VFR;% 每一天中的每一圈，按天分开
%     test_Con = test_VFR;
%     pre_Con = pre_VFR;
%     pre_VFR_2_Con=pre_VFR_2;
%     sample_VFR_2_Con=sample_VFR_2;
%     test_VFR_2_Con=test_VFR_2;
%     cd(file_outfd0);
%     save(file_output_mean{icon},'sample_all_Con','test_all_Con','pre_all_Con','sample_Con','test_Con','pre_Con',...
%         'mean_VFR_Con','pre_VFR_2_Con','sample_VFR_2_Con','test_VFR_2_Con');
% else
%     mean_VFR_AD = mean_VFR;
%     sample_all_AD = sample_VFR_all;
%     test_all_AD = test_VFR_all;
%     pre_all_AD = pre_VFR_all;
%     sample_AD = sample_VFR;
%     test_AD = test_VFR;
%     pre_AD = pre_VFR;
%     pre_VFR_2_AD=pre_VFR_2;
%     sample_VFR_2_AD=sample_VFR_2;
%     test_VFR_2_AD=test_VFR_2;
%     cd(file_outfd0);
%     save(file_output_mean{icon},'sample_all_AD','test_all_AD','pre_all_AD','sample_AD','test_AD','pre_AD',...
%         'mean_VFR_AD','pre_VFR_2_AD','sample_VFR_2_AD','test_VFR_2_AD');
% end
