%%
clear all
close all
plot_sign=1;
title_session={'pre','semple','test'};
title_trial={'Lap1','Lap2','Lap3','Lap4','Lap5'};
% file_out='H:\group data\VFR_image\trail_eight\VFR_all_before_reward_8trail';
% fig_out = 'H:\group data\VFR_image\trail_eight\';
% file_input_power = 'Data_Z_gamma_multitaper_0p5s_8trail.mat';
% file_output = 'VFR_before_reward_8trail';
%  pre_use=[2 3 4 5 7 8 9 10]; %1-5为cw，6-10为ccw，不用pre的第一个trail，即1和6
% N_trail=8; %仅使用sample-test的前八个trail

% file_out='H:\group data\VFR_image\all_trail\VFR_all_before_reward_all_trail_sam1_test0_1.mat';
% fig_out = 'H:\group data\VFR_image\all_trail\before_reward\';
% file_input_power = 'Data_Z_gamma_multitaper_all_trail_sam1_test0_1.mat'; % 仅保留sample正确，而且test有停下的trail
file_input_power = 'Data_Z_gamma_multitaper_all_trail_prerunning3.mat';% 2  = speedlimit0.68；3 = speedlimit0.68
file_output = 'Data_VFR_prerunning_all_trail_allgood.mat';
pre_use=[1 2 3 4 5 6 7 8 9 10]; % 圈数
D = [1 2 1 2 1 2 1 2 1 2];% 奇数为cw，偶数为ccw
trial_use = [1 1 2 2 3 3 4 4 5 5];% trial数
Nlaps = length(pre_use);
N_trail=20;
A_start=1;
file_input_speed = 'Data_angle_ontrack.mat';  % input the speed information

directories_allData_v0_allgood

resultFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
case1 = '-ontrack';
ncells = 3;nspkikes = 5;
case2 = sprintf('_%ucell_%uspk',ncells, nspkikes);
midmod = '_cycmid';

zscoremax = 5;
% use all tetrodes with cell firing
CSClist = CSClist_CA1;

Fs=2000;
% set for running speed bins
vr1=5.05;
vr2=35.05;
logvr = 3.35:0.15:6.2;
bin=0.5;
vr=2.^logvr*bin;
tickx=[0 1 2 3 4 5 6]; %横坐标
vr_center = mean([vr(1:end-1);vr(2:end)]);
[~,v1]=min(abs(vr-vr1));
[~,v2]=min(abs(vr-vr2));
%vr=[0,vr];
nbin=length(vr_center);
invh = 5;  % test the effect of invh in the code:test_20200319.m

g1 = 4;

ndirs = 1;
narea = size(CSClist,2);
ia = narea;
% VFR_exp = cell(narea,ndirs,isession);
fre_bin=196; %vfr矩阵的行数
v_bin=30;%vfr矩阵的列数

VFR0_ns={};%所有session，每一圈功率谱的平均功率谱
VFR={}; %每天，每个session中，每一圈的功率谱
VFR_session_all_day=cell(1,ndirs); %每个session的功率谱先连起来再求VFR，每一天的都存在这里
VFR_test_all_day=cell(1,2); %test的正确或者错误lap连起来再求VFR，每一天的都存这里，第一列正确，第二列错误
VFR_session_all_day_mean=cell(1,ndirs);%所有天的平均，（每个session的功率谱先连起来再求VFR）

for ns = A_start:isession
    VFR_exp={};%每个session中，每一圈的功率谱
    VFR_nd={};%每个session，每一圈功率谱的平均功率谱
    path_ns = path{ns};
    cd(path_ns);
    goodphase_ns = Seq_cutPhase{ns,1};%***
    subfolder1 = Directionfolder{1};
    subfolder2 = Phasecutfolder;
    outFolder = [resultFolder path_ns(13:end)];
    inFolder = outFolder;
    case3 = num2str(goodphase_ns);
    %     trackdata_ns = trackdata{ns};
    %     load(trackdata_ns,'n_prerunning','n_sample_test','ts_prerunning_CCW_end','ts_prerunning_CCW_start',...
    %         'ts_prerunning_CW_end','ts_prerunning_CW_start','ts_sample_end','ts_sample_start',...
    %         'ts_test_end','ts_test_start','ts_sample_reward','ts_test_reward','sign_correct_sample','sign_correct_test');
    
    if exist([inFolder,file_input_power],'file')>0
        
        load([inFolder,file_input_power],'T_start','Sz','Sz_vel','Vel_win','freq_gamma','win0','winstep0','param','TS_CW','TS_CCW','TS_running');
        [~,g] = min(abs(freq_gamma-g1));
        Sz_part={};
        Sz_vel_part={};
        Vel_win_part = {};
        T_part={};
        ind_alllap=cell(1,ndirs);
        ind_singlap=cell(1,ndirs);
        %         ind_test_cor = [];
        %         ind_test_err = [];
        
        for nd=1:ndirs % nd=1,prerunning; nd=2,semple; nd=3,test
            if nd==1
                for nl=1:Nlaps
                    nlap_use=pre_use(nl);
                    ntrial = trial_use(nl);
                    d = D(nl);
                    ind = find(T_start>=TS_running(nlap_use,1) & T_start<=TS_running(nlap_use,2));
                    %                         ind = find(T_start>=ts_start_stop_point{1,nd}(nlap_use,1) & T_start<=ts_start_stop_point{1,nd}(nlap_use,2));
                    ind_alllap{nd} = [ind_alllap{nd},ind];
                    ind_singlap{nd}{ntrial,d} = ind;
                    Sz_part{nd}{ntrial,d}= Sz{ia}(ind,:); %有幅度限制的功率谱,
                    Sz_vel_part{nd}{ntrial,d}= Sz_vel{ia}(ind,:); %有速度限制的功率谱,
                    Vel_win_part{nd}{ntrial,d} = Vel_win(ind,:);%对应的速度
                    T_part{nd}{ntrial,d}= T_start(ind); %对应的时间
                end
                %                 else
                %                     if nd==2
                %                         for nlap=1:n_sample_test
                %                             if sign_correct_sample(nlap)==1 && ~isnan(sign_correct_test(nlap))%仅仅保留sample正确，而且test有停下的trail，否则vfr为空
                %                                 ind = find(T_start>=ts_start_stop_point{1,nd}(nlap,1) & T_start<=ts_start_stop_point{1,nd}(nlap,2));
                %                             else
                %                                 ind = [];
                %                             end
                %                             ind_session{nd} = [ind_session{nd},ind];
                %                             Sz_part{ia,1}{nd}{nlap}= Sz{ia,1}(ind,:);
                %                             Sz_vel_part{ia,1}{nd}{nlap}= Sz_vel{ia,1}(ind,:); %有速度限制的功率谱,
                %                             Vel_win_part{nd}{nlap} = Vel_win(ind,:);%对应的速度
                %                             T_part{nd}{nlap}= T_start(ind);
                %                         end
                %                     else
                %                         for nlap=1:min([n_sample_test,length(ts_start_stop_point{1,nd})])
                %                             ind = [];
                %                             if sign_correct_sample(nlap)==1 && sign_correct_test(nlap)==1%仅仅保留sample正确，而且test有停下的trail，否则vfr为空
                %                                 ind = find(T_start>=ts_start_stop_point{1,nd}(nlap,1) & T_start<=ts_start_stop_point{1,nd}(nlap,2));
                %                                 ind_test_cor = [ind_test_cor,ind];
                %                             elseif sign_correct_sample(nlap)==1 && sign_correct_test(nlap)==0
                %                                 ind = find(T_start>=ts_start_stop_point{1,nd}(nlap,1) & T_start<=ts_start_stop_point{1,nd}(nlap,2));
                %                                 ind_test_err = [ind_test_err,ind];
                %                             end
                %                             ind_session{nd} = [ind_session{nd},ind];
                %                             Sz_part{ia,1}{nd}{nlap}= Sz{ia,1}(ind,:);
                %                             Sz_vel_part{ia,1}{nd}{nlap}= Sz_vel{ia,1}(ind,:); %有速度限制的功率谱,
                %                             Vel_win_part{nd}{nlap} = Vel_win(ind,:);%对应的速度
                %                             T_part{nd}{nlap}= T_start(ind);
                %                         end
                %                     end
            end
        end
        
        
        %==========pre，semple，test各圈平均的VFR===========%
        %         for nd=1:ndirs
        %             if nd==1
        %                 for nl=1:length(Vel_win_part{nd}) %1-5为cw，6-10为ccw
        %                     gammaZ_VFR
        %                 end
        %             else
        %                 for nl=1:min([n_sample_test,length(ts_start_stop_point{1,nd})])
        %                     gammaZ_VFR
        %                 end
        %             end
        %             VFR0_ns{ia}{1,nd}(:,:,ns) = nanmean(VFR{ia}{nd,ns}(:,:,:),3);
        %             VFR_nd{ia,nd} = nanmean(VFR{ia}{nd,ns}(:,:,:),3);
        %         end
        
        %====画出一天内pre，semple，test各圈平均的VFR，并保存===%
        %         if plot_sign
        %             cl = [-0.5,0.5];
        %
        %             ffa = figure('unit','centimeters','position',[3 5 13*3 10]);
        %             for nd=1:ndirs
        %                 subplot(1,3,nd);
        %                 imagesc(logvr(v1:v2), freq_gamma(g:end), VFR_nd{ia,nd}(g:end,v1:v2-1),cl);
        %                 colormap(jet(128)); % 图像默认是蓝色到黄色，加了这一句是蓝色到红色
        %                 colorbar;
        %                 axis xy; axis tight;
        %                 set(gca, 'FontSize',18);%改变字体大小
        %                 set(gca, 'XTick',tickx);
        %                 set(gca, 'XTickLabel',2.^tickx*bin);
        %                 xlabel('Running speed (cm/s)');
        %                 ylabel('Frequency (Hz)');
        %                 title([title_session{nd},'-lap-mean']);
        %                 set(gcf,'PaperPositionMode','auto') %把figure按照显示的大小保存
        %             end
        %             filenam = strcat(Fig_outfd0,'R',num2str(Ind_Rat(ns)),'-',date{ns},'-mean',num2str(ia));
        %             saveas(ffa,filenam,'bmp');
        %         end
        %=============同一个session内的prerunning，所有lap的CW和CCW点先连起来，再计算VFR===================%
        trial = unique(trial_use);
        ntrial = length(trial);
        VFR_session_all=cell(1,ndirs);VFR_trial = cell(1,ntrial) ;
        Sz_vel_session_all = cell(1,ndirs);Sz_vel_trial = cell(1,ntrial);
        vel_session_all = cell(1,ndirs);vel_trial = cell(1,ntrial);
        ind_trial = cell(1,ntrial);
        for nd = 1:ndirs
            for nt = 1:ntrial
                ind_trial{nt} = [ind_singlap{nd}{nt,1},ind_singlap{nd}{nt,2}];
                Sz_vel_trial{nt} = Sz_vel{1,1}(ind_trial{nt},:)';
                vel_trial{nt} = Vel_win(ind_trial{nt});
                
                ind11t = find(~isnan(mean(Sz_vel_trial{nt})));%移除nan
                Sz_vel_trial{nt}=Sz_vel_trial{nt}(:,ind11t);
                vel_trial{nt}=vel_trial{nt}(ind11t);
                
                ind22t = find(mean(Sz_vel_trial{nt})<=zscoremax);%移除zscore异常高的值
                Sz_vel_trial{nt}=Sz_vel_trial{nt}(:,ind22t);
                vel_trial{nt}=vel_trial{nt}(ind22t);
                
                for nb=1:nbin
                    VFR_trial{nt}(:,nb) = P_estimator(Sz_vel_trial{nt},vel_trial{nt},vr_center(nb),invh);
                end
                VFR_trial_all_day{nt}(:,:,ns)=VFR_trial{nt};
            end
            Sz_vel_session_all{nd} = Sz_vel{1,1}(ind_alllap{nd},:)';
            vel_session_all{nd} = Vel_win(ind_alllap{nd});
            
            ind11 = find(~isnan(mean(Sz_vel_session_all{nd})));%移除nan
            Sz_vel_session_all{nd}=Sz_vel_session_all{nd}(:,ind11);
            vel_session_all{nd}=vel_session_all{nd}(ind11);
            
            ind22 = find(mean(Sz_vel_session_all{nd})<=zscoremax);%移除zscore异常高的值
            Sz_vel_session_all{nd}=Sz_vel_session_all{nd}(:,ind22);
            vel_session_all{nd}=vel_session_all{nd}(ind22);
            
            for nb=1:nbin
                VFR_session_all{nd}(:,nb) = P_estimator(Sz_vel_session_all{nd},vel_session_all{nd},vr_center(nb),invh);
            end
            VFR_session_all_day{nd}(:,:,ns)=VFR_session_all{nd};
        end
        %============把test中正确的lap连起来再计算VFR，test中错误的也一样=======%
%         ind_test=cell(1,2); ind_test{1} = ind_test_cor; ind_test{2} = ind_test_err;
%         VFR_test=cell(1,2);  %test的正确和错误
%         Sz_vel_test=cell(1,2);
%         vel_test=cell(1,2);
%         for i = 1:length(ind_test)
%             VFR_test{i} = nan(length(freq_gamma),length(vr_center));
%             if ~isempty(ind_test{i})
%                 Sz_vel_test{i} = Sz_vel{1,1}(ind_test{i},:)';
%                 vel_test{i} = Vel_win(ind_test{i});
%                 
%                 ind11 = find(~isnan(mean(Sz_vel_test{i})));%移除nan
%                 Sz_vel_test{i}=Sz_vel_test{i}(:,ind11);
%                 vel_test{i}=vel_test{i};
%                 
%                 ind22 = find(mean(Sz_vel_test{i})<=zscoremax);%移除zscore异常高的值
%                 Sz_vel_test{i}=Sz_vel_test{i}(:,ind22);
%                 vel_test{i}=vel_test{i}(ind22);
%                 for nb=1:nbin
%                     VFR_test{i}(:,nb) = P_estimator(Sz_vel_test{i},vel_test{i},vr_center(nb),invh);
%                 end
%             end
%             VFR_test_all_day{i}(:,:,ns)=VFR_test{i};
%         end
        
        %===========画图（每圈先连起来CW和CCW功率谱，再计算VFR的）=======================%
        if plot_sign
            ffb = figure('unit','centimeters','position',[0 5 50.8 12]);
            for nt=1:ntrial
                subplot(1,ntrial,nt);
                imagesc(logvr(v1:v2), freq_gamma(g:end),VFR_trial{nt},[-0.5 0.5]);
                axis square
                colormap(jet(128)); % 图像默认是蓝色到黄色，加了这一句是蓝色到红色
                colorbar;
                axis xy; axis tight;
                set(gca, 'FontSize',16);%改变字体大小
                set(gca, 'XTick',tickx);
                set(gca, 'XTickLabel',2.^tickx*bin,'FontSize',14);
                xlabel('Running speed (cm/s)');
                ylabel('Frequency (Hz)');
                title(title_trial{nt});
                set(gcf,'PaperPositionMode','auto') %把figure按照显示的大小保存
%                 title(title_trial{nt});
            end
            filenam = strcat(resultFolder,'\VFR-result\Fig-VFR-Rat',num2str(Ind_Rat(ns)),'-',path_ns(21:end-1),'-all_trial3');
            saveas(ffb,filenam,'bmp');
            
%             %===================test正确错误的每一天的图============%
%             ffa = figure('unit','centimeters','position',[3 5 13*2 10]);
%             cl = [-0.4,0.4];%设置colorbar的范围
%             title_test={'correct','error'};
%             for i=1:length(VFR_test)
%                 subplot(1,2,i);
%                 imagesc(logvr(v1:v2), freq_gamma(g:end), VFR_test{i}(g:end,v1:v2-1),cl);
%                 colormap(jet(128)); colorbar;
%                 axis xy; axis tight;
%                 set(gca, 'FontSize',16);%改变字体大小
%                 set(gca, 'XTick',tickx);
%                 set(gca, 'XTickLabel',2.^tickx*bin);
%                 xlabel('Running speed (cm/s)');
%                 ylabel('Frequency (Hz)');
%                 title(['test-',title_test{i}]);
%                 set(gcf,'PaperPositionMode','auto') %把figure按照显示的大小保存
%             end
%             filenam = strcat(Fig_outfd0,'R',num2str(Ind_Rat(ns)),'-',date{ns},'-test-cor-err');
%             saveas(ffa,filenam,'bmp');
        end
        
        clear T_start Sz Sz_vel Vel_win
        clear data_pos_vel
        clear ffa ffb
        save([outFolder file_output],'VFR_trial','VFR_session_all','ntrial','freq_gamma');%保存当天的VFR
    end
end
save([resultFolder,'\VFR-result\Data_VFR_trials_allday_allgood.mat'],'VFR_trial_all_day','ntrial','freq_gamma')
%%
for nt = 1:ntrial
    VFR_trial_all_day_mean{nt} = mean(VFR_trial_all_day{nt}(:,:,[1,2,4:6,8,10:13,15,16]),3,'omitnan');%%[1,2,4,10,21,24]%(:,:,[1:11,16,20:22,25,26,28,29,30])(:,:,[2;4;5;6;8;10;11;16;20;22;25]')[1:11,16,20:22,25,26,28,29,30]
end
ffa = figure('unit','centimeters','position',[0 5 50.8 12]);
for nt=1:ntrial
    subplot(1,ntrial,nt);
    imagesc(logvr(v1:v2), freq_gamma(g:end),VFR_trial_all_day_mean{nt},[-0.5 0.5]);
    axis square
    colormap(jet(128)); % 图像默认是蓝色到黄色，加了这一句是蓝色到红色
    colorbar;
    axis xy; axis tight;
    set(gca, 'FontSize',16);%改变字体大小
    set(gca, 'XTick',tickx);
    set(gca, 'XTickLabel',2.^tickx*bin,'FontSize',14);
    xlabel('Running speed (cm/s)');
    ylabel('Frequency (Hz)');
    title(title_trial{nt});
    set(gcf,'PaperPositionMode','auto') %把figure按照显示的大小保存
    
end
cd('H:\neuralynx\gamma in sequence result\VFR-result')
saveas(gcf,'VFR_allgood_f2lap','epsc')
saveas(gcf,'VFR_allgood_f2lap','bmp')
%%
% % close all
% % 
% % %==========pre，semple，test各自总的平均的VFR===========%
% % VFR_lap_ns = {};%所有session的平均功率谱，样本是lap
% % for nd = 1: ndirs
% %     VFR_lap_ns{ia,nd} = nanmean(VFR0_ns{ia}{1,nd},3);
% %     VFR_session_all_day_mean{nd} = nanmean(VFR_session_all_day{nd},3);%
% % end
% % %=======test的正确和错误的平均值================%
% % VFR_test_ns={};
% % for i=1:2
% %     VFR_test_ns{i}=nanmean(VFR_test_all_day{i},3);
% % end
% % 
% % %====画出pre，semple，test所有天，所有圈平均的VFR，并保存===%
% % ffa = figure('unit','centimeters','position',[3 5 13*3 10]);
% % cl = [-0.4,0.4];%设置colorbar的范围
% % for nd=1:ndirs
% %     subplot(1,3,nd);
% %     imagesc(logvr(v1:v2), freq_gamma(g:end), VFR_lap_ns{ia,nd}(g:end,v1:v2-1),cl);
% %     colormap(jet(128)); colorbar;
% %     axis xy; axis tight;
% %     set(gca, 'FontSize',18);%改变字体大小
% %     set(gca, 'XTick',tickx);
% %     set(gca, 'XTickLabel',2.^tickx*bin);
% %     xlabel('Running speed (cm/s)');
% %     ylabel('Frequency (Hz)');
% %     title([title_session{nd},'-mean']);
% %     set(gcf,'PaperPositionMode','auto') %把figure按照显示的大小保存
% % end
% % filenam = strcat(Fig_outfd,'before-reward-',note,'-mean');
% % saveas(ffa,filenam,'bmp');
% % %     saveas(ffa,[filenam,'.eps'],'psc2');
% % %================================功率谱先连起来再算VFR，所有天的平均结果=====================================%
% % 
% % a1 = find(freq_gamma>45 & freq_gamma<55);
% % for nd=1:ndirs
% %     ffb = figure('unit','centimeters','position',[3 5 11.68 9.75]);
% %     a=VFR_session_all_day_mean{nd}(g:end,v1:v2-1);
% %     a(a1,:)=nan;
% %     imagesc(logvr(v1:v2), freq_gamma(g:end),a,cl);
% %     colormap(jet(128)); colorbar;
% %     axis xy; axis tight;
% %     set(gca, 'FontSize',18);%改变字体大小
% %     set(gca, 'XTick',tickx);
% %     set(gca, 'XTickLabel',2.^tickx*bin);
% %     xlabel('Running speed (cm/s)');
% %     ylabel('Frequency (Hz)');
% %     title([title_session{nd},'-all']);
% %     set(gcf,'PaperPositionMode','auto') %把figure按照显示的大小保存
% %     filenam = strcat(Fig_outfd,'before-reward-',note,'-all-',title_session{nd});
% %     saveas(ffb,filenam,'bmp');
% %     saveas(ffb,[filenam,'.eps'],'psc2');
% % end
% % 
% % 
% % %===================test正确错误的总图============%
% % ffa = figure('unit','centimeters','position',[3 5 13*2 10]);
% % cl = [-0.4,0.4];%设置colorbar的范围
% % title_test={'correct','error'};
% % for i=1:length(VFR_test_ns)
% %     subplot(1,2,i);
% %     imagesc(logvr(v1:v2), freq_gamma(g:end), VFR_test_ns{i}(g:end,v1:v2-1),cl);
% %     colormap(jet(128)); colorbar;
% %     axis xy; axis tight;
% %     set(gca, 'FontSize',18);%改变字体大小
% %     set(gca, 'XTick',tickx);
% %     set(gca, 'XTickLabel',2.^tickx*bin);
% %     xlabel('Running speed (cm/s)');
% %     ylabel('Frequency (Hz)');
% %     title(['test-',title_test{i}]);
% %     set(gcf,'PaperPositionMode','auto') %把figure按照显示的大小保存
% % end
% % filenam = strcat(Fig_outfd,'before-reward-',note,'-test-cor-err');
% % saveas(ffa,filenam,'bmp');
% % saveas(ffa,[filenam,'.eps'],'psc2');
% % %=================================================================%
% % if icon==1
% %     VFR_Con = VFR;%每一圈的功率谱
% %     VFR_session_Con = VFR0_ns;%把3个session中的各圈平均，数据按照每一天为单位保存
% %     VFR_mean_Con = VFR_lap_ns;%把所有天的lap平均
% %     VFR_session_all_day_Con = VFR_session_all_day;
% %     VFR_session_all_day_mean_Con = VFR_session_all_day_mean;
% %     save(file_out,'VFR_Con','VFR_mean_Con','VFR_session_Con','VFR_session_all_day_Con','VFR_session_all_day_mean_Con')%保存所有实验天的vfr，以及平均值
% % else
% %     VFR_AD = VFR;
% %     VFR_session_AD = VFR0_ns;
% %     VFR_mean_AD = VFR_lap_ns;
% %     VFR_session_all_day_AD = VFR_session_all_day;
% %     VFR_session_all_day_mean_AD = VFR_session_all_day_mean;
% %     save(file_out,'VFR_AD','VFR_mean_AD','VFR_session_AD','VFR_session_all_day_AD','VFR_session_all_day_mean_AD','-append')%保存所有实验天的vfr，以及平均值
% % end
% % close all
% 
% %%
% clear
% %频段和速度上平均
% % file_outfd0 = 'H:\group data\VFR_image\trail_eight';
% % file_input = 'VFR_before_reward_8trail.mat';
% % file_input_fre = 'Data_Z_gamma_multitaper_0p5s_8trail.mat';
% % file_output_mean = {'mean_vfr_con_before_reward_8trail';'mean_vfr_AD_before_reward_8trail'};
% file_outfd0 = 'H:\group data\VFR_image\all_trail';
% file_input = 'VFR_before_reward_all_trail_sam1_test0_1.mat';
% file_input_fre = 'Data_Z_gamma_multitaper_all_trail_sam1_test0_1.mat';
% file_output_mean = {'mean_vfr_con_before_reward_all_trail';'mean_vfr_AD_before_reward_all_trail'};
% 
% vr1=0.65;
% vr2=35.05;
% logvr = 0:0.2:6.2;
% bin=0.5;
% vr=2.^logvr*bin;
% vr_center = mean([vr(1:end-1);vr(2:end)]);
% [~,v1]=min(abs(vr_center-vr1));
% [~,v2]=min(abs(vr_center-vr2));
% 
% for icon = 1:2
%     %每一圈的VFR
%     pre_VFR={};    pre_VFR_all=[];
%     sample_VFR={}; sample_VFR_all=[];
%     test_VFR={};   test_VFR_all=[];
%     %每个session内每一圈连起来的VFR
%     pre_VFR_2=[];    pre_VFR_all_2=[];
%     sample_VFR_2=[]; sample_VFR_all_2=[];
%     test_VFR_2=[];   test_VFR_all_2=[];
%     
%     if icon == 1
%         cd('E:\matlab\bin\code for CT\directories');
%         directories_allData_v2%control
%         
%     else
%         cd('E:\matlab\bin\code for CT\directories');
%         directories_allData_v3 %AD
%         
%     end
%     isession = size(path,1);
%     ndirs = 3;
%     mean_VFR = cell(ndirs,1 );
%     mean_VFR_all = cell(ndirs, 20);%sample_test最多跑20次
%     theta_low = 6;
%     theta_high = 12;
%     gamma_slow_low = 25;
%     gamma_slow_high = 45;
%     gamma_fast_low = 55;
%     gamma_fast_high = 100;
%     for ns = 1:isession
%         path_ns = path{ns};
%         cd(path_ns);
%         VFR_ns = strcat(path_ns,file_input);
%         freq_ns = strcat(path_ns,file_input_fre);
%         load(VFR_ns, 'VFR_nd','VFR_exp','VFR_session_all');
%         load(freq_ns, 'freq_gamma');
%         trackdata_ns = trackdata{ns};
%         load(trackdata_ns,'n_prerunning','n_sample_test','sign_correct_test','sign_correct_sample');
%         load('Data_angle_ontrack.mat','n_test');
%         speed_lim = v1:v2;
%         
%         theta_num = find(freq_gamma>=theta_low & freq_gamma<=theta_high);
%         gamma_slow_num = find(freq_gamma>=gamma_slow_low & freq_gamma<=gamma_slow_high);
%         gamma_fast_num = find(freq_gamma>=gamma_fast_low & freq_gamma<=gamma_fast_high);
%         ind_fre={theta_num,gamma_slow_num,gamma_fast_num};
%         ntrail_session=[]; aa=[];
%         for nd= 1:ndirs
%             for nl=1:length(VFR_exp{1,1})
%                 aa(nd,nl) = ~isempty(VFR_exp{1,1}{nd,nl});
%             end
%             ntrail_session(nd) = sum(aa(nd,:));
%         end
%         
%         for nd = 1:ndirs
%             if nd == 1
%                 for nl = 1:ntrail_session(nd)
%                     for fre=1:length(ind_fre) %三个频段
%                         pre_VFR {ns,1}(nl,fre) = nanmean(nanmean(VFR_exp{1,1}{nd,nl}(ind_fre{fre},speed_lim)));
%                     end
%                 end
%             else
%                 for nl = 1:ntrail_session(nd) %8个lap
%                     theta = NaN;
%                     gamma_slow = NaN;
%                     gamma_fast = NaN;
%                     if ~isempty(VFR_exp{1,1}{nd,nl})
%                         if sign_correct_sample(nl) == 1&(sign_correct_test(nl)==0||sign_correct_test(nl)==1)
%                             theta = nanmean(nanmean(VFR_exp{1,1}{nd,nl}(theta_num,:)));
%                             gamma_slow = nanmean(nanmean(VFR_exp{1,1}{nd,nl}(gamma_slow_num,:)));
%                             gamma_fast = nanmean(nanmean(VFR_exp{1,1}{nd,nl}(gamma_fast_num,:)));
%                         end
%                     end
%                     if nd==2
%                         sample_VFR {ns,1}(nl,1)= theta;
%                         sample_VFR {ns,1}(nl,2)= gamma_slow;
%                         sample_VFR {ns,1}(nl,3)= gamma_fast;
%                     else
%                         test_VFR {ns,1}(nl,1)= theta;
%                         test_VFR {ns,1}(nl,2)= gamma_slow;
%                         test_VFR {ns,1}(nl,3)= gamma_fast;
%                     end
%                 end
%             end
%             theta_ns = nanmean(nanmean(VFR_nd{1,nd}(theta_num,:)));%每天实验的平均
%             gamma_slow_ns = nanmean(nanmean(VFR_nd{1,nd}(gamma_slow_num,:)));
%             gamma_fast_ns = nanmean(nanmean(VFR_nd{1,nd}(gamma_fast_num,:)));
%             mean_VFR{nd,1}(1,ns) = theta_ns;
%             mean_VFR{nd,1}(2,ns) = gamma_slow_ns;
%             mean_VFR{nd,1}(3,ns) = gamma_fast_ns;
%         end
%         sample_VFR_all = [sample_VFR_all;sample_VFR{ns,1}]; %把每天每一圈的连起来，做统计
%         test_VFR_all = [test_VFR_all;test_VFR{ns,1}];
%         pre_VFR_all = [pre_VFR_all;pre_VFR{ns,1}];
%         
%         %===========================先把功率谱连起来，再算VFR=======================%
%         for fre=1:length(ind_fre)
%             pre_VFR_2 (ns,fre) = nanmean(nanmean(VFR_session_all{1}(ind_fre{fre},speed_lim))); %每个session内每一圈连起来的VFR
%             sample_VFR_2 (ns,fre) = nanmean(nanmean(VFR_session_all{2}(ind_fre{fre},speed_lim))); %每个session内每一圈连起来的VFR
%             test_VFR_2 (ns,fre) = nanmean(nanmean(VFR_session_all{3}(ind_fre{fre},speed_lim))); %每个session内每一圈连起来的VFR
%         end
%     end
%     %把AD和con的数据分别存起来
%     if icon == 1
%         mean_VFR_Con = mean_VFR;
%         sample_all_Con = sample_VFR_all;%每一天中的每一圈，不按天分开
%         test_all_Con = test_VFR_all;
%         pre_all_Con = pre_VFR_all;
%         sample_Con = sample_VFR;% 每一天中的每一圈，按天分开
%         test_Con = test_VFR;
%         pre_Con = pre_VFR;
%         pre_VFR_2_Con=pre_VFR_2;
%         sample_VFR_2_Con=sample_VFR_2;
%         test_VFR_2_Con=test_VFR_2;
%         cd(file_outfd0);
%         save(file_output_mean{icon},'sample_all_Con','test_all_Con','pre_all_Con','sample_Con','test_Con','pre_Con',...
%             'mean_VFR_Con','pre_VFR_2_Con','sample_VFR_2_Con','test_VFR_2_Con');
%     else
%         mean_VFR_AD = mean_VFR;
%         sample_all_AD = sample_VFR_all;
%         test_all_AD = test_VFR_all;
%         pre_all_AD = pre_VFR_all;
%         sample_AD = sample_VFR;
%         test_AD = test_VFR;
%         pre_AD = pre_VFR;
%         pre_VFR_2_AD=pre_VFR_2;
%         sample_VFR_2_AD=sample_VFR_2;
%         test_VFR_2_AD=test_VFR_2;
%         cd(file_outfd0);
%         save(file_output_mean{icon},'sample_all_AD','test_all_AD','pre_all_AD','sample_AD','test_AD','pre_AD',...
%             'mean_VFR_AD','pre_VFR_2_AD','sample_VFR_2_AD','test_VFR_2_AD');
%     end
% end
