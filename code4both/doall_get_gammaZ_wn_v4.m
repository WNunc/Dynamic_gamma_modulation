% edited on 2021-03-17,wxl
% 先算出功率谱，然后截取每一个lap的窗口内的功率谱，拼接起来，再进行zscore
% use multitaper power spectrum in Chronux
% z-score TFR

% 2021-12-20 修改自doall_get_gammaZ_wxl_v2
% 仅使用sample-test的前八个trail，只使用pre的后四个trail，也就是pre，sample，test各自都是8圈
% zscore的范围包括：pre的后8圈，sample前8圈中正确的，test前8圈中正确/错误的，没停的不要，
% 而且只采用到奖励点之前的
% 2022-1-10 修改，还是用整圈的做tfr，仅仅保留sample正确，而且test有停下的trail
%
% 2022-8-8 修改自doall_get_gammaZ_wxl_v3
% 适配 gamma within the sequence 程序
%%
clear
close all
directories_allData_v0

ndirs=1;
% pre_use=[2 3 4 5 7 8 9 10]; %1-5为cw，6-10为ccw，不用pre的第一个trail，即1和6
% N_trail=8; %仅使用sample-test的前八个trail
%file_output = 'Data_Z_gamma_multitaper_0p5s_8trail.mat';
pre_use=[1 2 3 4 5 6 7 8 9 10]; %单数为cw，双数为ccw
N_trail=20;
% 仅保留 prerunning 
file_output = 'Data_Z_gamma_multitaper_all_trail_prerunning3.mat'; % 2  = speedlimit0.68；3 = speedlimit0.68
file_input_speed = 'Data_angle_ontrack.mat';  % input the speed information
% use all tetrodes with cell firing
CSClist = CSClist_CA1;
Fs = 2000;
% speedlimit = 0.75;%师兄的
speedlimit = 5.08; %%%%%%
speedmax = 35.03;
% set window width
win0 = 0.5;
winstep0 = 0.2;
% 生成multitaper参数
padval=2;  %additional padding orders of two in the spectrogram
param = struct('tapers',[3,5], 'Fs',Fs, 'pad', padval);
g1 = 4;
g2 = 100;
isession = size(path,1);
narea = size(CSClist,2);
resultFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
for ns = 1:isession
    
    path_ns = path{ns};cd(path_ns);disp(pwd)
    goodphase_ns = Seq_cutPhase{ns,1};%***
    subfolder1 = Directionfolder{1};
    subfolder2 = Phasecutfolder;
    outFolder = [resultFolder path_ns(13:end)];
    case3 = num2str(goodphase_ns);
    trackdata_ns = trackdata{ns};
    load(trackdata_ns,'n_prerunning','n_sample_test','ts_prerunning_CCW_end','ts_prerunning_CCW_start',...
        'ts_prerunning_CW_end','ts_prerunning_CW_start','ts_sample_end','ts_sample_start',...
        'ts_test_end','ts_test_start','ts_sample_reward','ts_test_reward','sign_correct_sample','sign_correct_test');
    Csclist_ca1 = CSClist{ns,1};
    
    file_input1 = strcat(path_ns,subfolder2{1},'data_AllLap_thetaphasecut_',case3,'.mat');
    TS_CW = load(file_input1, 'TS_running');%越过第1个点为start，越过第18个点为end（逆时针相反）
    TS_CW = TS_CW.TS_running;
    file_input2 = strcat(path_ns,subfolder2{2},'data_AllLap_thetaphasecut_',case3,'.mat');
    TS_CCW = load(file_input2, 'TS_running');%越过第1个点为start，越过第18个点为end（逆时针相反）
    TS_CCW = TS_CCW.TS_running;
    TS_running = [TS_CW;TS_CCW];
    TS_running = TS_running([1,6,2,7,3,8,4,9,5,10],:);
    
    load([subfolder1,file_input_speed],'data_angle_all');
    
    disp('Getting power info');
    
    % find out EEG time points jumping out of the range [-2000,2000]uV
    gap = 1;  % remove window length
    csclist_ns = cell2mat(CSClist(ns,:));
    removinds_win = [];
    Vel_win = [];
    win_keep = [];
    win_keep_vel = [];
    Sz_vel_csc = {};
    
    removinds=[];
    for ncsc = 1:length(csclist_ns)
        file0 = strcat(path_ns,'CSC',num2str(Csclist_ca1(ncsc)),'.ncs');
        [e, eTSi, eTS, Fs_raw] = loadCSC_new_cz_allrecording(file0);
        [removinds_all0,windows0] = fixeeg_cz(e,gap);
        if ~isempty(removinds_all0)
            removinds = [removinds,removinds_all0];
        end
    end
    removinds = unique(removinds);
    
    % convert to windows numbers
    win1 = 1:round(winstep0*Fs):length(e)-round(win0*Fs);
    win2 = round(win0*Fs):round(winstep0*Fs):length(e);
    nwin = min([length(win1),length(win2)]);
    removinds_win0 = [];
    
    for iw = 1:nwin
        ind0 = find((removinds >= win1(iw) & removinds <= win2(iw)));
        if ~isempty(ind0)
            removinds_win0 = [removinds_win0,iw];%想要去除的功率谱的时间下标
        end
        % get mean running speed in each window (脑电一个窗里的所有速度取平均值，作为此时间窗的速度指标)
        ind0 = find(data_angle_all(:,1) >= eTSi(win1(iw)) & data_angle_all(:,1) <= eTSi(win2(iw)));
        Vel_win(iw,1) = mean(data_angle_all(ind0,3));%每一个时间窗内的速度
    end
    
    % save the window numbers
    win_keep = setdiff(1:nwin,removinds_win0);%能用的脑电的功率谱下标
    ind = find(Vel_win >= speedlimit & Vel_win <= speedmax);
    win_keep_vel = intersect(win_keep,ind);%满足脑电幅度要求，行为学速度要求的功率谱的下标
    
    
    % get power for gamma rhythms
    T_start = [];
    Sz = cell(narea,1);
    Sz_vel = cell(narea,1);
    for ia = 1:narea
        csclist_ia = CSClist{ns,ia};
        if ~isempty(csclist_ia)
            Sz_na = [];
            Sz_vel_na = [];
            e=[];
            eTSi=[];
            if n_sample_test >N_trail
                n_sample_test=N_trail;
            end
            for ncsc = 1:length(csclist_ia)
                file0 = strcat(path_ns,'CSC',num2str(Csclist_ca1(ncsc)),'.ncs');
                [e, eTSi, eTS, Fs_raw] = loadCSC_new_cz_allrecording(file0);
                eeg_end=find(eTSi<TS_CCW(5,2),1,'last'); %仅对5圈pre进行分析
                e=e(1:eeg_end+2000*60);
                eTSi=eTSi(1:eeg_end+2000*60);
                S = [];
                S_vel = [];
                nwin = [];
                
                % get power spectrum on raw EEG
                [S0, ~, freq0] = mtspecgramc(e', [win0 winstep0], param); %200 ms bins with no overlap, as in Ahmed and Mehta 2012
                Nwin=round(Fs*win0); % number of samples in window
                Nstep=round(winstep0*Fs); % number of samples to step through
                winstart=1:Nstep:length(eTSi)-Nwin+1;
                T0 = eTSi(winstart)';
                T_start = T0;
                % all = [ind_all,ind];
                
                %==========================保留整圈（都是从路过第一个点为开始）=====================%
                ind_all=[];
                for nd=1:ndirs % nd=1,prerunning; nd=2,semple; nd=3,test
                    if nd==1
                        for nlap=1:length(pre_use)
                            nlap_use=pre_use(nlap);
                            ind = find(T_start>=TS_running(nlap_use,1) & T_start<=TS_running(nlap_use,2));
                            ind_all = [ind_all,ind];
                        end
%                     else
%                         if nd==2
%                             for nlap=1:n_sample_test %仅仅保留sample正确，而且test有停下的trail，否则vfr为空
%                                 if sign_correct_sample(nlap)==1 && ~isnan(sign_correct_test(nlap))
%                                     ind = find(T_start>=ts_start_stop_point{1,nd}(nlap,1) & T_start<=ts_start_stop_point{1,nd}(nlap,2));
%                                 else
%                                     ind = [];
%                                 end
%                                 ind_all = [ind_all,ind];
%                             end
%                         else
%                             for nlap=1:min([n_sample_test,length(ts_start_stop_point{1,nd})])
%                                 if sign_correct_sample(nlap)==1 && ~isnan(sign_correct_test(nlap))%仅仅保留sample正确，而且test有停下的trail，否则vfr为空
%                                     ind = find(T_start>=ts_start_stop_point{1,nd}(nlap,1) & T_start<=ts_start_stop_point{1,nd}(nlap,2));
%                                 else
%                                     ind = [];
%                                 end
%                                 ind_all = [ind_all,ind]; %所有的test（sample正确时候）
%                             end
%                         end
                    end
                end
                
                %=================================================================================================================================%
                ind = find(freq0>g1 & freq0<g2);
                S0 = S0(:,ind);
                freq_gamma = freq0(ind);
                
                
                S0_keep = nan(size(S0));
                win_keep_use = find(win_keep<=length(S0),1,'last');
                S0_keep(win_keep(1:win_keep_use),:) = S0(win_keep(1:win_keep_use),:);%满足幅值要求的eeg的功率谱
                
                S0_keep_vel = nan(size(S0));
                win_keep_vel_use = find(win_keep_vel<=length(S0),1,'last');
                S0_keep_vel(win_keep_vel(1:win_keep_vel_use),:) = S0(win_keep_vel(1:win_keep_vel_use),:);
                
                S = nan(size(S0));
                S(ind_all,:) = S0_keep(ind_all,:); % 8个pre，8个sample-test的时间点
                S_vel = nan(size(S0));
                S_vel(ind_all,:) = S0_keep_vel(ind_all,:);%满足速度要求的eeg的功率谱，有上限也有下限
                nwin = size(S0,1);
                
                if isempty(Sz_na) && isempty(Sz_vel_na)  %把每个导联的功率谱都加起来
                    Sz_na = nanzscore_cz(S,0,1);
                    Sz_vel_na = nanzscore_cz(S_vel,0,1);
                else
                    Sz_na = Sz_na+nanzscore_cz(S,0,1);
                    Sz_vel_na = Sz_vel_na+nanzscore_cz(S_vel,0,1);
                end
                Sz_vel_csc{ncsc} = nanzscore_cz(S_vel,0,1);%把每一导的功率谱存起来
            end
            Sz_na = Sz_na./length(csclist_ia);%一次实验的功率谱平均值
            Sz_vel_na = Sz_vel_na./length(csclist_ia);
            
            Sz{ia,1} = Sz_na(1:nwin,:);
            Sz_na(1:nwin,:) = [];
            Sz_vel{ia,1} = Sz_vel_na(1:nwin,:);
            Sz_vel_na(1:nwin,:) = [];
        end
    end
    
%     save([outFolder,file_output],'T_start','Sz','Sz_vel','Sz_vel_csc','Vel_win','freq_gamma','win0','winstep0','param','-v7.3');
    save([outFolder,file_output],'T_start','Sz','Sz_vel','Sz_vel_csc','Vel_win','freq_gamma','win0','winstep0','param','TS_CW','TS_CCW','TS_running');
    clear Sz Sz_vel T_start Sz_vel_csc
    fclose all;
    
    cd ../
end