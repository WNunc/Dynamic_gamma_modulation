% crate by WN on 2022/08/08
% 参考 ZN decode_replay_use_sequence_ZN_v7.m中的SWR检测
% 找到所有SWR，并找到rest期间（60s）的SWR的开始和结束
% 计算这些SWR发生期间相锁神经元和非相锁神经元的放电率（replay fire rate）
%%
clear
close all
directories_allData_v0
resuletFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
inputName = 'data_fireRate_exclude_cell.mat';
outputName = 'data_fireRate_exclude_cell_SWR_1min.mat';% 没有V的是3std，V2是2std,1min是用休息时检测的SWR,1std
peak0 = 1;
gap =1;% tov detect SWRs
sampfreq = 2000;% EEG sapmle rate
pre_use=[1 2 3 4 5 6 7 8 9 10]; % 圈数
D = [1 2 1 2 1 2 1 2 1 2];% 奇数为cw，偶数为ccw
trial_use = [1 1 2 2 3 3 4 4 5 5];% trial数

nlaps = length(pre_use);
ntrial = length(unique(trial_use));
for ns = 1:isession
    path_ns = path{ns};
    disp(path_ns);cd(path_ns);
    outFolder = [resuletFolder path_ns(13:end)];
    inFolder = outFolder;
    % 载入排除掉的spike
    if ~exist([inFolder inputName],'file')
        disp('escape this day')
        continue
    end
    load([inFolder inputName])
    % 载入trackdata
    trackdata_ns = trackdata{ns};
    load(trackdata_ns,'Ts_prerunning')
    %% 取SWRs的时间点
    % 载入EEG
    Csclist = CSClist_CA1{ns};
    nCSC= length(Csclist);
    all_removeeg = {};
    if ~isempty(Csclist) %确定csclist是否为空
        %所有remove的eeg点
        EEG = cell(0);
        time_EEG=cell(0);
        
        for ncsc = 1:nCSC
            filnam = strcat(path_ns,'CSC',num2str(Csclist(ncsc)),'.ncs');
            if(exist(filnam)~=2)
                continue
            end
            % read EEG
            [sample,tt]=loadCSC_new_cz(filnam);
            EEG{ncsc}=sample;
            time_EEG{ncsc}=tt;
            % kick out time points
            [removinds,~] = fixeeg_cz(sample,gap);  % find out-of-range EEG
            removeeg=unique(removinds);
            all_removeeg{ncsc}=removeeg;
        end
        choseeeg = 2;
        time_EEG_chosen = time_EEG{choseeeg};%挑选其中一导的时间
    end
    
    %     aaa = [EEG{1:7}];figure;
    %     for i = 1:length(RippleOnsetIndex)
    %         plot(aaa(RippleOnsetIndex(i):RippleOffsetIndex(i),:))
    %         pause(0.2)
    %     end
    
    % 找到prerunning休息期间的SWR
    TS_pre_rest = cell(ntrial,2);
    RippleIndex_rest = {};
    RippleIndex_trial = {};
    rate_fg = {};rate_nfg = {};
    RATE_fg = {};RATE_nfg = {};
    spk_fg = cell(5,2);spk_nfg = cell(5,2);
    EEGind_rest = cell(5,2);
    RATE_trial = [];RATE_trial0 = [];
    for nt = 1:ntrial
        for D = 1:2 %CW和CCW
            % rest的开始结束
            if D == 1
                rest_strat = Ts_prerunning{1}(nt,2)./1000000;
                rest_end = Ts_prerunning{2}(nt,1)./1000000;
            elseif nt < 5
                rest_strat = Ts_prerunning{2}(nt,2)./1000000;
                rest_end = Ts_prerunning{1}(nt+1,1)./1000000;
            else
                rest_strat = Ts_prerunning{2}(nt,2)./1000000;
                rest_end = rest_strat + 30;
            end
            % rest期间的ripples
            TS_pre_rest{nt,D} = [rest_strat,rest_end];
            
            ind_r = find(time_EEG_chosen>rest_strat  &  time_EEG_chosen<rest_end);
            EEGind_rest{nt,D} = ind_r;
        end
        % 1min EEG
        ind_R = [EEGind_rest{nt,1};EEGind_rest{nt,2}];
        EEG_R = cell2mat(EEG);
        EEG_R = EEG_R(ind_R,:);
        [smp,nEEG] = size(EEG_R);
        nEEG = ones(1,nEEG);
        EEG_R = mat2cell(EEG_R,smp,nEEG);
        EEGtime_R = time_EEG_chosen(ind_R);
        % 挑选休息期间的1minEEG 检测SWR
        [RippleOnsetIndex, RippleOffsetIndex] = DetectRipples_v4(EEG_R,sampfreq);%SWR时长大于15ms,std=2 对所有导找SWRs
        SWR_TS = [EEGtime_R(RippleOnsetIndex),EEGtime_R(RippleOffsetIndex)];
%         aaa = [EEG_R{1:7}];figure;
%         for i = 1:length(RippleOnsetIndex)
%             plot(aaa(RippleOnsetIndex(i):RippleOffsetIndex(i),:))
%             pause(0.05)
%         end
        % ripples期间的两类神经元的放电率
        for D = 1:2
            ind_swr = find(SWR_TS(:,1)>TS_pre_rest{nt,D}(1) & SWR_TS(:,2)<TS_pre_rest{nt,D}(2));
            RippleIndex_rest{nt,D} = ind_swr;
            for c = 1:size(SPK{1,D},1)
                Num_spk_fg = 0;Num_spk_nfg = 0;interval = 0;
                for ts = 1:length(ind_swr)
                    rippleTs_onset =  SWR_TS(ind_swr(ts),1);
                    rippleTs_offset = SWR_TS(ind_swr(ts),2);
                    interval = interval + (rippleTs_offset - rippleTs_onset);
                    
                    spk_fg{nt,D}(ts,c) = length(find( SPK{1,D}{c,1}>rippleTs_onset &SPK{1,D}{c,1}<rippleTs_offset));
                    spk_nfg{nt,D}(ts,c) = length(find( SPK{2,D}{c,1}>rippleTs_onset &SPK{2,D}{c,1}<rippleTs_offset));
                    
                    Num_spk_fg = Num_spk_fg + spk_fg{nt,D}(ts,c);
                    Num_spk_nfg = Num_spk_nfg + spk_nfg{nt,D}(ts,c);
                end
                rate_fg{nt,D}(c) = Num_spk_fg/interval;
                rate_nfg{nt,D}(c) = Num_spk_nfg/interval;
            end
        end
        
        RATE_fg{nt} = [rate_fg{nt,1},rate_fg{nt,2}];
        RATE_nfg{nt} = [rate_nfg{nt,1},rate_nfg{nt,2}];
        A = RATE_fg{nt};B = RATE_nfg{nt};
        
        RATE_trial(nt,1) = mean(A(A~=0),'omitnan');
        RATE_trial(nt,2) = mean(B(B~=0),'omitnan');
        RATE_trial0(nt,1) = mean(A,'omitnan');
        RATE_trial0(nt,2) = mean(B,'omitnan');
    end
    
    save([outFolder,outputName],'time_EEG_chosen','RippleIndex_rest',...
        'rate_fg','rate_nfg','RATE_trial','RATE_trial0','spk_fg','spk_nfg')
    
    
end

%%
clear
inputName = 'data_fireRate_exclude_cell_SWR_1min.mat';
resuletFolder = 'H:\neuralynx\gamma in sequence result';
directories_allData_v0_allgood
nsn = 0;
RATE_fg = [];RATE_nfg = [];
for ns = [1,2,4,10,14]%[1,2,4:6,8,10:13,15,16]%1:isession %
    nsn = nsn + 1;
    path_ns = path{ns};
    disp(path_ns);
    outFolder = [resuletFolder path_ns(13:end)];
    inFolder = outFolder;
    load([inFolder,inputName])
    RATE_fg = [RATE_fg, cell2mat(rate_fg)];
    RATE_nfg = [RATE_nfg, cell2mat(rate_nfg)];
    RATE_trial(isnan(RATE_trial)) = 0;
    RATE_allday(:,:,nsn) = RATE_trial;
end
RATE_fg (:,mean(RATE_fg,'omitnan')==0) = [];
RATE_nfg (:,mean(RATE_nfg,'omitnan')==0) = [];
RATE_fg (RATE_fg==0) = nan;
RATE_nfg (RATE_nfg==0) = nan;


Rate_FG = squeeze(RATE_allday(:,1,:));
Rate_NFG = squeeze(RATE_allday(:,2,:));

Rate_FG_m = mean(Rate_FG,1);
Rate_NFG_m = mean(Rate_NFG,1);