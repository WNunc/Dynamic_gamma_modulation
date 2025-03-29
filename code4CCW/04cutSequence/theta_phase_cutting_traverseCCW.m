%% cut theta phase on spike min phase
% histogram use all spike
% call the thetawindows_all_CT_v4_EEG_wn get all spikes' phase

clear

close all
directories_allData_v0
subfolder = 'Tseq-CCW\';
result_stat = {};
for ns = 14 %1:isession
    path_ns = path{ns};
    cd(path_ns);
    disp(path_ns)
    csclist = csclist_CA1_mid{ns};%用bundle中心的导
    csclistind = find(CSClist_CA1{ns}==csclist);
    ttlist = [path_ns 'TTList_dCA1_pyr.txt'];
    load(trackdata{ns}, 'Ang_RewardLoc_ontrack')
    load([path{ns} subfolder 'Data_angle_ontrack.mat'])
    mkdir('CuttingPhaseCCW')
    % Load scores
    Scores_input = [path_ns subfolder 'BayesData_CircMap_wn_dt40ms_5ms_1cells_1.mat'];
    load(Scores_input)
    scores = Scores;
    nseg = 1;%哪个部分？ 1=pre-running, 2=sample, 3=test, 4=post
    
    for nl = 1:5
        lap_ts = scores{nseg}{nl,6};%行为学时间
        lap_pos = scores{nseg}{nl,5}(scores{nseg}{nl,4});%实际位置
        lap_speed = scores{nseg}{nl,8}(1,:);
        lap_ts = lap_ts(1:end-20);
        lap_pos = lap_pos(1:end-20);
        lap_speed = lap_speed(1:end-20);
        [TSlimit, speed, speed_ang]= TSfind_laplimit(data_angle,Ang_RewardLoc_ontrack,lap_speed,lap_pos,lap_ts,nseg,nl);
%         TS_running(nl,1:2) = TSlimit;
    end
    
    [windows,~,name_tet,ph1,ph2] = thetawindows_all_CT_v4_EEG_wn(csclist,ttlist,path_ns,TS_running);
    % min_phase = ph1;% all session
    min_phase = ph2;
    disp(['the minimum spike firing phase is ' num2str(min_phase)])
    start_phase = min_phase;
    mth = csclistind;% bundle中心的导对应的编号
    thetacycle = [];
    traverseIND = 0;
    phasestep = 5;
    
    while traverseIND < (360/phasestep)
        cut_phase = start_phase +  traverseIND * phasestep;
        cut_phase = mod(cut_phase,360)
        traverseIND = traverseIND + 1;
        
        for nl = 1:5 % nl = 2;%哪一圈？
            lap_ts = scores{nseg}{nl,6};%行为学时间
            lap_ang = scores{nseg}{nl,5};%位置（角度）
            lap_pxn = scores{nseg}{nl,3};%解码位置（Pxn）
            lap_pos = scores{nseg}{nl,5}(scores{nseg}{nl,4});%实际位置
            lap_EEG = scores{nseg}{nl,11};%EEG数据
            lap_EEG_ts = scores{nseg}{nl,10};%EEG时间
            lap_theta = scores{nseg}{nl,18};%scores中的theta
            lap_thetaphs = scores{nseg}{nl,21};%scores中的theta相位
            lap_timepoint_track = size(lap_pxn,2);%路径总时间点
            lap_timepoint_EEG = size(lap_EEG,2);%EEG总时间点
            lap_speed = scores{nseg}{nl,8}(1,:);
            %% 运动时的时间
            limit = TS_running(nl,:);
            %% 用spike发放最少的相位切割theta
            % 对应的时间点
            
            % [~,I] = min(abs(lap_thetaphs(4,:)-max_phase),[],1);
            [~,I]=findpeaks(-abs(lap_thetaphs(mth,:)-cut_phase),'MinPeakHeight',-5);
            
            % 对应的位置
            IND = [];
            for i= 1:length(I)
                [~,IND(i)] = min(abs(lap_ts - lap_EEG_ts(I(i))));
            end
            
            thetacycle.ind{nl} = IND;% 全部的index
            thetacycle.ts{nl} = lap_ts(IND); % 全部的时间点
            thetacycle.EEG_ts{nl} = lap_EEG_ts(I); % 全部的EEG时间点
            for i = 1:length(I)-1
                thetacycle.cycle_ind{nl}(i,:) = [IND(i),IND(i+1)]; %每个theta周期开始和结束的index
                thetacycle.cycle_ts{nl}(i,:) = [lap_ts(IND(i)),lap_ts(IND(i+1))];%每个theta周期开始和结束的时间点
            end
        end
%     save(['CuttingPhaseCCW\data_AllLap_thetaphasecut_' num2str(cut_phase) '.mat'],'thetacycle','TS_running',...
%         'ph1','ph2','min_phase','name_tet','cut_phase','mth');
    end
end




%% function
function [TSlimit, speed, speed_ang]= TSfind_laplimit(data_angle,Ang_RewardLoc_ontrack,lap_speed,lap_pos,lap_ts, nsegment,nlap)
% input:
% data_angle ---- from file named Data_angle_ontrack.mat
% Ang_RewardLoc_ontrack ---- from file named like date_CT_tracking.mat
% output:
% TSlimit ---- the start and the end timgstamp of one lap
% speed ---- raw speed from data_angle
% speed_ang ---- raw angle speed from data_angle
speed = data_angle{nsegment}{nlap}(:,3);% cm/s
speed_ang = data_angle{nsegment}{nlap}(:,4);% rad/s
ind_ts_start = max(find(lap_speed'>=10 & lap_pos>=Ang_RewardLoc_ontrack(18)));
ind_ts_end = max(find(lap_speed'>=10 & lap_pos>=Ang_RewardLoc_ontrack(1)));
TSlimit = [lap_ts(ind_ts_start), lap_ts(ind_ts_end)];
end