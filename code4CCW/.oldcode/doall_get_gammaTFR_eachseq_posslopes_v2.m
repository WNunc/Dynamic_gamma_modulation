% Edit on 05/18/2019
% is a part of doall_get_gammaTFR_eachseq.m

 
% for positive slopes only
% ע�⣡ ˳ʱ�룬���ǿ���ֱ���������б�ʣ���ʱ���󸺵�б��
%�����Ҫ����ʱ����˳ʱ�룬������data on trackҪ�仯��ÿ��cell��ratemapҲӦ����Ӧ�仯��reward locationҲ����Ӧ�仯������sample-test�Ľ�����Ͳ���ͬһ���ˣ�
% ��¼ÿ��trail�У�reward֮ǰ�м���sequence��reward֮���м�����
% 20211214�汾���������ģ������������ĵ��ƶ����룬��ͼ������������������ʾ��λ��
% Nseq_everyday0{ns,nseg}(nl,1) = length(onset);%ÿһ��ÿһȦ��������sequence����
% Nseq_everyday1{ns,nseg}(nl,1) %ÿһ��ÿһȦ��TFR��nan��sequence����
clc
clear
sign_plot = 0; %�Ƿ�Ҫ��ͼ
file_cell_num = 'F:\grope data\cell\cell_cut_num.mat';
for icon=1:2
    if icon==1
        directories_allData_v2
        % load('F:\grope data\sequence\Ind_day.mat','Con_ind_day');
        % ind_day = Con_ind_day;
        data_folder ='F:\grope data\sequence\Control\';
        fig_dir = 'F:\grope data\sequence\Control\';
        file_input = 'BayesData_CircMap_v1_dt40ms_10ms_1cells_1spk_v5.mat';% use the place field in each session as different decoder
        fig_folder = 'Figure_Pxn_sequence_diff_dicoder_ds\'; % 4 rats, zscore gamma over all sessions, remove jumping-out points, only use eeg with max theta power
        load(file_cell_num,'Ncell_Con2');
        
    else
        directories_allData_v3
        % load('F:\grope data\sequence\Ind_day.mat','AD_ind_day');
        % ind_day = AD_ind_day;
        data_folder ='F:\grope data\sequence\AD\';
        fig_dir = 'F:\grope data\sequence\AD\';
        file_input = 'BayesData_CircMap_v1_dt40ms_10ms_1cells_1spk_v4.mat';% use the place field in each session as different decoder
        fig_folder = 'Figure_sequence_20211214\';
        load(file_cell_num,'Ncell_AD');
    end
    file_input_TFR = 'Data_TFR_gamma_maxtheta.mat'; % only use the tetrode with max theta power
    file_spike_input1 = 'Cells_ALLLaps_v2_vel_0.mat';  % used to get all spikes
    file_spike_input2 = 'Cells_eachsegment_wxl_2.mat';  % used to sort place field
    file_output1 = 'group_20211214.mat';
    file_output1 = strcat(data_folder,file_output1);
    a_illustration ='20211214�汾������б����ϵ�·������ͼ�еĺ���';
    % file_input_TFR = 'Data_TFR_gamma_SR.mat'; % only use the tetrode with max theta power
    TTList0 = 'TTList_dCA1_pyr.txt';
    dt = .04;
    step = .01;
    ncell_threshold = 25;% Only use data with more than 25 cells in a session
    
    % Default parameters
    maxJump_thr = 30;  % in position bins
    timeWin = 6;     % in time bin
    timeStep = 1;     % in time bin
    Distance_thr = 0; % in position bins
    jump_prop_thr = 0;
    
    Nseq = 0;
    Fs = 2000;
    Pxn_all = {};
    TFRz_all = {};
    data_info = [];
    para_all = [];
    
    
    if sign_plot ==1
        ffa = figure('Units','normalized','Position',[0 0 0.6 1]);
    end
    Diam_inner=100;
    Nseq_everyday0={};
    Nseq_everyday1={};
    
    for ns = 1:isession
        % for J =1:6
        %     ns=ind_day(J);
        path_ns = path{ns};
        cd(path_ns);
        disp(path_ns)
        
        csclist_ns = CSClist_CA1{ns};
        
        if exist(file_input,'file')>0 && ~isempty(csclist_ns)
            load(file_input,'scores');
            load(file_input_TFR);
            trackdata_ns = trackdata{ns};
            load(trackdata_ns,'Ind_rewardloc','Ind_rewardloc_sample','Ind_rewardloc_test','Ind_rewardloc_posttest',...
                'Sign_correct_sample','Sign_correct_test','Sign_correct_posttest',...
                'Ang_RewardLoc_ontrack','n_prerunning','ts_sample_reward','ts_test_reward','N_loc');
            Ts_sample_reward{1}=ts_sample_reward;
            Ts_test_reward{1}=ts_test_reward;
            S1 = load(file_spike_input1);  % used to get all spikes
            spikes = S1.spikes(:,1);
            Ncell = size(spikes,1);
            mapAxis = S1.mapAxis;
            S2 = load(file_spike_input2);  % used to get all spikes
            fieldProp_eachsegment = S2.fieldProp_eachsegment;
            COM_allcell = nan(Ncell,length(fieldProp_eachsegment));
            COM_cell_sort=[];
            ind_cell_sort=[];
            for nseg = 1:length(fieldProp_eachsegment)
                for nc = 1:Ncell
                    if ~isempty(fieldProp_eachsegment{nseg}{nc})
                        COM_allcell(nc,nseg) = fieldProp_eachsegment{nseg}{nc}(1,1).x_COM;  % x_COM of the 1st place field��ÿ��cell��λ�����ڹ���ϵ�bin��
                    end
                end
                [COM_cell_sort(:,nseg),ind_cell_sort(:,nseg)] = sort(COM_allcell(:,nseg));
            end
            
            %         Ind_rewardloc_all = [{nan(n_prerunning,1)},Ind_rewardloc_sample,Ind_rewardloc_test,Ind_rewardloc_posttest];
            %         Sign_correct_all = [{nan(n_prerunning,1)},Sign_correct_sample,Sign_correct_test,Sign_correct_posttest];
            Ind_rewardloc_all = [{nan(n_prerunning,1)},{nan(n_prerunning,1)},Ind_rewardloc_sample,Ind_rewardloc_test]; %ȥ����post test
            Sign_correct_all = [{nan(n_prerunning,1)},{nan(n_prerunning,1)},Sign_correct_sample,Sign_correct_test];
            Time_stop_all = [{nan(n_prerunning,1)},{nan(n_prerunning,1)},Ts_sample_reward,Ts_test_reward];
            %         a{1}=scores{1}(1:5,:);
            %         a{2}=scores{1}(6:10,:);a{2}(:,1)={1;2;3;4;5};
            %         a{3}=scores{2};
            %         a{4}=scores{3};
            %         scores=a; clear a  % pre ˳ʱ�룬��ʱ���
            
            %         %====================����cell����һ��������һ��======================%
            %         csclist_ns=gettetnumbers(TTList0);
            %         [e, eTSi, eTS, Fs_raw] = loadCSC_new_cz_allrecording(strcat(pwd,'\CSC',num2str(csclist_ns),'.ncs'));
            %
            %         eeg_t=eegfilt(e',Fs,0,12);%���� theta����
            %         eeg_th=eegfilt(eeg_t,Fs,8,0);
            %         eeg_t=eegfilt(e',Fs,0,55);%���� theta����
            %         eeg_sg=eegfilt(eeg_t,Fs,25,0);
            %         eeg_t=eegfilt(e',Fs,0,100);%���� theta����
            %         eeg_fg=eegfilt(eeg_t,Fs,55,0);
            %         clc;
            %         disp(path_ns)
            %         %=======================================================================%
            
    para=1/2;%parameter��ϵ��
        Ang_RewardLoc_ontrack_zone=nan(N_loc,2);%����һ��N_loc��2�����飬�������2����ʾÿ����������2����������߽���ұ߽磩
        for i=2:N_loc
            Ang_RewardLoc_ontrack_zone(i,1)=Ang_RewardLoc_ontrack(i)-(Ang_RewardLoc_ontrack(i)-Ang_RewardLoc_ontrack(i-1))*para;%��i���ǰһ������е㣬�ǵ�i�����ʼ�߽�
        end
        Ang_RewardLoc_ontrack_zone(1,1)=2*Ang_RewardLoc_ontrack(1)-Ang_RewardLoc_ontrack_zone(2,1);%��һ�������ʼ�߽�ĽǶ���ǰ����һ�㣬��ÿ����ķ�Χ������
        for i=1:N_loc-1
            Ang_RewardLoc_ontrack_zone(i,2)=Ang_RewardLoc_ontrack(i)+(Ang_RewardLoc_ontrack(i+1)-Ang_RewardLoc_ontrack(i))*para;%��i��ͺ�һ������е㣬�ǵ�i��Ľ�β�߽�
        end
        Ang_RewardLoc_ontrack_zone(N_loc,2)=2*Ang_RewardLoc_ontrack(N_loc)-Ang_RewardLoc_ontrack_zone(N_loc-1,2);%���һ����ı߽�Ƕ���������һ�㣬��ÿ����ķ�Χ������
       
            for nseg = 1:size(scores,2)
                score_nseg = scores{nseg};
                for nl = 1:size(score_nseg,1)
                    t_nl = score_nseg{nl,2};
                    pxn_nl = score_nseg{nl,3};
                    if isnan(t_nl(1))
                        continue
                    else
                        % Select continuous swipe ahead sequences
                        if nseg == 1|| nseg==3
                            [onset, offset, para,Centroid] = DetectSequenceEvents_pos(score_nseg{nl,3},mapAxis,...
                                maxJump_thr,timeWin,timeStep,Distance_thr,jump_prop_thr); %���ÿһ��sequence�Ŀ�ʼ�ͽ���bin
                        else
                            [onset, offset, para,Centroid] = DetectSequenceEvents_neg(score_nseg{nl,3},mapAxis,...
                                maxJump_thr,timeWin,timeStep,Distance_thr,jump_prop_thr);
                        end
                        Sign_correct_nl = Sign_correct_all{nseg}(nl);
                        Ind_rewardloc_nl = Ind_rewardloc_all{nseg}(nl);
                        Time_stop_nl = Time_stop_all{nseg}(nl)/10^6;
                        
                        Nseq0 = length(onset);
                        %                 Nseq_everyday0{J,nseg}(nl,1) =Nseq0;
                        Nseq_everyday0{ns,nseg}(nl,1) =Nseq0;
                        Nseq1 = 0;
                        if ~isnan(Sign_correct_nl)
                            Nseq_before_stop = 0;
                            Nseq_after_stop = 0;
                        else
                            Nseq_before_stop = nan;
                            Nseq_after_stop = nan;
                        end
                        % for each sequence
                        for i = 1:Nseq0
                            ind1 = onset(i); ind2 = offset(i);
                            
                            % get Pxn
                            pxn_nseq = score_nseg{nl,3}(:,ind1:ind2);
                            % get TFR
                            ind11 = find(abs(TFR_gamma{2,1}-score_nseg{nl,6}(ind1)) < 1/Fs);
                            if ~isempty(ind11)
                                [~,ind11] = min(abs(TFR_gamma{2,1}-score_nseg{nl,6}(ind1)));
                            end
                            ind22 = find(abs(TFR_gamma{2,1}-score_nseg{nl,6}(ind2)-(dt-step)) < 1/Fs);
                            if ~isempty(ind22)
                                [~,ind22] = min(abs(TFR_gamma{2,1}-score_nseg{nl,6}(ind2)-(dt-step)));
                            end
                            
                            if isempty(ind11) || isempty(ind22)
                                % if TFR has NaN in this sequence
                                continue
                            else
                                Nseq = Nseq+1;
                                Nseq1 = Nseq1+1;
                                
                                % get number of active cells and spikes
                                nspk = nan(Ncell,1);
                                for nc = 1:Ncell
                                    ind = find(spikes{nc,1} >= score_nseg{nl,6}(ind1) &...
                                        spikes{nc,1} < score_nseg{nl,6}(ind2)+dt-step);
                                    nspk(nc) = length(ind);
                                end
                                ind = find(nspk>0);
                                ncell_nseq = length(ind); %��һ��sequence�У��м���cell���˵磬�����˼���spike
                                nspk_nseq = sum(nspk);
                                
                                % get sequence fit
                                bin2use = find(~isnan(sum(pxn_nseq)));
                                % Cir_reg(���ʣ�Բ���90��λ��bin��ʱ�䣬ʹ�õ�sequence bin)
                                [~,calphase,~,~,slope] = ...
                                    Cir_reg(pxn_nseq,score_nseg{nl,5},score_nseg{nl,6}(ind1:ind2)-score_nseg{nl,6}(ind1),bin2use);
                                [~,posbin_ind] = histc(calphase,score_nseg{nl,5});
                                replaySC = replayScore_cir(pxn_nseq,posbin_ind'+1);
                                
                                % get Pxn for this sequence
                                [~, decodedbin] = max(pxn_nseq);
                                Pxn_all(Nseq,:) = [{score_nseg{nl,6}(ind1:ind2)},...
                                    {pxn_nseq},...
                                    {decodedbin},...
                                    {score_nseg{nl,5}(decodedbin)'},...
                                    {score_nseg{nl,4}(ind1:ind2)},...
                                    {score_nseg{nl,5}(score_nseg{nl,4}(ind1:ind2))},...
                                    {score_nseg{nl,8}(1,ind1:ind2)},...
                                    {score_nseg{nl,8}(4,ind1:ind2)},...
                                    {posbin_ind'+1}...
                                    {score_nseg{nl,5}(posbin_ind+1)}];
                                
                                % get TFR for each sequence
                                TFR_gamma{1,2} = squeeze(TFR_gamma{1,2});
                                TFR_gamma{2,2} = squeeze(TFR_gamma{2,2});
                                tfr_nseq = TFR_gamma{2,2}(:,ind11:ind22);
                                TFRz_all(Nseq,:) = [{TFR_gamma{2,1}(1,ind11:ind22)},...
                                    {tfr_nseq}];
                                
                                % get info for each sequence
                                seq_loc0 = {score_nseg{nl,5}(score_nseg{nl,4}(ind2))};%����seqʱ�����ﴦ���ĸ���
                               seq_loc = find(seq_loc0>Ang_RewardLoc_ontrack_zone(:,1)&seq_loc0<Ang_RewardLoc_ontrack_zone(:,2));
                                if ~isnan(Sign_correct_nl)
                                    data_info(Nseq,:) = [Ind_Rat(ns),...
                                        ns,...
                                        nseg,...
                                        nl,...
                                        Ind_rewardloc,...
                                        Sign_correct_nl,...
                                        Ind_rewardloc_nl-Ind_rewardloc,...
                                        Ang_RewardLoc_ontrack(Ind_rewardloc_nl),...
                                        Ang_RewardLoc_ontrack(Ind_rewardloc),...
                                        ncell_nseq,...
                                        nspk_nseq,...
                                        slope,...
                                        replaySC,...
                                        Time_stop_nl,...
                                        seq_loc];
                                else
                                    data_info(Nseq,:) = [Ind_Rat(ns),...
                                        ns,...
                                        nseg,...
                                        nl,...
                                        Ind_rewardloc,...
                                        Sign_correct_nl,...
                                        nan,...
                                        nan,...
                                        Ang_RewardLoc_ontrack(Ind_rewardloc),...
                                        ncell_nseq,...
                                        nspk_nseq,...
                                        slope,...
                                        replaySC,...
                                        nan];
                                end
                                
                                % in sample lap, use correct sign in its paired test lap
                                if nseg == 3
                                    Sign_correct_nl = Sign_correct_all{4}(nl);
                                    Ind_rewardloc_nl = Ind_rewardloc_all{4}(nl);
                                    data_info(Nseq,6) = Sign_correct_nl;
                                    if ~isnan(Sign_correct_nl)
                                        data_info(Nseq,7) = Ind_rewardloc_nl-Ind_rewardloc;
                                    else
                                        data_info(Nseq,7) = nan;
                                    end
                                end
                                
                                if  mean(Pxn_all{Nseq,1})< Time_stop_nl
                                    Nseq_before_stop = Nseq_before_stop+1;
                                elseif mean(Pxn_all{Nseq,1})> Time_stop_nl
                                    Nseq_after_stop = Nseq_after_stop+1;
                                end
                                
                                % get para for each sequence
                                para_all(Nseq,:) = para(i,:);
                                Centroid_all(Nseq,:) = Centroid(i,:);
                                %===================================== plot each sequence ========================================%
                                if sign_plot
                                    figure_eachseq
%                                     if Nseq==113
%                                         a=1
%                                     end;
                                end
                            end
                        end
                        %                 Nseq_everyday1{J,nseg}(nl,1) =Nseq1;
                        %                 Nseq_everyday1{J,nseg}(nl,2) =Nseq_before_stop;
                        %                 Nseq_everyday1{J,nseg}(nl,3) =Nseq_after_stop;
                        Nseq_everyday1{ns,nseg}(nl,1) =Nseq1;
                        Nseq_everyday1{ns,nseg}(nl,2) =Nseq_before_stop;
                        Nseq_everyday1{ns,nseg}(nl,3) =Nseq_after_stop;
                    end
                end
            end
            clear scores TFR_gamma score_nseg tfr_nseq S1
        end
        Nseq_everyday(ns,1) =Nseq; %һ��֮��һ��ɸѡ������sequence
    end
    % for J=1:6
    %     for nseg=1:4
    %             mean_nseg0(J,nseg) = mean(Nseq_everyday0{J,nseg});
    %             mean_nseg1(J,nseg) = mean(Nseq_everyday1{J,nseg}(:,1));
    %             mean_nseg_before_stop(J,nseg) = mean(Nseq_everyday1{J,nseg}(:,2));
    %             mean_nseg_after_stop(J,nseg) = mean(Nseq_everyday1{J,nseg}(:,3));
    %     end
    % end
    
    save(file_output1,'-v7.3');
end