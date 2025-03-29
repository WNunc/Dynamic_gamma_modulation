% do Bayesian decoding on all running laps in 4 types:
% pre-running, sample trials, test trials, and post-test trials
% modified from v4:
% this script should be run after doall_decoding_CT_v4.m having run
% 用spike数量筛选cell（10个每圈）
clear
disp('++++++start singlelap decoding++++++')
disp('exclude not fast gamma phase locking cell')
% directories_allData_v2
directories_allData_v0
%%
minspikes = 1;
mincells = 1;
dt = .04;% 40ms time bin
step = .005;% 5ms step length
peak0 = 1; % firing rate threshold to remove place cells
spk_threshold = 10;
TTList0 = 'TTList_dCA1_pyr.txt';
subfolder = 'Tseq\';
file_input1 = [subfolder 'Cells_allsegment_v1_vel_0.mat'];  % use to get all spikes
file_input2 = [subfolder 'Cells_allsegment_v1_vel_5.mat'];  % use to get all spikes
file_input3 = [subfolder 'Cells_singleLap_v2_vel_5.mat'];  % use itself as the decoder
file_input4 = [subfolder 'data_phaselocking_TSlap_vel0.mat'];  % ** phase locking feature trial TS
file_input5 = [subfolder 'data_phaselocking_TSlap_vel0_f2lap.mat'];  % ** phase locking feature trial TS first 2 lap
% file_input5 = [subfolder 'data_phaselocking_spkmin.mat'];  % phase locking feature running TS
fprintf(1,'all spike:\t%s\n%s\ndecoder:\t%s\nphase lock:\t%s\n%s\n',file_input1,file_input2,file_input3,file_input4,file_input5)
subfix = '-ontrack_dsfg';%
lockat = {'firstlap','alllap','f2lap'};

for L = 1:3 %1 = 第一圈相锁 ；2 = 全部圈相锁； 3 = 前两圈相锁
    for ns = 1:isession
        path_ns = path{ns};
        cd(path_ns);
        Ncell = getnumberofcells_cz_v1(TTList0);
        trackdata_ns = trackdata{ns};
        load(file_input4)
        load(file_input5)
        if Ncell < 1
            continue
        end
        disp(pwd)
        for nt = 1% this file for pre_running only
            switch nt
                case 1
                    laps=5;
                case 2
                    laps=8;
                case 3
                    laps=8;
                case 4
                    laps=6;
            end
            % Load Spikes data
            S1 = load(file_input1);  % used to get all spikes
            spikes = S1.spikes;
            S2 = load(file_input2);  % used to get segment ratemap
            Ratemap_seg = S2.Ratemap_seg{nt}; 
            % Remove the place cells whose peak firing rate<1Hz in all laps
            peak_all = max(Ratemap_seg); 
            ind = find(peak_all >= peak0);
            spikes = spikes(ind,:);
            % find cell on track as detecting
            X = cell2struct(S2.fieldProp_seg{nt}(ind),'placefield',1);
            load(trackdata_ns,'Ang_RewardLoc_ontrack')
            COM = [];
            for ncell  = 1:length(ind)
                COM(ncell) = X(ncell).placefield(1).x_COM;
            end
            
            cind_ot = find(COM>Ang_RewardLoc_ontrack(1) & COM<Ang_RewardLoc_ontrack(18));
            cellnumontrack = length(cind_ot);
            if L == 1
                fgamma_rayleighP = [plFeature_singlap{1,3}.rayleighP];
            elseif L == 2
                fgamma_rayleighP = [plFeature_alllap{1,3}.rayleighP];
            else
                fgamma_rayleighP = [plFeature_f2lap{1,3}.rayleighP];
            end
            fgamma_rayleighP_ontrack = fgamma_rayleighP(ind);
            fgamma_rayleighP_ontrack = fgamma_rayleighP_ontrack(cind_ot);
            
            fgPhase_singlap_ot = spikePhase_singlap(ind,3,:);
            fgPhase_singlap_ot = fgPhase_singlap_ot(cind_ot,1,:);
            fgPhase_singlap_ot = squeeze(fgPhase_singlap_ot);
            [spkNum_singlap,spkNum_ind] = spknumfinder(fgPhase_singlap_ot,spk_threshold);
            cind_ot_nonfg = find(fgamma_rayleighP_ontrack>=0.05 & spkNum_ind'==1); % ontrack的ind % 相锁神经元的ontrack的ind spkNum_ind'==1 |(spkNum_alllap'>50) 
            cind_ot_nonfg = cind_ot(cind_ot_nonfg)
            cellnumontrack_nonfg = length(cind_ot_nonfg);
            cind_ot_exnonfg = cind_ot(~ismember(cind_ot,cind_ot_nonfg))
            cellnumontrack_exnonfg = length(cind_ot_exnonfg);
            
            fprintf(1,['pre-running ',lockat{L},':\t%u cells ontrack and %u cell not phase locking in fgamma\n'],...
                cellnumontrack,cellnumontrack_nonfg)
            S3 = load(file_input3);  % used singlelap ratemap to make the decoder
            if cellnumontrack_fg == 0
                fprintf(1,'on fast gamma phase in this session (%s)\n', path_ns);
                continue
            end
            
            for nl = 1:laps
                Ratemap_SingleLaps_vel = S3.Ratemap_singlelap{nl,nt};
                Ratemap_SingleLaps_vel = Ratemap_SingleLaps_vel(:,ind);
                % use Ratemap_Ratemap_seg_vel as decoder
                [scores1,Ts_start_stop] = decodewindows_CT_v3(trackdata_ns,spikes(cind_ot_exnonfg,:),Ratemap_SingleLaps_vel(:,cind_ot_exnonfg),mincells,minspikes,[dt,step],nt,nl);%
                % save scores
                % scores = scores{1,nt};
                if ~isempty(scores1)
                    file_old=strcat('scores',num2str(nt)','-',num2str(nl),'-ontrack_v5.mat');
                    load ([subfolder file_old],'scores')
                    scores(:,1:9) = scores1{1,nt};
                end
                
                Ts_start_stop=Ts_start_stop{1,1};
                if nt ~= 1
                    A(:,1)=Ts_start_stop(:,1);
                    A(:,2)=Ts_start_stop(:,3);
                    Ts_reward=Ts_start_stop(:,2);
                    Ts_start_stop=A;
                    %                 file_output=strcat('scores',num2str(nt)','-',num2str(nl),subfix,'.mat');
                    file_output=strcat('scores',num2str(nt)','-',num2str(nl),subfix,'_TSlap_vel0',lockat{L},'v2.mat');% **
                    save([subfolder file_output],'scores','Ts_start_stop','Ts_reward','cind_ot','cind_ot_nonfg','cind_ot_exnonfg');
                else
                    %                 file_output=strcat('scores',num2str(nt)','-',num2str(nl),subfix,'.mat');
                    file_output=strcat('scores',num2str(nt)','-',num2str(nl),subfix,'_TSlap_vel0',lockat{L},'v2.mat');% **
                    save([subfolder file_output],'scores','Ts_start_stop','cind_ot','cind_ot_nonfg','cind_ot_exnonfg');
                end
            end
        end
    end
end


function [N,L] = spknumfinder(spike, threshold)
% 找出来每一圈数量不超过threshold的spike
% 输出每一圈spike数量和每一圈满足要求的cell的逻辑值
[i,j] = size(spike);
for I = 1:i
    for J = 1:j
        N(I,J) = length(spike{I,J});
        
        if N(I,J)>=threshold
            L(I,J) = 1;
        else
            L(I,J) = 0;
        end
    end
end
L = sum(L,2);
L = L == j;
end