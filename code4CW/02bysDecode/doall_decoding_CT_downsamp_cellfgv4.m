% do Bayesian decoding on all running laps in 4 types:
% pre-running, sample trials, test trials, and post-test trials
% modified from v4:
% this script should be run after doall_decoding_CT_v4.m having run
% 用平均放电率筛选cell（和fgcell的匹配）前提条件是运行:
% \..\code4both\doall_cell_order_2typecell.m 得到平均放电率
% 匹配相锁神经元和非相锁神经元数量
clear
disp('++++++start singlelap decoding++++++')
disp('exclude not fast gamma phase locking cell')
% directories_allData_v2
% directories_allData_v0
directories_allData_v0_allgood
%%
minspikes = 1;
mincells = 1;
dt = .04;% 40ms time bin
step = .005;% 5ms step length
peak0 = 1; % firing rate threshold to remove place cells
spk_threshold = 10;
TTList0 = 'TTList_dCA1_pyr.txt';
subfolder = 'Tseq\';
resuletFolder = 'H:\neuralynx\gamma in sequence result';

file_input1 = [subfolder 'Cells_allsegment_v1_vel_0.mat'];  % use to get all spikes
file_input2 = [subfolder 'Cells_allsegment_v1_vel_5.mat'];  % use to get all spikes
file_input3 = [subfolder 'Cells_singleLap_v2_vel_5.mat'];  % use itself as the decoder
file_input4 = [subfolder 'data_phaselocking_TSlap_vel0.mat'];  % ** phase locking feature trial TS
file_input5 = [subfolder 'data_phaselocking_TSlap_vel0_f2lap.mat'];  % ** phase locking feature trial TS first 2 lap
file_input6 = [];  % cell feature and cell ind

% file_input5 = [subfolder 'data_phaselocking_spkmin.mat'];  % phase locking feature running TS
fprintf(1,'all spike:\t%s\n%s\ndecoder:\t%s\nphase lock:\t%s\n%s\n',file_input1,file_input2,file_input3,file_input4,file_input5)
subfix = '-ontrack_dsfg';%
lockat = {'firstlap','alllap','f2lap'};
D = 1;%1 = CW

for L = 2 %1:3 %1 = 第一圈相锁 ；2 = 全部圈相锁； 3 = 前两圈相锁
    for ns = [1,2,4:6,8,10:13,15,16]%1:isession
        path_ns = path{ns};
        cd(path_ns);
        outfolder = [resuletFolder path_ns(13:end)];
        file_input6 = [outfolder,'data_feature_cell.mat'];
        load(file_input6)
        Ncell = getnumberofcells_cz_v1(TTList0);
        trackdata_ns = trackdata{ns};
        load(file_input4)
        load(file_input5)
        if Ncell < 1
            continue
        end
        disp(pwd)
        cind_ot = [];cind_ot_fg = [];
        cind_ot_nonfg = [];cind_ot_exnonfg = [];
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
            cind_ot_fg = cind(D).fg;
            cellnumontrack_fg = length(cind_ot_fg);
            cind_ot_exfg = cind(D).nfg0;
            fgFR = fg_meanFR{D};
            fgFRmean = mean(fg_meanFR{D});
            fgFRmax = max(fg_meanFR{D});
            fgFRmin = min(fg_meanFR{D});
            nfgFR = nfg_meanFR0{D};
            tempind = find(nfgFR>fgFRmin & nfgFR<fgFRmax);
            nfgFR0 = nfgFR(tempind);
            cind_ot_nonfg0 = cind_ot_exfg(tempind);
            cellnumontrack_nonfg = length(cind_ot_nonfg0);
            
%             cellnumontrack = length(cind_ot);
%             if L == 1
%                 fgamma_rayleighP = [plFeature_singlap{1,3}.rayleighP];
%             elseif L == 2
%                 fgamma_rayleighP = [plFeature_alllap{1,3}.rayleighP];
%             else
%                 fgamma_rayleighP = [plFeature_f2lap{1,3}.rayleighP];
%             end
%             fgamma_rayleighP_ontrack = fgamma_rayleighP(ind);
%             fgamma_rayleighP_ontrack = fgamma_rayleighP_ontrack(cind_ot);
            
            fgPhase_singlap_ot = spikePhase_singlap(ind,3,:);
            fgPhase_singlap_ot = fgPhase_singlap_ot(cind_ot,1,:);
            fgPhase_singlap_ot = squeeze(fgPhase_singlap_ot);
            [spkNum_singlap,spkNum_ind] = spknumfinder(fgPhase_singlap_ot,spk_threshold);
            spkNum_alllap = sum(spkNum_singlap,2);
%             cind_ot_fg = find(fgamma_rayleighP_ontrack<0.05 & spkNum_ind'==1);%find(fgamma_rayleighP_ontrack<0.05 & ) % ontrack的ind
%             cind_ot_fg = cind_ot(cind_ot_fg)
%             cellnumontrack_fg = length(cind_ot_fg);
%             
%             cind_ot_nonfg = find(fgamma_rayleighP_ontrack>=0.05 & spkNum_ind'==1); % ontrack的ind % 相锁神经元的ontrack的ind spkNum_ind'==1 |(spkNum_alllap'>50)
%             cind_ot_nonfg = cind_ot(cind_ot_nonfg)
%             cellnumontrack_nonfg = length(cind_ot_nonfg);
            
            if cellnumontrack_fg<cellnumontrack_nonfg% 如果fgcell数量少，那就从nonfg里随机取
                num = 1;
                tempind = sort(randperm(cellnumontrack_nonfg,cellnumontrack_fg),'ascend');
                cind_ot_nonfg = cind_ot_nonfg0(tempind);
                nfgFR = nfgFR0(tempind);
                while ttest(nfgFR,fgFR)
                    tempind = sort(randperm(cellnumontrack_nonfg,cellnumontrack_fg),'ascend');
                    cind_ot_nonfg = cind_ot_nonfg0(tempind);
                    nfgFR = nfgFR(tempind);
                    num = num+1;
                end
                cellnumontrack_nonfg = length(cind_ot_nonfg);
                cind_ot_exnonfg = cind_ot(~ismember(cind_ot,cind_ot_nonfg))
                cellnumontrack_exnonfg = length(cind_ot_exnonfg);
            elseif cellnumontrack_fg>cellnumontrack_nonfg % 从其余的放电率不太高的cell里挑几个
                Cdiff = cellnumontrack_fg-cellnumontrack_nonfg;
                tempind = find(~ismember(cind_ot_exfg,cind_ot_nonfg0));
                cellnumontrack_NONFG = length(tempind);
                ind_append = cind_ot_exfg(tempind);
                [B,I] = sort(nfgFR(tempind),'descend');
                ind_append = ind_append(I(1:Cdiff));
                cind_ot_nonfg = [cind_ot_nonfg0,ind_append];
                cellnumontrack_nonfg = length(cind_ot_nonfg);
                cind_ot_exnonfg = cind_ot(~ismember(cind_ot,cind_ot_nonfg))
                cellnumontrack_exnonfg = length(cind_ot_exnonfg);
                %                 rand('seed',ns)
                %                 tempind = sort(randperm(cellnumontrack_NONFG,Cdiff),'ascend');
                %                 cind_ot_NONFG = cind_ot_NONFG(tempind);
                %                 cind_ot_nonfg = [cind_ot_nonfg,cind_ot_NONFG]
                %                 cellnumontrack_nonfg = length(cind_ot_nonfg);
                %                 cind_ot_exnonfg = cind_ot(~ismember(cind_ot,cind_ot_nonfg))
                %                 cellnumontrack_exnonfg = length(cind_ot_exnonfg);
            else% 如果相等，直接用
                cind_ot_nonfg = cind_ot_nonfg0;
                cind_ot_exnonfg = cind_ot(~ismember(cind_ot,cind_ot_nonfg0))
                cellnumontrack_exnonfg = length(cind_ot_exnonfg);
            end
            
            spknumontrack_fg(ns) = sum(spkNum_alllap(ismember(cind_ot,cind_ot_fg)))
            spknumontrack_exfg(ns) = sum(spkNum_alllap(ismember(cind_ot,cind_ot_exfg)))
            spknumontrack_nonfg(ns) = sum(spkNum_alllap(ismember(cind_ot,cind_ot_nonfg)))
            
            fprintf(1,['pre-running ',lockat{L},':\t%u cells ontrack and %u cell phase locking in fgamma\n\texclude %u cell\n'],...
                cellnumontrack,cellnumontrack_fg,cellnumontrack_nonfg)
            S3 = load(file_input3);  % used singlelap ratemap to make the decoder
            if cellnumontrack_fg == 0
                fprintf(1,'on fast gamma phase in this session (%s)\n', path_ns);
                continue
            end
%             cd(outfolder)
%             datestr(now,'HH_mm_ss')
%             save('data_cellind_matchFR_D1.mat','cind_ot','cind_ot_fg','cind_ot_nonfg','cind_ot_exnonfg')
            %             for nl = 1:laps
            %                 Ratemap_SingleLaps_vel = S3.Ratemap_singlelap{nl,nt};
            %                 Ratemap_SingleLaps_vel = Ratemap_SingleLaps_vel(:,ind);
            %                 % use Ratemap_Ratemap_seg_vel as decoder
            %                 [scores1,Ts_start_stop] = decodewindows_CT_v3(trackdata_ns,spikes(cind_ot_exnonfg,:),Ratemap_SingleLaps_vel(:,cind_ot_exnonfg),mincells,minspikes,[dt,step],nt,nl);%
            %                 % save scores
            %                 % scores = scores{1,nt};
            %                 if ~isempty(scores1)
            %                     file_old=strcat('scores',num2str(nt)','-',num2str(nl),'-ontrack_v5.mat');
            %                     load ([subfolder file_old],'scores')
            %                     scores(:,1:9) = scores1{1,nt};
            %                 end
            %
            %                 Ts_start_stop=Ts_start_stop{1,1};
            %                 if nt ~= 1
            %                     A(:,1)=Ts_start_stop(:,1);
            %                     A(:,2)=Ts_start_stop(:,3);
            %                     Ts_reward=Ts_start_stop(:,2);
            %                     Ts_start_stop=A;
            %                     %                 file_output=strcat('scores',num2str(nt)','-',num2str(nl),subfix,'.mat');
            %                     file_output=strcat('scores',num2str(nt)','-',num2str(nl),subfix,'_TSlap_vel0',lockat{L},'v3.mat');% **
            %                     save([subfolder file_output],'scores','Ts_start_stop','Ts_reward','cind_ot','cind_ot_nonfg','cind_ot_exnonfg');
            %                 else
            %                     %                 file_output=strcat('scores',num2str(nt)','-',num2str(nl),subfix,'.mat');
            %                     file_output=strcat('scores',num2str(nt)','-',num2str(nl),subfix,'_TSlap_vel0',lockat{L},'v3.mat');% **
            %                     save([subfolder file_output],'scores','Ts_start_stop','cind_ot','cind_ot_nonfg','cind_ot_exnonfg');
            %                 end
            %             end
        end
    end
end


function [N,L] = spknumfinder(spike, threshold)
% 找出来每一圈数量超过数量threshold的spike
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