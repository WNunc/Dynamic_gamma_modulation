% 搭配doall_getThetaCycle_sgDomina_v2.m使用
% 但是根据spike重新计算相位
% 但是根据相位差分布画出差小于0的结果

clear
close all
n_std = 0;
directories_allData_v0_allgood
resuletFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
Dx = {'-cw','-ccw'};
case1 = '-ontrack';
ncells = 3;nspkikes = 5;
case2 = sprintf('_%ucell_%uspk',ncells, nspkikes);
midmod = '_cycmid';
nsn = 0;

peak0 = 1; % firing rate threshold to remove place cells

fileinput1 = ['sgamma_dominant_thetacyc_IND_std' num2str(n_std) '_v2.mat'];% sg dominant theta cycle ind
fileinput2 = [];% good theta cycle in each lap
% fileinput3 = 'data_phaselocking_TSlap_vel0.mat';% theta cycle
fileinput4 = 'Cells_allsegment_v1_vel_0.mat';% spike info
fileinput5 = 'Cells_allsegment_v1_vel_5.mat';% ratemap info
fileinput6 = [];% decoded info
fileinput7 = [];% all theta cycle in each lap
fileoutput1 = 'sgamma_phase_eachcell.mat';% sgamma
outputFolder = ['H:\neuralynx\sgamma result v3 std-' num2str(n_std) '\'];
mkdir(outputFolder)
nt = 1;nlap = 5;
for ns = 3:isession%[1,2,4:6,8,10:13,15,16]%
    path_ns = path{ns};
    CSCList0 = CSClist_CA1{ns};
    disp(path_ns)
    cd(path_ns)
    outFolder = [resuletFolder path_ns(13:end)];
    SPK = cell(1,nlap);SPKPhase = cell(1,nlap);
    CSCnum = [];
    CellIND = []; % 存放ontrack的cell的ind，将来能和fg相锁对应
    Cell_TheSG = [];% 存每个cell画图需要的变量
    LFP = cell(nlap,4);
    % merge the theta info on two direction
    for D = 1:2 % CW和CCW两个方向
        goodphase_ns = Seq_cutPhase{ns,D};%***
        case3 = num2str(goodphase_ns);% directories_allData_v0_allgood 俩方向相位一样
        % input folder
        inputFolder = [path_ns Directionfolder{D}];
        %% 导入spike
        % Load Spikes data
        S1 = load([inputFolder fileinput4]);  % used to get all spikes
        spikes = S1.spikes;
        TTList0 = S1.TTList0;
        [CSC0,cscnum]= creatCSCList(TTList0,CSCList0);
        cscnum = cscnum';
        S2 = load([inputFolder fileinput5]);  % used to get ratemap
        Ratemap_seg = S2.Ratemap_seg{nt};
        % Remove the place cells whose peak firing rate<1Hz in all laps
        peak_all = max(Ratemap_seg); % use S2 for detecting
        ind = find(peak_all >= peak0);
        % SPK{D} = spikes(ind,1); % spike TS
        SPK0 = spikes(ind,1);% spike TS
% %% 导入spikephase
% load([inputFolder fileinput3]);
% % SPKPhase{D} =spikePhase(ind,2); % sgamma phase
% SPKPhase0 = spikePhase(ind,2);% sgamma phase
        % 每个cell对应的csc
        CSCnum0 = cscnum(ind);
        CSCnum1 = [];
        %% 导入sg波形
        for nl = 1:nlap
            fileinput6 = strcat('scores',num2str(nt)','-',num2str(nl),'-ontrack_v5.mat');
            load([inputFolder fileinput6])
            LFP{nl,1} = [LFP{nl,1},scores{nl,10}];% TS of LFP
            LFP{nl,2} = [LFP{nl,2},scores{nl,17}];% bandpass sgamma
            LFP{nl,3} = [LFP{nl,3},scores{nl,18}];% bandpass theta
            LFP{nl,4} = [LFP{nl,4},scores{nl,20}];% bandpass sgamma phase
            % SPK0 = SPK0(cind_ot);SPKPhase0 = SPKPhase0(cind_ot);
            % 按lap找到对应的spike，和spike对应的相位
            SPK1 = {};SPKPhase1 = {};
            for ii = 1:length(SPK0)
                iiind = SPK0{ii,1}>Ts_start_stop(nl,1) & SPK0{ii,1}<Ts_start_stop(nl,2);
                SPK1{ii,1} = SPK0{ii,1}(iiind);
                [spikesEeg] = SpikeTStoEEGind(SPK1{ii,1}, LFP{nl,1});
                SPKPhase1{ii,1} = LFP{nl,4}(CSCnum0(ii),spikesEeg)';
            end
            
            SPK{1,nl} = [SPK{1,nl};SPK1(cind_ot)];
            SPKPhase{1,nl} = [SPKPhase{1,nl};SPKPhase1(cind_ot)];
        end
        CSCnum1 = CSCnum0(cind_ot);
        CSCnum = [CSCnum;CSCnum1];
        CellIND = [CellIND,cind_ot];
    end
    icell = length(SPK{1});
    %% 导入 good theta ind
    load([outFolder,fileinput1])
    fileinput7 = ['data_theta_seq_info_AllLap',case1,case2,case3,midmod,'_v5-both.mat'];
    load([outFolder,fileinput7])
    for nl = 1:nlap
        fileinput2 = ['data_theta_seq_info',case1,'_lap',num2str(nl),case2,case3,midmod,'_v5-both.mat'];
        load([outFolder,fileinput2],'ThetaGood')
        thcyc_sl = ind_th_sl{nl};
        nthc_sl = length(thcyc_sl);
        thcyc_al = ind_th_al{nl};% 5圈zscore sgpower
        nthc_al = length(thcyc_al);
        %% 单圈
        if ~isempty(thcyc_sl)
            for nth = 1:nthc_sl
                indth = thcyc_sl(nth);
                if indth > size(theta_INFO{nl},1)
                    continue
                end
                thcyc_onset = theta_INFO{nl}{indth,3}(1);
                thcyc_offset = theta_INFO{nl}{indth,3}(2);
                thseq = theta_INFO{nl}{indth,8};
                thseq_t = thcyc_onset:0.005:thcyc_offset;
                thseq_s = 0.5*2*pi/90:2*pi/90:2*pi-0.5*2*pi/90;
                thseq_cyc = [thcyc_onset,thcyc_offset];
                for nc = 1:icell
                    spkind = find(SPK{nl}{nc}>thcyc_onset & SPK{nl}{nc}<thcyc_offset);
                    if isempty(spkind)
                        continue
                    else
                        if ~ismember(theta_INFO{nl}{indth,1},[ThetaGood{:,1}])
                            c = 1;% 不在满足条件的sequence
                            figname = sprintf('B Lap-%u thetaCyc-%u cell-%u',nl, indth, nc);
                        else
                            c = 2;% 在满足条件的sequence
                            figname = sprintf('A Lap-%u thetaCyc-%u cell-%u',nl, indth, nc);
                        end
                        sz = ones(1,2*length(spkind))*88;
                        t = [SPK{nl}{nc}(spkind);SPK{nl}{nc}(spkind)]';
                        p = [SPKPhase{nl}{nc}(spkind);SPKPhase{nl}{nc}(spkind)+360];
                        itt = LFP{nl,1}>thcyc_onset&LFP{nl,1}<thcyc_offset;
                        SGTT = LFP{nl,1}(itt);
                        SGwave = LFP{nl,2}(CSCnum(nc),itt)+LFP{nl,3}(CSCnum(nc),itt);
                        SGphase = LFP{nl,4}(CSCnum(nc),itt);
                        t_sgcyc = findpeaks(SGphase);
                        t_sgcyc = t_sgcyc.loc;
                        t_sgcyc = SGTT(t_sgcyc(1:end-1));
                        [sgcycind] = get_sgphase(SGTT, SGphase,[thcyc_onset,thcyc_offset],SPK{nl}{nc}(spkind));
                        
                        % single cell SGphase in each lap
                        Cell_TheSG{nc,nl}{nth,1}= indth;
                        Cell_TheSG{nc,nl}{nth,2}= thseq;
                        Cell_TheSG{nc,nl}{nth,3}= thseq_s;
                        Cell_TheSG{nc,nl}{nth,4}= thseq_t;
                        Cell_TheSG{nc,nl}{nth,5}= thseq_cyc;
                        Cell_TheSG{nc,nl}{nth,6}= [t;p';sz]; % spike time/phase scatter size
                        Cell_TheSG{nc,nl}{nth,7}= [SGTT;SGwave;SGphase;...
                            LFP{nl,2}(CSCnum(nc),itt);LFP{nl,3}(CSCnum(nc),itt)]; 
                            % slow gamma time; sg+theta; sgphase;
                            % sg; theta
                        Cell_TheSG{nc,nl}{nth,8}= t_sgcyc; % sg cycle
                        Cell_TheSG{nc,nl}{nth,9}= sgcycind; 
                        Cell_TheSG{nc,nl}{nth,10}= c; % 是否是满足条件的theta
                        Cell_TheSG{nc,nl}{nth,11}= figname; % 是否是满足条件的theta
                        
                    end
                end
            end
        end
%         
%         if ~isempty(thcyc_al)
%             for nth = 1:nthc_al
%                 indth = thcyc_al(nth);
%                 thcyc_onset = ThetaGood{indth,3}(1);
%                 thcyc_offset = ThetaGood{indth,3}(2);
%                 for nc = 1:icell
%                     
%                     
%                     
%                     
%                 end
%             end
%         end
    end
    % 画图
    plot_cellPhs_in_sgDomThCyc_allLap
    save([outputFolder path_ns(13:end) fileoutput1],'Cell_TheSG','CellIND')
    close all
    
end



% function
%% 找到cell对应的CSC
function [CSC0,CSCnum]= creatCSCList(TTList0,CSCList0)
% open TTlist
    celllist = TTList0;
    fid=fopen(celllist);
    if (fid == -1)
        warning([ 'Could not open tfile ' celllist]);
    else
        % read the file names from the t-file list
        TT0 = ReadFileList(celllist);
        numCells0 = length(TT0);
        if numCells0==1 && max(TT0{1}==-1)
            % no cells in ttlist
            Ncell=0;
            error('no cells in ttlist')
        else
            Ncell=numCells0;
        end
    end
    
    % create the CSC file list
    CSC0 = cell(Ncell,1);
    for ii = 1:Ncell
        TTind = strfind(  TT0{ii}, '_' );
        CSC0{ii,1} = ['CSC' TT0{ii}(3:TTind-1) '.ncs'] ;
        CSC0{ii,2} = TT0{ii}(3:TTind-1);
        CSC0{ii,3} = str2double(CSC0{ii,2});
    end
    [~,CSCnum] = ismember([CSC0{:,3}],CSCList0);
end


%% 找到离spike发生最近的eeg时间点
function [spikesEeg] = SpikeTStoEEGind(spikeTS, eegTS)
spikesEeg = [];
%eliminate spikes that occur earlier than the first time stamp of the EEG
for k = 1:length(spikeTS)
    % Find closest eeg timestamp to the current spike timestamp
    tDiff = (eegTS-spikeTS(k)).^2;
    [~,eegTS_ind] = min(tDiff);
    spikesEeg = [spikesEeg,eegTS_ind];
end
spikesEeg(spikesEeg<1) = 1;
end



