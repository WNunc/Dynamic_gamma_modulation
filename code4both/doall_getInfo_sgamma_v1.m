% 搭配doall_getThetaCycle_sgDomina_v1.m使用

clear
close all

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

fileinput1 = 'sgamma_dominant_thetacyc_IND.mat';% good theta cycle ind
fileinput2 = [];% good theta cycle in each lap
fileinput3 = 'data_phaselocking_TSlap_vel0.mat';% theta cycle
fileinput4 = 'Cells_allsegment_v1_vel_0.mat';% spike info
fileinput5 = 'Cells_allsegment_v1_vel_5.mat';% ratemap info
fileinput6 = [];% decoded info

fileoutput1 = 'sgamma_dominant.mat';% sgamma
outputFolder = 'H:\neuralynx\sgamma result\';

nt = 1;nlap = 5;
for ns = 1:isession
    
    path_ns = path{ns};
    CSCList0 = CSClist_CA1{ns};
    disp(path_ns)
    cd(path_ns)
    outFolder = [resuletFolder path_ns(13:end)];
    SPK = cell(1,nlap);SPKPhase = cell(1,nlap);
    CSCnum = [];
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
        %% 导入spikephase
        load([inputFolder fileinput3]);
        % SPKPhase{D} =spikePhase(ind,2); % sgamma phase
        SPKPhase0 = spikePhase(ind,2);% sgamma phase
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
%             SPK0 = SPK0(cind_ot);SPKPhase0 = SPKPhase0(cind_ot);
            SPK1 = {};SPKPhase1 = {};
            for ii = 1:length(SPK0)
                iiind = SPK0{ii,1}>Ts_start_stop(nl,1) & SPK0{ii,1}<Ts_start_stop(nl,2);
                SPK1{ii,1} = SPK0{ii,1}(iiind);
                SPKPhase1{ii,1} = SPKPhase0{ii,1}(iiind);
            end
            
            SPK{1,nl} = [SPK{1,nl};SPK1(cind_ot)];
            SPKPhase{1,nl} = [SPKPhase{1,nl};SPKPhase1(cind_ot)];
        end
        CSCnum1 = CSCnum0(cind_ot);
        CSCnum = [CSCnum;CSCnum1];
        
    end
    icell = length(SPK{1});
    %% 导入 good theta ind
    load([outFolder,fileinput1])
    for nl = 1:nlap
        fileinput2 = ['data_theta_seq_info',case1,'_lap',num2str(nl),case2,case3,midmod,'_v5-both.mat'];
        load([outFolder,fileinput2],'ThetaGood')
        thcyc_sl = ind_th_sl{nl};
        nthc_sl = length(thcyc_sl);
        thcyc_al = ind_th_al{nl};
        nthc_al = length(thcyc_al);
        %% 单圈
        if ~isempty(thcyc_sl)
            for nth = 1:nthc_sl
                indth = thcyc_sl(nth);
                thcyc_onset = ThetaGood{indth,3}(1);
                thcyc_offset = ThetaGood{indth,3}(2);
                thseq = ThetaGood{indth,8};
                thseq_t = thcyc_onset:0.005:thcyc_offset;
                thseq_s = 0.5*2*pi/90:2*pi/90:2*pi-0.5*2*pi/90;
                for nc = 1:icell
                    spkind = find(SPK{nl}{nc}>thcyc_onset & SPK{nl}{nc}<thcyc_offset);
                    if isempty(spkind)
                        continue
                    else
                        figname = sprintf('Lap-%u thetaCyc-%u cell-%u',nl, indth, nc);
                        sz = ones(1,2*length(spkind))*88;
                        t = [SPK{nl}{nc}(spkind);SPK{nl}{nc}(spkind)]';
                        p = [SPKPhase{nl}{nc}(spkind),SPKPhase{nl}{nc}(spkind)+360];
                        itt = LFP{nl,1}>thcyc_onset&LFP{nl,1}<thcyc_offset;
                        SGTT = LFP{nl,1}(itt);
                        SGwave = LFP{nl,2}(CSCnum(nc),itt)+LFP{nl,3}(CSCnum(nc),itt);
                        SGphase = mod(LFP{nl,4}(CSCnum(nc),itt)+180,360);
                        t_sgcyc = findpeaks(SGphase);
                        t_sgcyc = t_sgcyc.loc;
                        t_sgcyc = SGTT(t_sgcyc(1:end-1));
                        
                        FFA = figure(1);
                        set(FFA,'Position',[771 180 395 629])
                        subplot(5,1,1:2)
                        imagesc(thseq_t,thseq_s,thseq);axis xy
                        ylabel('position(rad)')
                        colormap(jet)
                        xticklabels([])
                        xlim([thcyc_onset,thcyc_offset])
                        set(gca,'FontSize',9);
                        title(figname,'FontSize',14)
                        subplot(5,1,3:4)
                        scatter(t,p,sz, 'k|','LineWidth',1.5);%
                        hold on
                        stem(t_sgcyc,720*ones(1,length(t_sgcyc)),'r','LineWidth',1,'Marker','none','ShowBaseLine','off')
                        hold off
                        xlim([thcyc_onset,thcyc_offset])
                        
                        ylim([0 720])
                        yticks([0,360,720])
                        yticklabels([0,360,720])
                        xticklabels([])
                        ylabel('gamma_s phase(deg)')
                        set(gca,'FontSize',9);
                        
                        subplot(5,1,5)
                        plot(SGTT,SGwave)
                        maxy = get(gca,'YLim');
                        % set(gca,'ytick',[])
                        hold on
                        stem(t_sgcyc,maxy(2)*ones(1,length(t_sgcyc)),'r','LineWidth',1,'Marker','none','ShowBaseLine','off')
                        stem(t_sgcyc,maxy(1)*ones(1,length(t_sgcyc)),'r','LineWidth',1,'Marker','none','ShowBaseLine','off')
                        hold off
                        xlim([thcyc_onset,thcyc_offset])
                        xlabel('time(s)')
                        axis tight
                        ylabel('amp(μV)')
                        set(gca,'FontSize',9);
                        mkdir([outputFolder path_ns(13:end)])
                        saveas(FFA,[outputFolder path_ns(13:end) figname '.png'])
                        saveas(FFA,[outputFolder path_ns(13:end) figname '.epsc'])
                        clear FFA
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
%         
%         
    end
end

% 找到cell对应的CSC
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


% 按lap找到对应的spike

