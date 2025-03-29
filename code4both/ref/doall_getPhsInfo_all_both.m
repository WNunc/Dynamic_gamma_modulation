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
tcolor = {'black','red'};
peak0 = 1; % firing rate threshold to remove place cells

% fileinput1 = ['sgamma_dominant_thetacyc_IND_std' num2str(n_std) '_v2.mat'];% sg dominant theta cycle ind
% fileinput2 = [];% good theta cycle in each lap
% fileinput3 = 'data_phaselocking_TSlap_vel0.mat';% theta cycle
fileinput4 = 'Cells_allsegment_v1_vel_0.mat';% spike info
fileinput5 = 'Cells_allsegment_v1_vel_5.mat';% ratemap info
fileinput6 = [];% decoded info
%fileinput7 = [];% all theta cycle in each lap
fileoutput1 = 'data_info_allphase_both.mat';% phase info
outputFolder = ['H:\neuralynx\phase precession rseult\'];
mkdir(outputFolder)
nt = 1;nlap = 5;
for ns = 1:isession
    path_ns = path{ns};
    CSCList0 = CSClist_CA1{ns};
    disp(path_ns)
    cd(path_ns)
    outFolder = [resuletFolder path_ns(13:end)];
    SPK = cell(1,nlap);SPKPhase = cell(3,nlap);SPKPosition = cell(1,nlap);
    CSCnum = [];
    NC_ot = 0;
    M = [];% spike个数 cell x lap
    LFP = cell(nlap,4);
    
    thetaCSC = find(CSCList0 == csclist_CA1_mid{ns});
    
    SPK_all = cell(1,nlap);% all spike
    % final output spike in field info
    SPK_in = cell(1,nlap);
    
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
        fieldProp_seg = S2.fieldProp_seg;
        mapAxis = S2.mapAxis;
        % Remove the place cells whose peak firing rate<1Hz in all laps
        peak_all = max(Ratemap_seg); % use S2 for detecting
        ind = find(peak_all >= peak0);
        % SPK{D} = spikes(ind,1); % spike TS
        SPK0 = spikes(ind,1:2);% spike TS
        fieldProp_seg = fieldProp_seg{nt}(ind);
        
        % 每个cell对应的csc
        CSCnum0 = cscnum(ind);
        CSCnum1 = [];
        ml = [];%spike num on lap
        %% 导入theta、slowgamma、fastgamma相位
        for nl = 1:nlap
            fileinput6 = strcat('scores',num2str(nt)','-',num2str(nl),'-ontrack_v5.mat');
            load([inputFolder fileinput6])
            LFP{nl,1} = [LFP{nl,1},scores{nl,10}];% TS of LFP
            LFP{nl,2} = [LFP{nl,2},scores{nl,19}];% bandpass fgamma phase
            LFP{nl,3} = [LFP{nl,3},scores{nl,20}];% bandpass sgamma phase
            LFP{nl,4} = [LFP{nl,4},scores{nl,21}(thetaCSC,:)];% bandpass theta phase
            
            % 按lap找到对应的spike，和spike对应的三个相位
            
            SPKts = {};SPKphs = {};SPKpos = {};
            for ii = 1:size(SPK0,1)
                iiind = SPK0{ii,1}>Ts_start_stop(nl,1) & SPK0{ii,1}<Ts_start_stop(nl,2);
                SPKts{ii,1} = SPK0{ii,1}(iiind);
                
                [spikesEeg] = SpikeTStoEEGind(SPKts{ii,1}, LFP{nl,1});
                SPKphs{ii,1} = LFP{nl,2}(CSCnum0(ii),spikesEeg)';% fg phs
                SPKphs{ii,2} = LFP{nl,3}(CSCnum0(ii),spikesEeg)';% sg phs
                SPKphs{ii,3} = LFP{nl,4}(1,spikesEeg)';% th phs
                SPKpos{ii,1} = SPK0{ii,2}(iiind);
            end
            fieldP = fieldProp_seg(cind_ot);
            st = SPKts(cind_ot);
            sphs1 = SPKphs(cind_ot,1);%fg
            sphs2 = SPKphs(cind_ot,2);%sg
            sphs3 = SPKphs(cind_ot,3);%th
            spos = SPKpos(cind_ot);
            spike_all = [spos,sphs1,sphs2,sphs3,st];
            
            SPK_all{1,nl} = [SPK_all{1,nl};spike_all];
            
            % find spike in the place field
            NC_ot0 = length(cind_ot);
            NC_ot = NC_ot + NC_ot0;
            m = zeros(NC_ot0,1);
            for nc = 1:NC_ot0
                
                sT = st{nc};
                sPhs1 = sphs1{nc};
                sPhs2 = sphs2{nc};
                sPhs3 = sphs3{nc};
                sP = spos{nc};
                
                if isempty(fieldP{1,nc})
                    continue
                else
                    start=fieldP{nc}(1).startBin;
                    stop=fieldP{nc}(1).stopBin;
                end
                
                %                 switch D
                %                     case 1
                %                         ixx = sP>mapAxis(start) & sP<mapAxis(stop);
                %                         pos = sP(ixx);
                %                         nspk_in = length(pos);
                %                         phs1 = sPhs1(ixx);
                %                         phs2 = sPhs2(ixx);
                %                         phs3 = sPhs3(ixx);
                %                         ts = sT(ixx);
                %                     case 2
                %                         ixx = sP>mapAxis(stop) & sP<mapAxis(start);
                %                         pos = sP(ixx);
                %                         nspk_in = length(pos);
                %                         phs1 = sPhs1(ixx);
                %                         phs2= sPhs2(ixx);
                %                         phs3 = sPhs3(ixx);
                %                         ts = sT(ixx);
                %                 end
                ixx = sP>mapAxis(start) & sP<mapAxis(stop);
                pos = sP(ixx);
                nspk_in = length(pos);
                phs1 = sPhs1(ixx);
                phs2 = sPhs2(ixx);
                phs3 = sPhs3(ixx);
                ts = sT(ixx);
                
                m(nc) = nspk_in;
                spk_temp(nc,:) = {pos,phs1,phs2,phs3,ts};
                
            end
            SPK_in{1,nl} = [SPK_in{1,nl};spk_temp];
            ml = [ml,m];
            clear spk_temp
        end
        CSCnum1 = CSCnum0(cind_ot);
        CSCnum = [CSCnum;CSCnum1];
        M = [M;ml];
    end
    
    NC = length(SPK_in{1});
    ffn = fullfile(outputFolder,path_ns(13:end));
    mkdir(ffn)
    save([ffn,fileoutput1],'SPK_in','SPK_all','M','NC','CSCnum','TTList0')
    
    
end




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