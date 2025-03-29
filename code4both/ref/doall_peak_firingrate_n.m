% MNP=min peak firing rate
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
fileinput4 = 'Cells_singlelap_v2_vel_0.mat';% spike info
fileinput5 = 'Cells_allsegment_v1_vel_5.mat';% ratemap info
fileinput6 = [];% decoded info
%fileinput7 = [];% all theta cycle in each lap
fileoutput1 = 'data_PeakFRate_singlelap.mat';% peak firerate single lap
outputFolder = ['H:\neuralynx\phase precession rseult\'];
mkdir(outputFolder)
nt = 1;nlap = 5;
%%
for ns = 1:isession
    path_ns=path{ns};
    disp(path_ns);
    PFR=[];
    for nt = 1%1:3
        switch nt
            case 1
                laps=5;
            case 2
                laps=8;
            case 3
                laps=8;
        end
        %NC=S.Ncell;
        for D = 1:2 % CW和CCW两个方向
            goodphase_ns = Seq_cutPhase{nt,D};%***
            case3 = num2str(goodphase_ns);% directories_allData_v0_allgood 俩方向相位一样
            PFR0 = [];
            % input folder
            inputFolder = [path_ns Directionfolder{D}];
            %% 导入spike
            % Load Spikes data
            S1 = load([inputFolder fileinput4]);  % used to get all spikes
            spikes = S1.spikes;
            FProp=S1.fieldProp_singlelap;
            S2 = load([inputFolder fileinput5]);  % used to get ratemap
            Ratemap_seg = S2.Ratemap_seg{nt};
            % Remove the place cells whose peak firing rate<1Hz in all laps
            peak_all = max(Ratemap_seg); % use S2 for detecting
            ind = find(peak_all >= peak0);
            
            % load cind_ot
            fileinput6 = strcat('scores',num2str(nt)','-1','-ontrack_v5.mat');
            load([inputFolder fileinput6],'cind_ot')
            NC = length(cind_ot);
            
            for nc = 1:NC
                for nl = 1:laps
                    FP = FProp{nl,nt}(1,ind);
                    FP = FP(1,cind_ot);
                    if ~isempty(FP{1,nc})
                        PFR0{nc,nt}(1,nl)=FP{1,nc}(1).peakRate;
                    else
                        PFR0{nc,nt}(1,nl)=0;
                    end
                end
            end
            PFR = [PFR;PFR0];
        end
        
        %
        NP=[];
        for nt = 1%1:3
            for nc = 1: length(PFR)
                NP(nc,nt)=min(PFR{nc,nt});
            end
        end
        
        MNP=[];
        for nc = 1:length(NP)
            MNP(nc,1)=min(NP(nc,:));
        end
        
    end%
    file_output=fullfile(outputFolder,path_ns(13:end),fileoutput1);
    save(file_output,'PFR','MNP');
end