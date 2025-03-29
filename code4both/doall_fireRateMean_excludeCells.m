% crate by WN on 2022/08/08
% 统计相锁神经元和非相锁神经元之间的平均放电率的差异
%%
clear
close all
directories_allData_v0
resuletFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
case1 = {'-ontrack_exfg','-ontrack_dsfg'}; %exclude non-fgamma
ncells = 3;nspkikes = 5;
case2 = sprintf('_%ucell_%uspk',ncells, nspkikes);
midmod = '_cycmid';

peak0 = 1;
gap =1;% tov detect SWRs
sampfreq = 2000;% EEG sapmle rate
nlaps = 5;
nseg = 1;
lockat = {'firstlap','alllap','f2lap'};L = 2;%用全部圈相锁的数据
fg_avgFR_all = {};nfg_avgFR_all = {};
outputName = 'data_fireRate_exclude_cell.mat';
for ns = [2;4;5;6;8;10;11;16;20;22;25]'%1:isession %[1,2,4,10,21,24]%[1:11,16,20:22,25,26,28:30]%
    path_ns = path{ns};
    disp(path_ns);cd(path_ns);
    outFolder = [resuletFolder path_ns(13:end)];
    inFolder = outFolder;
    SPK = {};
    for D = 1:2
        goodphase_ns = Seq_cutPhase{ns,D};%***
        subfolder1 = Directionfolder{D};
        subfolder2 = Phasecutfolder{D};
        case3 = num2str(goodphase_ns);
        
        file_input1 = strcat(path_ns,subfolder1,...% 排除 fastgamma 神经元
            'scores',num2str(nseg),'-',num2str(1),case1{1},'_TSlap_vel0',lockat{L},'v2.mat');
        FGcell_ind = load(file_input1);
        cind_ot_fg = FGcell_ind.cind_ot_fg;
        file_input2 = strcat(path_ns,subfolder1,...% 排除 非 fastgamma 神经元
            'scores',num2str(nseg),'-',num2str(1),case1{2},'_TSlap_vel0',lockat{L},'v3.mat');
        if ~exist(file_input2,'file')
            saveflag = 0
            break
        end
        NFGcell_ind = load(file_input2);
        saveflag = 0;%%%
        cind_ot_nonfg = NFGcell_ind.cind_ot_nonfg;
        disp(strcat(num2str(FGcell_ind.cind_ot == NFGcell_ind.cind_ot)))%double check 一下之前做的解码对不对，都是 1 就对了
        % Load Spikes data
        file_input3 = strcat(path_ns,subfolder1,'Cells_allsegment_v1_vel_0.mat');% load spike
        file_input4 = strcat(path_ns,subfolder1,'Cells_allsegment_v1_vel_5.mat');% load rate map
        
        S1 = load(file_input3);  % used to get all spikes
        spikes = S1.spikes;
        S2 = load(file_input4);  % used to get segment ratemap
        Ratemap_seg = S2.Ratemap_seg{nseg};% 只有prerunning
        % Remove the place cells whose peak firing rate<1Hz in all laps
        peak_all = max(Ratemap_seg);
        ind = find(peak_all >= peak0);
        spikes = spikes(ind,:);
        SPK{1,D} = spikes(cind_ot_fg,:);% 被排除的相锁神经元
        SPK{2,D} = spikes(cind_ot_nonfg,:);% 被排除的非相锁神经元
        % 统计平均放电率
        fieldProp = S2.fieldProp_seg{nseg};
        fieldProp = fieldProp(ind);
        
        for c = 1:length(cind_ot_fg)
            fg_avgFR = fieldProp{cind_ot_fg(c)}(1).avgRate;
            nfg_avgFR = fieldProp{cind_ot_nonfg(c)}(1).avgRate;
            fg_avgFR_all{ns,D}(c) = fg_avgFR;
            nfg_avgFR_all{ns,D}(c) = nfg_avgFR;
        end
        fg_avgFR_mean(ns) = mean([fg_avgFR_all{ns,:}]);
        nfg_avgFR_mean(ns) = mean([nfg_avgFR_all{ns,:}]);
        
    end
    if saveflag~=0
        save([outFolder,outputName],'SPK','fg_avgFR','nfg_avgFR')
    end
end
save([resuletFolder,'\fireRate-result\Data_avgFR_all.mat'],...
    'fg_avgFR_all','nfg_avgFR_all','fg_avgFR_mean','nfg_avgFR_mean')
fg_avgFR_mean = fg_avgFR_mean(fg_avgFR_mean~=0)';
nfg_avgFR_mean = nfg_avgFR_mean(nfg_avgFR_mean~=0)';
avgFR_mean = [fg_avgFR_mean,nfg_avgFR_mean];

