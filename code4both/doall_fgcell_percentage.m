% crate by WN on 2022/08/28
% 统计相锁神经元和非相锁神经元占on track 神经元的比例
%%
clear
close all
directories_allData_v0_allgood
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
nsn = 0;
for ns = [1,2,4:6,8,10:13,15,16]
    nsn = nsn+1;
    path_ns = path{ns};
    disp(path_ns);cd(path_ns);
    outFolder = [resuletFolder path_ns(13:end)];
    inFolder = outFolder;
    SPK = {};
    
    cellNum_ot = 0;
    cellNum_ot_fg = 0;
    for D = 1:2
        goodphase_ns = Seq_cutPhase{ns,D};%***
        subfolder1 = Directionfolder{D};
        subfolder2 = Phasecutfolder{D};
        case3 = num2str(goodphase_ns);
        
        file_input1 = strcat(path_ns,subfolder1,...% 排除 fastgamma 神经元
            'scores',num2str(nseg),'-',num2str(1),case1{1},'_TSlap_vel0',lockat{L},'v2.mat');
        FGcell_ind = load(file_input1);
        cind_ot_fg = FGcell_ind.cind_ot_fg;
        cind_ot = FGcell_ind.cind_ot;
        cellNum_ot = cellNum_ot + length(cind_ot);
        cellNum_ot_fg = cellNum_ot_fg + length(cind_ot_fg);
        
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
    end
    ratio(nsn) = cellNum_ot_fg/cellNum_ot;
end
