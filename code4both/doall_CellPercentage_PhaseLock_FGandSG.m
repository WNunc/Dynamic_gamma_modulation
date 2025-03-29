clear
close all
directories_allData_v0_allgood
resuletFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
case1 = {'-ontrack_exfg','-ontrack_exsg'}; %exclude
ncells = 3;nspkikes = 5;
case2 = sprintf('_%ucell_%uspk',ncells, nspkikes);
midmod = '_cycmid';
peak0 = 1;
nlaps = 5;
nseg = 1;

numbins = 90; % Number of bins
bin_ang = 2*pi/numbins;
mapAxis = bin_ang/2:bin_ang:(2*pi-bin_ang/2);

lockat = {'firstlap','alllap','f2lap'};L = 2;%用全部圈相锁的数据
outputName = 'data_PhsLockCell_FGandSG';

nsn = 0;
FG_num = 0;
SG_num = 0;
FGSG_num = 0;
for ns = [1,2,4:6,8,10:13,15,16]%1:isession%
    nsn = nsn+1;
    path_ns = path{ns};
    disp(path_ns);cd(path_ns);
    outFolder = [resuletFolder path_ns(13:end)];
    cd(outFolder)
    inFolder = outFolder;
    cind = [];% 存放解码排除掉的cell ind （ontrack的）
    for D = 1:2
        goodphase_ns = Seq_cutPhase{ns,D};%***
        subfolder1 = Directionfolder{D};
        subfolder2 = Phasecutfolder{D};
        case3 = num2str(goodphase_ns);
        
        file_input1 = strcat(path_ns,subfolder1,...% 排除 fastgamma 相锁神经元
            'scores',num2str(nseg),'-',num2str(1),case1{1},'_TSlap_vel0',lockat{L},'v2.mat');
        FGcell_ind = load(file_input1);
        cind_ot_fg = FGcell_ind.cind_ot_fg
        
        file_input2 = strcat(path_ns,subfolder1,...% 排除 slowgamma 相锁神经元
            'scores',num2str(nseg),'-',num2str(1),case1{2},'_TSlap_vel0',lockat{L},'v2.mat');
        SGcell_ind = load(file_input2);
        saveflag = 0;%%%
        cind_ot_sg = SGcell_ind.cind_ot_sg
        FG_num = FG_num+length(cind_ot_fg);
        SG_num = SG_num+length(cind_ot_sg);
        FGSG_num = FGSG_num+length(intersect(cind_ot_fg,cind_ot_sg));
    
    end
end

pie([368 7 11 102]/488,[0 0 1 1])


