% 计算fg相锁和非相锁的神经元的空间信息率
clear
close all
directories_allData_v0_allgood
resuletFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
case1 = {'-ontrack_exfg','-ontrack_dsfg'}; %exclude
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
fg_avgFR_all = {};nfg_avgFR_all = {}; % 域内平均
fg_pkFR_all = {};nfg_pkFR_all = {}; % 峰值
fg_szFR_all = {};nfg_szFR_all = {}; % 位置域大小
fg_meanFR_all = {};nfg_meanFR_all = {}; % 整个track的平均
outputName = 'data_spatialinformation_exclude_cell';
ratemap_all_fg = [];
ratemap_all_nfg = [];
COM_all_fg = [];
COM_all_nfg = [];
COM0_all_nfg = [];
nsn = 0;
rtMap = cell(3,1);
for ns = [1,2,4:6,8,10:13,15,16]%1:isession%
    nsn = nsn+1;
    path_ns = path{ns};
    disp(path_ns);cd(path_ns);
    outFolder = [resuletFolder path_ns(13:end)];
    cd(outFolder)
    inFolder = outFolder;
    
    for D = 1:2
        goodphase_ns = Seq_cutPhase{ns,D};%***
        subfolder1 = Directionfolder{D};
        subfolder2 = Phasecutfolder{D};
        case3 = num2str(goodphase_ns);
        
        file_input1 = strcat(path_ns,subfolder1,...% 排除 fastgamma 神经元
            'scores',num2str(nseg),'-',num2str(1),case1{1},'_TSlap_vel0',lockat{L},'v2.mat');
        FGcell_ind = load(file_input1);
        cind_ot_fg = FGcell_ind.cind_ot_fg;
        cind_ot_exfg = FGcell_ind.cind_ot_exfg;
        file_input2 = strcat(path_ns,subfolder1,...% 排除 非 fastgamma 神经元
            'scores',num2str(nseg),'-',num2str(1),case1{2},'_TSlap_vel0',lockat{L},'v3.mat');
        if ~exist(file_input2,'file')
            saveflag = 0;
            break
        end
        NFGcell_ind = load(file_input2);
        saveflag = 0;%%%
        cind_ot_nonfg = NFGcell_ind.cind_ot_nonfg;
        disp(strcat(num2str(FGcell_ind.cind_ot == NFGcell_ind.cind_ot)))%double check 一下之前做的解码对不对，都是 1 就对了
        
%         %%%%%%%%%%%%%%%%%%%%匹配放电率的cellind%%%%%%%%%%%%%%%%%%%%
%         load(['data_cellind_matchFR_D' num2str(D) '.mat'])
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Load Spikes data
        file_input3 = strcat(path_ns,subfolder1,'Cells_allsegment_v1_vel_0.mat');% load spike
        file_input4 = strcat(path_ns,subfolder1,'Cells_allsegment_v1_vel_5.mat');% load rate map
        
        % Load Position data
        file_input5 = strcat(path_ns,subfolder1,'Data_angle_ontrack.mat');
        load(file_input5);
        seg = 1;
        angle_ontrack = data_angle{seg};
        angleall = vertcat(angle_ontrack{:});
        % angleall0 = angleall(:,2);
        if D==1
            angleall0 = angleall(angleall(:,4)>0,2);
        else
            angleall0 = angleall(angleall(:,4)<0,2);
        end
        
        S1 = load(file_input3);  % used to get all spikes
        spikes = S1.spikes;
        S2 = load(file_input4);  % used to get segment ratemap
        Ratemap_seg = S2.Ratemap_seg{nseg};% 只有prerunning
        % Remove the place cells whose peak firing rate<1Hz in all laps
        peak_all = max(Ratemap_seg);
        ind = find(peak_all >= peak0);
        peakFR = peak_all(ind);
        Ratemap_seg = Ratemap_seg(:,ind);
        Ratemap_seg_norm = Ratemap_seg;
        Ratemap_seg_norm = Ratemap_seg_norm./peakFR;
        rtMap{1} = [rtMap{1},Ratemap_seg(:,cind_ot_fg)];% 被排除的相锁神经元
        rtMap{2} = [rtMap{2},Ratemap_seg(:,cind_ot_nonfg)];% 被排除的非相锁神经元
        rtMap{3} = [rtMap{3},Ratemap_seg(:,cind_ot_exfg)];% 被排除的非相锁神经元
        ratemap_all_fg = [ratemap_all_fg,Ratemap_seg_norm(:,cind_ot_fg)];
        ratemap_all_nfg = [ratemap_all_nfg,Ratemap_seg_norm(:,cind_ot_nonfg)];
        
        % 统计在全部lap中的放电率,以及com
        fieldProp = S2.fieldProp_seg{nseg};
        fieldProp = fieldProp(ind);
        fg_COM_all = [];nfg_COM_all = [];
        for c = 1:length(cind_ot_fg)
            % spatial information
            [fg_SI1,fg_SI2] = SpatialInfo(Ratemap_seg(:,cind_ot_fg(c)),angleall0,mapAxis);
            [nfg_SI1,nfg_SI2] = SpatialInfo(Ratemap_seg(:,cind_ot_nonfg(c)),angleall0,mapAxis);
            fg_SI_all{nsn,D}(c,1:2) = [fg_SI1,fg_SI2];
            nfg_SI_all{nsn,D}(c,1:2) = [nfg_SI1,nfg_SI2];
        end
        
        for c = 1:length(cind_ot_exfg)
            % spatial information
            [nfg_SI10,nfg_SI20] = SpatialInfo(Ratemap_seg(:,cind_ot_exfg(c)),angleall0,mapAxis);
            nfg_SI0_all{nsn,D}(c,1:2) = [nfg_SI10,nfg_SI20];
        end
        
    end
    fg_SI = fg_SI_all(nsn,:);
    nfg_SI = nfg_SI_all(nsn,:);
    nfg_SI0 = nfg_SI0_all(nsn,:);
    
    save([outFolder,'data_feature_cell.mat'],...
        'fg_SI','nfg_SI','nfg_SI0','-append');
    
    fg_SI_mean(nsn,1:2) = mean(vertcat(fg_SI_all{nsn,:}));
    nfg_SI_mean(nsn,1:2) = mean(vertcat(nfg_SI_all{nsn,:}));
    nfg_SI0_mean(nsn,1:2) = mean(vertcat(nfg_SI0_all{nsn,:}));
end

FG_SI = vertcat(fg_SI_all{:});
NFG_SI = vertcat(nfg_SI_all{:});
NFG_SI0 = vertcat(nfg_SI0_all{:});