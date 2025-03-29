% 分析用来做解码器的神经元是不是都有位置域
% 导入shuffle位置域之后的结果，ratemap和spatial information
% 导入真实的ratemap和spatial information
% 比较第一圈的结果

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
lockat = {'firstlap','alllap','f2lap'};L = 2;%用全部圈相锁的数据
peak0 = 1;
nlaps = 5;
nseg = 1;
lap = 1;
nsn = 0;
for ns = 1:isession
    path_ns = path{ns};
    cd(path_ns);
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
        
        file_input1 = strcat(path_ns,subfolder1,...% 排除 fastgamma 神经元
            'scores',num2str(nseg),'-',num2str(1),case1{1},'_TSlap_vel0',lockat{L},'v2.mat');
        FGcell_ind = load(file_input1);
        cind_ot = FGcell_ind.cind_ot;
        cind_ot_fg = FGcell_ind.cind_ot_fg;
        cind_ot_exfg = FGcell_ind.cind_ot_exfg;
        file_input2 = strcat(path_ns,subfolder1,...% 排除 非 fastgamma 神经元
            'scores',num2str(nseg),'-',num2str(1),case1{2},'_TSlap_vel0',lockat{L},'v3.mat');
        if ~exist(file_input2,'file')
            saveflag = 0;
            break
        end
        NFGcell_ind = load(file_input2);
        saveflag = 1;%%%
        cind_ot_nonfg = NFGcell_ind.cind_ot_nonfg;
        disp(strcat(num2str(FGcell_ind.cind_ot == NFGcell_ind.cind_ot)))%double check 一下之前做的解码对不对，都是 1 就对了
        % Load Spikes data
        
        file_input3 = strcat(path_ns,subfolder1,'Cells_allsegment_v1_vel_5.mat');% load rate map
        S1 = load(file_input3);  % used to get segment ratemap
        Ratemap_seg = S1.Ratemap_seg{nseg};% 只有prerunning
        % Remove the place cells whose peak firing rate<1Hz in all laps
        peak_all = max(Ratemap_seg);
        ind = find(peak_all >= peak0);
        
        file_input4 = [path_ns subfolder1 'Cells_singleLap_v2_vel_5.mat'];%实际的singlelap ratemap （decoder）
        file_input5 = [path_ns subfolder1 'Cells_singleLap_shuffle_vel_5_v3.mat'];%shuffled singlelap ratemap
        S2 = load(file_input4); 
        S3 = load(file_input5); 
        mapAxis = S3.mapAxis;
        Ratemap_singlelap = S2.Ratemap_singlelap;
        Ratemap_singlelap_sh = S3.ratemap_shuffled;    
        SInfo_temp = S3.Spatialinfo_shuffled;
        for lt = 1:nlaps
            Ratemap_singlelap{lt,nseg} = Ratemap_singlelap{lt,nseg}(:,ind);
            Ratemap_singlelap{lt,nseg} = Ratemap_singlelap{lt,nseg}(:,cind_ot);
            SInfo{lt,1} = SInfo_temp{lt}{101}(:,ind);
            SInfo{lt,1} = SInfo_temp{lt}{101}(:,cind_ot);
            for nshf = 1:100
                Ratemap_singlelap_sh{lt,nseg}{nshf} = Ratemap_singlelap_sh{lt,nseg}{nshf}(:,ind);
                Ratemap_singlelap_sh{lt,nseg}{nshf} = Ratemap_singlelap_sh{lt,nseg}{nshf}(:,cind_ot);
                SInfo_sh{lt,1}{nshf} = SInfo_temp{lt}{nshf}(:,ind);
                SInfo_sh{lt,1}{nshf} = SInfo_temp{lt}{nshf}(:,cind_ot);
            end
        end
        
        shuffle_mean_SI = [];
        shuffle_std_SI = [];
        shuffle_95_SI = [];
        for Ncell = 1:length(cind_ot)
            cm = [];
            cSI = [];
            for iii = 1:100
                cm(:,iii) = Ratemap_singlelap_sh{lap}{iii}(:,Ncell);
                cSI(:,iii) = SInfo_sh{lap}{iii}(:,Ncell);
            end
            
            shuffle_mean_rm = mean(cm');
            shuffle_std_rm = std(cm,[],2);
            CI95down = shuffle_mean_rm-1.96*shuffle_std_rm'/sqrt(100);
            CI95up = shuffle_mean_rm+1.96*shuffle_std_rm'/sqrt(100);
            
            shuffle_mean_SI(:,Ncell) = mean(cSI,2);
            shuffle_std_SI(:,Ncell) = std(cSI,[],2);
            shuffle_95_SI(:,Ncell) = prctile(cSI,95,2);
            
            % shuffle_sort = sort(cSI,2);
            % shuffle_95_SI(:,Ncell) = shuffle_sort(:,95);
            figure
            plot(mapAxis,cm,'Color',[0.85,0.85,0.85])
            title('100 x shuffled place field')
            xlabel('position(rad)')
            ylabel('firingrate(Hz)')
            xlim([0,2*pi])
            hold on
            plot(mapAxis,Ratemap_singlelap{lap,nseg}(:,Ncell),'Color',[1,0,0],'LineWidth',3)
            area = fill([mapAxis fliplr(mapAxis)],[CI95down fliplr(CI95up)],[0.15, 0.15, 0.15],'edgealpha', '0', 'facealpha', '0.8');
            %plot(mapAxis,shuffle_mean_rm,'w--','LineWidth',3)
            % cell2mat(Spatialinfo_shuffled(:))
        end
        
    end
end