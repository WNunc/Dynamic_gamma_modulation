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

outputName = 'data_ratemapOrder_exclude_cell';

ratemap_all = [];
COM_all = [];

for ns = [1,2,4:6,8,10:13,15,16]%1:isession
    path_ns = path{ns};
    disp(path_ns);cd(path_ns);
    outFolder = [resuletFolder path_ns(13:end)];
    inFolder = outFolder;
    rtMap = {};
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
        cind_ot = NFGcell_ind.cind_ot;
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
        peakFR = peak_all(ind);        
        Ratemap_seg = Ratemap_seg(:,ind);
        Ratemap_seg_norm = Ratemap_seg;
        Ratemap_seg_norm = Ratemap_seg_norm./peakFR;
%         rtMap{1,D} = Ratemap_seg(:,cind_ot_fg);% 被排除的相锁神经元
%         rtMap{2,D} = Ratemap_seg(:,cind_ot_nonfg);% 被排除的非相锁神经元
        ratemap_all = [ratemap_all,Ratemap_seg_norm(:,cind_ot)];

        fieldProp = S2.fieldProp_seg{nseg};
        fieldProp = fieldProp(ind);
        % 统计com
        cell_COM_all = [];
        for c = 1:length(cind_ot)
            COM = fieldProp{cind_ot(c)}(1).x_COM;
            cell_COM_all(c) = COM;
        end
        COM_all = [COM_all,cell_COM_all];

    end    
end

[COM_sort,ind_sort] = sort(COM_all);
%% 
ratemap_sort = ratemap_all(:,ind_sort);
Ncell_sort = length(COM_all);
ffa = figure('Units','normalized','Position',[0.3 0.11 0.25 0.5]);%figure的参数设置
imagesc(mapAxis,1:Ncell_sort,ratemap_sort');axis xy
set(gca,'XTick',0:pi:2*pi);
set(gca,'XTickLabel',{'0','pi','2pi'});
ylim([1,Ncell_sort]);
set(gca,'YTick',[1,Ncell_sort]);
xlabel('Angle on the track (rad)')
ylabel('Cell ID')
set(gca,'fontsize',14);
colorbar
colormap(jet);
cd('H:\neuralynx\gamma in sequence result\fieldProp\1115')
save('fieldProp_allCell_allsession,mat','COM_all','COM_sort','ind_sort','ratemap_all','mapAxis')
saveas(gcf,'placefieldOrder_allCell_f2lap2.png')
saveas(gcf,'placefieldOrder_allCell_f2lap2','epsc')

