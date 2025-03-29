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
% 全部
% fg相锁          匹配数量的非fg相锁   所有非相锁神经元
fg_avgFR_all = {};nfg_avgFR_all = {}; nfg_avgFR0_all = {};% 域内平均
fg_pkFR_all = {};nfg_pkFR_all = {}; nfg_pkFR0_all = {};% 峰值
fg_szFR_all = {};nfg_szFR_all = {}; nfg_szFR0_all = {};% 位置域大小
fg_meanFR_all = {};nfg_meanFR_all = {}; nfg_meanFR0_all = {};% 整个track的平均
outputName = 'data_ratemapOrder_exclude_cell';
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
    cind = [];% 存放解码排除掉的cell ind （ontrack的）
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
            saveflag = 0
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
            % 域内平均放电率
            fg_avgFR = fieldProp{cind_ot_fg(c)}(1).avgRate;
            nfg_avgFR = fieldProp{cind_ot_nonfg(c)}(1).avgRate;
            fg_avgFR_all{nsn,D}(c) = fg_avgFR;
            nfg_avgFR_all{nsn,D}(c) = nfg_avgFR;
            % 峰值放电率
            fg_pkFR = fieldProp{cind_ot_fg(c)}(1).peakRate;
            nfg_pkFR = fieldProp{cind_ot_nonfg(c)}(1).peakRate;
            fg_pkFR_all{nsn,D}(c) = fg_pkFR;
            nfg_pkFR_all{nsn,D}(c) = nfg_pkFR;
            % 位置域大小
            fg_SZ = fieldProp{cind_ot_fg(c)}(1).size;
            nfg_SZ = fieldProp{cind_ot_nonfg(c)}(1).size;
            fg_SZ_all{nsn,D}(c) = fg_SZ;
            nfg_SZ_all{nsn,D}(c) = nfg_SZ;
            % 质心位置
            fg_COM = fieldProp{cind_ot_fg(c)}(1).x_COM;
            nfg_COM = fieldProp{cind_ot_nonfg(c)}(1).x_COM;
            fg_COM_all(c) = fg_COM;
            nfg_COM_all(c) = nfg_COM;
            % 整体平均放电率
            fg_meanFR = mean(Ratemap_seg(:,cind_ot_fg(c)));
            nfg_meanFR = mean(Ratemap_seg(:,cind_ot_nonfg(c)));
            fg_meanFR_all{nsn,D}(c) = fg_meanFR;
            nfg_meanFR_all{nsn,D}(c) = nfg_meanFR;
        end
        COM_all_fg = [COM_all_fg,fg_COM_all];
        COM_all_nfg = [COM_all_nfg,nfg_COM_all];
        
        for c = 1:length(cind_ot_exfg)
            % 域内平均放电率
            nfg_avgFR0 = fieldProp{cind_ot_exfg(c)}(1).avgRate;
            nfg_avgFR0_all{nsn,D}(c) = nfg_avgFR0;
            % 峰值放电率
            nfg_pkFR0 = fieldProp{cind_ot_exfg(c)}(1).peakRate;
            nfg_pkFR0_all{nsn,D}(c) = nfg_pkFR0;
            % 位置域大小
            nfg_SZ0 = fieldProp{cind_ot_exfg(c)}(1).size;
            nfg_SZ0_all{nsn,D}(c) = nfg_SZ0;
            % 质心位置
            nfg_COM0 = fieldProp{cind_ot_exfg(c)}(1).x_COM;
            nfg_COM0_all(c) = nfg_COM0;
            % 整体平均放电率
            nfg_meanFR0 = mean(Ratemap_seg(:,cind_ot_exfg(c)));
            nfg_meanFR0_all{nsn,D}(c) = nfg_meanFR0;
        end
        COM0_all_nfg = [COM0_all_nfg,nfg_COM0_all];
    
    cind(D).fg = cind_ot_fg;
    cind(D).nfg = cind_ot_nonfg;
    cind(D).nfg0 = cind_ot_exfg;
    end
    
    fg_avgFR = fg_avgFR_all(nsn,:);
    nfg_avgFR = nfg_avgFR_all(nsn,:);
    nfg_avgFR0 = nfg_avgFR0_all(nsn,:);
    fg_pkFR = fg_pkFR_all(nsn,:);
    nfg_pkFR = nfg_pkFR_all(nsn,:);
    nfg_pkFR0 = nfg_pkFR0_all(nsn,:);
    fg_SZ = fg_SZ_all(nsn,:);
    nfg_SZ = nfg_SZ_all(nsn,:);
    nfg_SZ0 = nfg_SZ0_all(nsn,:);
    fg_meanFR = fg_meanFR_all(nsn,:);
    nfg_meanFR = nfg_meanFR_all(nsn,:);
    nfg_meanFR0 = nfg_meanFR0_all(nsn,:);
    
%     save([outFolder,'data_feature_cell.mat'],'cind',...
%         'fg_avgFR','nfg_avgFR','nfg_avgFR0',...
%         'fg_pkFR','nfg_pkFR','nfg_pkFR0',...
%         'fg_SZ','nfg_SZ','nfg_SZ0',...
%         'fg_meanFR','nfg_meanFR','nfg_meanFR0');
    
    fg_avgFR_mean(nsn) = mean([fg_avgFR_all{nsn,:}]);
    nfg_avgFR_mean(nsn) = mean([nfg_avgFR_all{nsn,:}]);
    nfg_avgFR0_mean(nsn) = mean([nfg_avgFR0_all{nsn,:}]);
    fg_pkFR_mean(nsn) = mean([fg_pkFR_all{nsn,:}]);
    nfg_pkFR_mean(nsn) = mean([nfg_pkFR_all{nsn,:}]);
    nfg_pkFR0_mean(nsn) = mean([nfg_pkFR0_all{nsn,:}]);
    fg_SZ_mean(nsn) = mean([fg_SZ_all{nsn,:}]);
    nfg_SZ_mean(nsn) = mean([nfg_SZ_all{nsn,:}]);
    nfg_SZ0_mean(nsn) = mean([nfg_SZ0_all{nsn,:}]);
    fg_meanFR_mean(nsn) = mean([fg_meanFR_all{nsn,:}]);
    nfg_meanFR_mean(nsn) = mean([nfg_meanFR_all{nsn,:}]);
    nfg_meanFR0_mean(nsn) = mean([nfg_meanFR0_all{nsn,:}]);
    
end
FG_avgFR = [fg_avgFR_all{:}];
NFG_avgFR = [nfg_avgFR_all{:}];
NFG_avgFR0 = [nfg_avgFR0_all{:}];
FG_pkFR = [fg_pkFR_all{:}];
NFG_pkFR = [nfg_pkFR_all{:}];
NFG_pkFR0 = [nfg_pkFR0_all{:}];
FG_SZ = [fg_SZ_all{:}];
NFG_SZ = [nfg_SZ_all{:}];
NFG_SZ0 = [nfg_SZ0_all{:}];
FG_meanFR = [fg_meanFR_all{:}];
NFG_meanFR = [nfg_meanFR_all{:}];
NFG_meanFR0 = [nfg_meanFR0_all{:}];
%%
[COM_fg_sort,ind_fg_sort] = sort(COM_all_fg);
[COM_nfg_sort,ind_nfg_sort] = sort(COM_all_nfg);
ratemap_fg_sort = ratemap_all_fg(:,ind_fg_sort);
ratemap_nfg_sort = ratemap_all_nfg(:,ind_nfg_sort);
Ncell_sort = length(COM_all_fg);
ffa = figure('Units','normalized','Position',[0.3 0.18 0.5 0.5]);%figure的参数设置
subplot(1,2,1)
imagesc(mapAxis,1:Ncell_sort,ratemap_fg_sort');axis xy
set(gca,'XTick',0:pi:2*pi);
set(gca,'XTickLabel',{'0','pi','2pi'});
ylim([1,Ncell_sort]);
set(gca,'YTick',[1,Ncell_sort]);
xlabel('Angle on the track (rad)')
ylabel('Cell ID')
set(gca,'fontsize',16);
axis xy
subplot(1,2,2)
imagesc(mapAxis,1:Ncell_sort,ratemap_nfg_sort');axis xy
set(gca,'XTick',0:pi:2*pi);
set(gca,'XTickLabel',{'0','pi','2pi'});
ylim([1,Ncell_sort]);
set(gca,'YTick',[1,Ncell_sort]);
xlabel('Angle on the track (rad)')
ylabel('Cell ID')
set(gca,'fontsize',16);
axis xy
colorbar('Position',...
    [0.933506951303143 0.109259259259259 0.0222222222222223 0.814814814814815]);
colormap(jet);

figure
h1 = histogram(COM_fg_sort,0:2*pi/18:2*pi);
h1.FaceColor = '#F24444';
h1.FaceAlpha = 1;
hold on
h2 = histogram(COM_nfg_sort,0:2*pi/18:2*pi,'Visible','off');
xbin = h2.BinEdges+h2.BinWidth*.5;
h2v = -h2.Values;
B = bar(xbin(1:end-1),h2v,1);
B.FaceColor = '#F2CA50';
hold off
xlim([0,2*pi])
ylim([-20,20])
xticks([0,pi,2*pi]);
xticklabels({'0','π','2π'});
xlabel('Angle on the track (rad)')
yticks([-20,-15,-10,-5,0,5,10,15,20]);
yticklabels({'20','15','10','5','0','5','10','15','20'});
ylabel('Cell count')
set(gca,'FontSize',14,'FontName','Arial')




[h,p,ks2stat] = kstest2(COM_fg_sort',COM_nfg_sort')
ffb = figure('Position', [361 271 464 536]);
H1 = cdfplot([0,COM_fg_sort]);
H1.LineWidth = 1.5;
H1.Color = '#F24444';

hold on
H2 = cdfplot([0,COM_nfg_sort]);
H2.LineWidth = 1.5;
H2.Color = '#F2CA50';
hold off
xlim([0,2*pi])
ylabel('Nomalizied cell number')
xlabel('Angle on the track(rad)')
grid off
title('CDF of cell')
legend({'fgcell','nfgcell'},'Location','northwest','Box','off')
set(gca,'FontSize',15)

cd('H:\neuralynx\gamma in sequence result\fieldProp\1115')
save('fieldProp_exCell_f2lap,mat','FG_avgFR','fg_avgFR_all','fg_avgFR_mean','NFG_avgFR','nfg_avgFR_all','nfg_avgFR_mean',...
    'FG_pkFR','fg_pkFR_all','fg_pkFR_mean','NFG_pkFR','nfg_pkFR_all','nfg_pkFR_mean',...
    'FG_SZ','fg_SZ_all','fg_SZ_mean','NFG_SZ','nfg_SZ_all','nfg_SZ_mean')
save('placefieldOrder_allCell_f2lap2.mat','ratemap_fg_sort','ratemap_nfg_sort','COM_all_fg','COM_all_nfg','mapAxis')
save('excell_KS-test_stat.mat','h','p','ks2stat')
saveas(ffa,'placefieldOrder_exCell_f2lap2.png')
saveas(ffa,'placefieldOrder_exCell_f2lap2','epsc')
saveas(ffa,'placefieldOrder_exCell_f2lap2.fig')

saveas(ffb,'CDF_exCell_f2lap2.png')
saveas(ffb,'CDF_exCell_f2lap2','epsc')
saveas(ffb,'FCF_exCell_f2lap2.fig')

%% 例子
rtMap{1, 2} = rtMap{1, 1}(:,ind_fg_sort)';
rtMap{2, 2} = rtMap{2, 1}(:,ind_nfg_sort)';
figure
subplot(1,2,1)
imagesc(rtMap{1, 2})
colormap('jet')
subplot(1,2,2)
imagesc(rtMap{2, 2})
colormap('jet')

figure('Position',[358 449 500 400]);
subplot(2,1,1)
plot(rtMap{1, 2}(22,:),'Color','#F24444','LineWidth',2)
ylabel('Fire rate (Hz)')
set(gca,'FontSize',12,'Color','none'); 
box off;axis tight;xticks([]);ylim([0,30])
legend({'nfgcell'},'Box','off')
subplot(2,1,2)
plot(rtMap{2, 2}(21,:),'Color','#F2CA50','LineWidth',2)
box off;axis tight;ylim([0,30])
xticks([1,90]);xticklabels({'0','2π'});xlabel('Position (rad)')
ylabel('Fire rate (Hz)')
set(gca,'FontSize',12,'Color','none'); 
legend({'fgcell'},'Box','off')
saveas(gcf,'example1.png')
saveas(gcf,'example1','epsc')

figure('Position',[358 449 500 400]);
subplot(2,1,1)
plot(rtMap{1, 2}(81,:),'Color','#F24444','LineWidth',2)
ylabel('Fire rate (Hz)')
set(gca,'FontSize',12,'Color','none'); 
box off;axis tight;xticks([]);ylim([0,25]);
subplot(2,1,2)
plot(rtMap{2, 2}(81,:),'Color','#F2CA50','LineWidth',2)
set(gca,'FontSize',12,'Color','none'); 
box off;axis tight;ylim([0,25]);
xticks([1,90]);xticklabels({'0','2π'});xlabel('Position (rad)')
ylabel('Fire rate (Hz)')
xticks([0,90]);xticklabels({0,'2π'});
saveas(gcf,'example2.png')
saveas(gcf,'example2','epsc')

figure('Position',[358 449 500 400]);
subplot(2,1,1)
plot(rtMap{1, 2}(40,:),'Color','#F24444','LineWidth',2)
ylabel('Fire rate (Hz)')
set(gca,'FontSize',12,'Color','none'); 
box off;axis tight;xticks([]);ylim([0,16]);
subplot(2,1,2)
plot(rtMap{2, 2}(40,:),'Color','#F2CA50','LineWidth',2)
set(gca,'FontSize',12,'Color','none'); 
box off;axis tight;ylim([0,16]);
xticks([1,90]);xticklabels({'0','2π'});xlabel('Position (rad)')
ylabel('Fire rate (Hz)')
saveas(gcf,'example3.png')
saveas(gcf,'example3','epsc')

figure('Position',[358 449 500 400]);
subplot(2,1,1)
plot(rtMap{1, 2}(55,:),'Color','#F24444','LineWidth',2)
ylabel('Fire rate (Hz)')
set(gca,'FontSize',12,'Color','none'); 
box off;axis tight;xticks([]);ylim([0,30]);
subplot(2,1,2)
plot(rtMap{2, 2}(58,:),'Color','#F2CA50','LineWidth',2)
set(gca,'FontSize',12,'Color','none'); 
box off;axis tight;ylim([0,30]);
xticks([1,90]);xticklabels({'0','2π'});xlabel('Position (rad)')
ylabel('Fire rate (Hz)')
saveas(gcf,'example4.png')
saveas(gcf,'example4','epsc')
