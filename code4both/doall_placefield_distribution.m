%% 神经元位置域分布


clear
close all
directories_allData_v0_up
resuletFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
case1 = '-ontrack';
ncells = 3;nspkikes = 5;
case2 = sprintf('_%ucell_%uspk',ncells, nspkikes);
midmod = '_cycmid';
TPsweep = 7;
SPsweep = 10;
peak0 = 1;
nseg = 1;
rMap_fg_alls = [];rMap_nfg_alls = [];
for ns = 1:isession
    path_ns = path{ns};
    cd(path_ns)
    disp(path_ns);cd(path_ns);
    trackdata_ns = trackdata{ns};
    load(trackdata_ns,'Ang_RewardLoc_ontrack');
    
    rMap_fg_alll = [];rMap_nfg_alll = [];
    for D = 1:2
        subfolder1 = Directionfolder{D};
        fileinput1 = [subfolder1 'Cells_allsegment_v1_vel_0.mat'];
        S1 = load(fileinput1);
        fileinput2 = [subfolder1 'Cells_allsegment_v1_vel_5.mat'];
        S2 = load(fileinput2);
        fileinput3 = [subfolder1,'scores1-1-ontrack_exfg_TSlap_vel0alllap.mat'];
        load(fileinput3,'cind_ot','cind_ot_fg')
        fileinput4 = [subfolder1,'scores1-1-ontrack_dsfg_TSlap_vel0alllap.mat'];
        load(fileinput4,'cind_dsamp')
        cind_ot_ind_nfg = ~ismember(cind_ot,cind_dsamp);
        cind_ot_nfg = cind_ot(cind_ot_ind_nfg);
        
        spikes = S1.spikes; % use S1 spike
        ratemaps = S2.Ratemap_seg{nseg};
        peak_all = max(ratemaps); % use S2 ratemap
        ind = find(peak_all >= peak0);
        % 处理一下得到解码用到的spike 和 ratemap
        Spikes = spikes(ind,:);% 这个
        %         Spikes = Spikes(cind_ot,:);
        Ratemaps = ratemaps(:,ind);% 这个
        %         Ratemaps = Ratemaps(:,cind_ot);
        
        % 相锁神经元的ratemap
        Ratemaps_fg = Ratemaps(:,cind_ot_fg);
        % 非相锁神经元（随机选出来的神经元）的ratemap
        Ratemaps_nfg = Ratemaps(:,cind_ot_nfg);
        % 每个session 中 all lap 的ratemap
        rMap_fg_alll = [rMap_fg_alll,Ratemaps_fg];
        rMap_nfg_alll = [rMap_nfg_alll,Ratemaps_nfg];
        
        
    end
    % 所有session的ratemap
    rMap_fg_alls = [rMap_fg_alls,rMap_fg_alll];
    rMap_nfg_alls = [rMap_nfg_alls,rMap_nfg_alll];
    
    rMap_fg_alll_mean = mean(rMap_fg_alll,2);
    rMap_nfg_alll_mean = mean(rMap_nfg_alll,2);
    rMap_fg_alll_sem = std(rMap_fg_alll,0,2)./size(rMap_fg_alll,2);
    rMap_nfg_alll_sem = std(rMap_nfg_alll,0,2)./size(rMap_nfg_alll,2);
    
    mapAxis = S1.mapAxis;
    
    [~,diff1] = min(abs(mapAxis - Ang_RewardLoc_ontrack(1)));
    [~,diff6] = min(abs(mapAxis - Ang_RewardLoc_ontrack(6)));
    [~,diff12] = min(abs(mapAxis - Ang_RewardLoc_ontrack(12)));
    [~,diff18] = min(abs(mapAxis - Ang_RewardLoc_ontrack(18)));
    pos_part = [diff6, diff12];
    
    figure%每个神经元的ratemap
    subplot(2,2,1)
    plotRatemapDist(rMap_fg_alll,pos_part);
    xlim([diff1,diff18])
    title('fgcell')
    subplot(2,2,2)
    plotRatemapDist(rMap_nfg_alll,pos_part);
    xlim([diff1,diff18])
    title('non-fgcell')
    subplot(2,2,3)
    plotRatemapDist(rMap_fg_alll_mean,pos_part);
    xlim([diff1,diff18])
    title('mean')
    subplot(2,2,4)
    plotRatemapDist(rMap_nfg_alll_mean,pos_part);
    xlim([diff1,diff18])
    title('mean')
    
    %     field_alllap 变量 画出来分布
    %     分布1.直接画线图
    %     2.画平均后的bar图
    %     在三个位置上画竖线
    
end
rMap_fg_alls_mean = mean(rMap_fg_alls,2);
rMap_nfg_alls_mean = mean(rMap_nfg_alls,2);
rMap_fg_alls_sem = std(rMap_fg_alls,0,2)/sqrt(size(rMap_fg_alls,2));
rMap_nfg_alls_sem = std(rMap_nfg_alls,0,2)/sqrt(size(rMap_nfg_alls,2));
% field_alllap 变量 画出来分布
% 分布1.直接画线图
% 2.画平均后的bar图
% 在三个位置上画竖线
figure%每个神经元的ratemap
subplot(3,1,1)
plotRatemapDist(rMap_fg_alls,pos_part);
xlim([diff1,diff18])
title('fgcell')
subplot(3,1,2)
plotRatemapDist(rMap_nfg_alls,pos_part);
xlim([diff1,diff18])
title('non-fgcell')
subplot(3,1,3)
plotRatemapDist(rMap_fg_alls_mean,pos_part);
errorbar(rMap_fg_alls_mean,rMap_fg_alls_mean)
hold on
errorbar(rMap_nfg_alls_mean,rMap_nfg_alls_mean)
maxy = get(gca,'YLim');
x = pos_part;
y = ones(1,size(pos_part,2)).*maxy(2);
stem(x,y,'r','Marker','none','LineWidth',2);
hold off
axis tight
xlim([diff1,diff18])
title('mean')
% subplot(2,2,4)
% plotRatemapDist(rMap_nfg_alls_mean,pos_part);
% 
% xlim([diff1,diff18])
% title('mean')
ratio = rMap_fg_alls_mean./rMap_nfg_alls_mean
find(ratio>1.5)

function plotRatemapDist(ratemap,pos)
plot(ratemap)
maxy = get(gca,'YLim');
x = pos;
y = ones(1,size(pos,2)).*maxy(2);
hold on
stem(x,y,'r','Marker','none','LineWidth',2);
hold off
axis tight
end