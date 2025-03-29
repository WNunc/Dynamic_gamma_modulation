% 找到每个theta cycle的线速度
% 统计各圈的速度变化

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
nlap = 5;
outputFolder = 'H:\neuralynx\sgamma result  v2\';
vmean_all = [];vmean_good_all = [];
for ns = [1,2,4:6,8,10:13,15,16]%1:isession%
    nsn = nsn + 1;
    path_ns = path{ns};
    disp(path_ns)
    cd(path_ns)
    outFolder = [resuletFolder path_ns(13:end)];
    D = 1;
    goodphase_ns = Seq_cutPhase{ns,D};%***
    case3 = num2str(goodphase_ns);% directories_allData_v0_allgood 俩方向相位一样
    fileinput1 = ['data_theta_seq_info_AllLap',case1,case2,case3,midmod,'_v5-both.mat'];
    load([outFolder,fileinput1])
    % 在track上的速度
    v_mean = [];
    for nl = 1:nlap
        vmean(nl) = mean([theta_INFO{nl}{:,11}],'omitnan');
        %满足要求的sequence的速度
        fileinput2 = ['data_theta_seq_info',case1,'_lap',num2str(nl),case2,case3,midmod,'_v5-both.mat'];
        load([outFolder,fileinput2],'ThetaGood')
        vmean_good(nl) = mean([ThetaGood{:,11}],'omitnan');
        
    end
    
    
    vmean_all(:,nsn) = vmean;
    vmean_good_all(:,nsn) = vmean_good;
    
end


figure
vmean_all_mean = mean(vmean_all,2);
vmean_all_std = std(vmean_all,0,2);
errorbar(1:nlap,vmean_all_mean,vmean_all_std/sqrt(nsn),'k', 'LineWidth',2)
xlim([0.5,5.5]);ylim([0,40]);
xlabel('Lap Num');ylabel('speed(cm/s)')
title('on-track speed')
axis square
set(gca,'FontSize',16)

figure
vmean_good_all_mean = mean(vmean_good_all,2);
vmean_good_all_std = std(vmean_good_all,0,2);
errorbar(1:nlap,vmean_good_all_mean,vmean_good_all_std/sqrt(nsn),'k', 'LineWidth',2)
xlim([0.5,5.5]);ylim([0,40]);
xlabel('Lap Num');ylabel('speed(cm/s)')
title('goodcyc speed')
axis square
set(gca,'FontSize',16)

