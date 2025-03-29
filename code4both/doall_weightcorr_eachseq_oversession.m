% 计算每个session中不分圈的全部sequence在3种解码条件下的WC值
% 按session统计WC

clear
close all
% directories_allData_v0
directories_allData_v0_allgood
resuletFolder = 'H:\neuralynx\gamma in sequence result';

subfix = '-ontrack'; % v5
subfix1 = '-ontrack_exsg'; %exclude sgamma
subfix2 = '-ontrack_exfg'; %exclude fgamma
subfix3 = '-ontrack_dssg'; %down sample sgamma
subfix4 = '-ontrack_dsfg'; %down sample fgamma
ncells = 3;nspkikes = 5;
ncases = sprintf('_%ucell_%uspk', ncells, nspkikes);
midmod = '_cycmid';
lockat = {'firstlap','alllap','f2lap'};


for L = 2
    nsn = 0;
    wc1_all = cell(1,5);wc2_all = cell(1,5);
    WC1_all = [];WC2_all = [];
    for ns = [1,2,4:6,8,10:13,15,16]%1:isession %[1;2;4;5;6;8;10;11;16;20;22;25]'% [1:11,16,20:22,25,26,28,29,30]%[1:4,6,7]% [1;2;5;6;8;10;11;13;20;22;23;24;27;28]'%[1:8,10:15,22:31]%[1:9,11,13:16]%
        
        path_ns = path{ns};
        cd(path_ns)
        disp(path_ns)
        goodphase_ns = Seq_cutPhase{ns,2};
        case3 = num2str(goodphase_ns);
        
        inFolder = [resuletFolder path_ns(13:end)];
        outFolder = inFolder;
        nsn = nsn + 1;
        load([outFolder 'weight_corr_eachseq.mat'],'wc1','wc2');
        wc1 = wc1';
        wc2 = wc2';
        %col1 = raw, col2 = exfgcell, col3 = exnfgcell
        wc1_all = cell2mat(wc1);
        wc2_all = cell2mat(wc2);
        WC1_all = [WC1_all;wc1_all];
        WC2_all = [WC2_all;wc2_all];
        wc1_mean = mean(wc1_all);
        wc2_mean = mean(wc2_all);
        
        WC1(nsn,:) = wc1_mean;
        WC2(nsn,:) = wc2_mean;
    end
end


figure(1)
color = {'flat','#F24444','#F2CA50'};
for x = 1:3
[bar,er] = barwitherror(x,WC1(:,x));
bar.FaceColor = color{x};
hold on
end
hold off
xticks([1,2,3])
xticklabels({'raw','FG-cell','NFG-cell'})
ylabel(['Weight correlation'])
set(gca,'FontSize',14)
ylim([0,0.2])

% saveas(gcf,'bar_eachseqWC.png')
% saveas(gcf,'bar_eachseqWC','epsc')

figure(2)
color = {'flat','#F24444','#F2CA50'};
for x = 1:3
[bar,er] = barwitherror(x,WC1_all(:,x));
bar.FaceColor = color{x};
hold on
end
hold off
xticks([1,2,3])
xticklabels({'raw','FG-cell','NFG-cell'})
ylabel(['Weight correlation'])
set(gca,'FontSize',14)
ylim([0,0.25])

saveas(gcf,'bar_eachseqWC_session.png')
saveas(gcf,'bar_eachseqWC_session','epsc')




