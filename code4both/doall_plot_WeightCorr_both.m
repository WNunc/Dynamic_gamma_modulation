%%
clear
close all
directories_allData_v0
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
% 第一圈相锁，全部圈相锁，前两圈相锁
% 最后一版结果用的是全部圈相锁


for L = 2 %1:2
    nsn = 0;
    for ns = [1;2;4;5;6;8;10;11;13;14;15;16;20;22;24;25;27]'%1:isession%[1;2;4;5;6;8;10;11;13;14;15;16;20;22;23;24;26;27;28]'%
        
        path_ns = path{ns};
        cd(path_ns)
        disp(path_ns)
        goodphase_ns = Seq_cutPhase{ns,2};
        case3 = num2str(goodphase_ns);
        
        inFolder = [resuletFolder path_ns(13:end)];
%         outFolder = inFolder;
        name = path_ns(14:end);
        name(name == '\') = '-';
        outFolder = ['I:\sequence\avg\' name ];
        
        TP = 68;TPmid = 0.5*TP;
        SP = 36;SPmid = 0.5*SP;
        TPsweep = 7;
        SPsweep = 10;
        Qsize = sprintf('_%utbins_%uxbins',TPsweep,SPsweep);
        tp1 = [TPmid-TPsweep:TPmid];
        tp2 = [TPmid+1:TPmid+1+TPsweep];
        sp1 = [SPmid-SPsweep:SPmid];
        sp2 = [SPmid+1:SPmid+1+SPsweep];
        %%
        PD1 = [];PD2 = [];
        sd1 = [];sd2 = [];
        sem1 = [];sem2 = [];
        Q = {};
        WC1 = nan(5,3);WC2 = nan(5,3);
        wc1 = cell(1,5);wc2 = cell(1,5);
        WC1_m = [];WC2_m = [];WC1_sem = [];WC2_sem = [];
        for    nl = 1:5
            file_input1 = [inFolder,'data_theta_seq_info',subfix,'_lap',num2str(nl),ncases,case3,midmod,'_v5-both.mat'];
            Allcell = load(file_input1,'TSS','TSS_mean1','TSS_mean2');
            file_input2 = [inFolder,'data_theta_seq_info',subfix2,'_lap',num2str(nl),ncases,case3,midmod,lockat{L},'_v5-bothv2.mat'];
            file_input3 = [inFolder,'data_theta_seq_info',subfix4,'_lap',num2str(nl),ncases,case3,midmod,lockat{L},'_v5-bothv3.mat'];
            if ~exist(file_input2,'file') || ~exist(file_input3,'file')
                break
            end
            Exfgcell = load(file_input2, 'TSS','TSS_mean1','TSS_mean2');
            Exnonfgcell = load(file_input3, 'TSS','TSS_mean1','TSS_mean2');
            fprintf(1,'input file:\n%s\n%s\n%s\n',file_input1,file_input2,file_input3)
            
            ThetaSeqm1 = Allcell.TSS_mean1([sp1,sp2],[tp1,tp2]);
            ThetaSeqm2 = Allcell.TSS_mean2([sp1,sp2],[tp1,tp2]);
            WC1(nl,1) = nanweightcorr(ThetaSeqm1,TPsweep,SPsweep);
            WC2(nl,1) = nanweightcorr(ThetaSeqm2,TPsweep,SPsweep);
            
            ThetaSeqm1 = Exfgcell.TSS_mean1([sp1,sp2],[tp1,tp2]);
            ThetaSeqm2 = Exfgcell.TSS_mean2([sp1,sp2],[tp1,tp2]);
            WC1(nl,2) = nanweightcorr(ThetaSeqm1,TPsweep,SPsweep);
            WC2(nl,2) = nanweightcorr(ThetaSeqm2,TPsweep,SPsweep);
            
            ThetaSeqm1 = Exnonfgcell.TSS_mean1([sp1,sp2],[tp1,tp2]);
            ThetaSeqm2 = Exnonfgcell.TSS_mean2([sp1,sp2],[tp1,tp2]);
            WC1(nl,3) = nanweightcorr(ThetaSeqm1,TPsweep,SPsweep);
            WC2(nl,3) = nanweightcorr(ThetaSeqm2,TPsweep,SPsweep);
%             
%             % check TSS size
%             thetaSeq1 = Allcell.TSS; S1 = size(thetaSeq1);
%             thetaSeq2 = Exfgcell.TSS; S2 = size(thetaSeq2);
%             thetaSeq3 = Exnonfgcell.TSS; S3 = size(thetaSeq3);
%             if  S1(1) ~= S2(1) || S1(1) ~= S3(1)
%                 disp('something wrong')
%                 continue
%             end
%             
%             
%             for tSi = 1:S1(1)
%                 wc1{nl}(tSi,1) = nanweightcorr(thetaSeq1{tSi,1},TPsweep,SPsweep);
%                 wc2{nl}(tSi,1) = nanweightcorr(thetaSeq1{tSi,2},TPsweep,SPsweep);
%                 wc1{nl}(tSi,2) = nanweightcorr(thetaSeq2{tSi,1},TPsweep,SPsweep);
%                 wc2{nl}(tSi,2) = nanweightcorr(thetaSeq2{tSi,2},TPsweep,SPsweep);
%                 wc1{nl}(tSi,3) = nanweightcorr(thetaSeq3{tSi,1},TPsweep,SPsweep);
%                 wc2{nl}(tSi,3) = nanweightcorr(thetaSeq3{tSi,2},TPsweep,SPsweep);
%             end
%             WC1_m(nl,:) = mean(wc1{nl});WC1_sem(nl,:) = std(wc1{nl})./sqrt(S1(1));
%             WC2_m(nl,:) = mean(wc2{nl});WC2_sem(nl,:) = std(wc2{nl})./sqrt(S1(1));
        end
        nlap = 1:5;
        type = 'Pearson';%'Spearman';'Pearson'
        [R1,P1] = corr(WC1(nlap,1),nlap','type',type)
        [R2,P2] = corr(WC2(nlap,1),nlap','type',type)
%         [R1,P1] = corrcoef(WC1(nlap,1),nlap')
%         [R2,P2] = corrcoef(WC2(nlap,1),nlap')
%         % 画单个sequence的WC
%         ff3 = figure('name',path_ns);
%         set(ff3,'Position',[960 55 960 420]);
%         subplot(1,2,1)
%         errorbar(WC1_m,WC1_sem,'LineWidth',2)
%         box off
%         xlim([0.5 5.5])
%         ylim([-0.1 0.3])
%         xlabel('Lap Num')
%         ylabel('Weight Correlation')
%         title('not control speed')
%         set(gca,'FontSize',18)
%         subplot(1,2,2)
%         errorbar(WC2_m,WC2_sem,'LineWidth',2)
%         box off
%         xlim([0.5 5.5])
%         ylim([-0.1 0.3])
%         xlabel('Lap Num')
%         title('control speed')
%         set(gca,'FontSize',18)
        
        nsn = nsn+1;
        WC1_all(:,:,nsn) = WC1;
        WC2_all(:,:,nsn) = WC2;
        %画平均之后的WC
        nlap = 1:5;
        ff1 = figure('name',path_ns);
        set(ff1,'Position',[960 55 960 420]);
        subplot(1,2,1)
        plot(WC1,'LineWidth',2)
        box off
        xlim([0.5 5.5])
        ylim([-0.1 0.3])
        xlabel('Lap Num')
        ylabel('Weight Correlation')
        set(gca,'FontSize',18)
        title('not control speed',['R = ' num2str(R1) ',P =' num2str(P1)])

        subplot(1,2,2)
        plot(WC2,'LineWidth',2)
        box off
        xlim([0.5 5.5])
        ylim([-0.1 0.3])
        xlabel('Lap Num')
        set(gca,'FontSize',18)
        title('not control speed',['R = ' num2str(R1) ',P =' num2str(P1)])
        legend({'allcell','exfgcell','exnonfgcell'},'Position',[0.85,0.75,0.13,0.20])
        file_output1 =[type '-Weight_Correlation' ncases Qsize case3 midmod lockat{L} '.png'];
%         saveas(ff1,[outFolder file_output1])
%         save([outFolder,'Weight_Correlation' ncases Qsize case3 midmod lockat{L} '.mat'],'WC1','WC2')

    end
    WC1_allmean = mean(WC1_all,3,'omitnan');
    WC2_allmean = mean(WC2_all,3,'omitnan');
    WC1_allsem = std(WC1_all,0,3,'omitnan')./sqrt(size(WC1_all,3));
    WC2_allsem = std(WC2_all,0,3,'omitnan')./sqrt(size(WC2_all,3));
    
    % 画总的平均的WC
    ff2 = figure;
    set(ff2,'Position',[960 55 960 420]);
    subplot(1,2,1)
    errorbar(WC1_allmean,WC1_allsem,'LineWidth',2)
    box off
    xlim([0.5 5.5])
    ylim([0 0.3])
    xlabel('Lap Num')
    ylabel('Weight Correlation')
    title('not control speed')
    set(gca,'FontSize',18)
    subplot(1,2,2)
    errorbar(WC2_allmean,WC2_allsem,'LineWidth',2)
    box off
    xlim([0.5 5.5])
    ylim([0 0.3])
    xlabel('Lap Num')
    title('control speed')
    set(gca,'FontSize',18)
    legend({'allcell','exfgcell','exnonfgcell'},'Position',[0.85,0.75,0.13,0.20])
%     saveas(ff2,['H:\neuralynx\gamma in sequence result\Weight Correlation_' lockat{L} '_both.png'])
%     save(['H:\neuralynx\gamma in sequence result\Weight Correlation_' lockat{L} '_both.mat'],...
%         'WC1_all','WC2_all')
end
