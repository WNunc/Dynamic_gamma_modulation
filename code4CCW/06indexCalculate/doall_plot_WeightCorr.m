5%%
clear
close all
directories_allData_v2
subfolder = 'Tseq-CCW\';
subfix = '-ontrack'; % v5
subfix1 = '-ontrack_exsg'; %exclude sgamma
subfix2 = '-ontrack_exfg'; %exclude fgamma
subfix3 = '-ontrack_dssg'; %down sample sgamma
subfix4 = '-ontrack_dsfg'; %down sample fgamma
ncells = 3;nspkikes = 5;
% ncases = sprintf('_%ubar_%ucell_%uspk',nbars, ncells, nspkikes);
ncases = sprintf('_%ucell_%uspk', ncells, nspkikes);
midmod = '_cycmid';
lockat = {'firstlap','alllap'};


for L = 1:2
    nsn = 0;
    for ns = 1:isession
        
        path_ns = path{ns};
        cd(path_ns)
        disp(path_ns)
        goodphase_ns = Seq_cutPhase{ns,2};%CCW是第二列
        case3 = num2str(goodphase_ns);
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
        for    nl = 1:5
            file_input1 = [path_ns,subfolder,'data_theta_seq_info',subfix,'_lap',num2str(nl),ncases,case3,midmod,'_v5.mat'];
            Allcell = load(file_input1,'TSS_mean1','TSS_mean2');
            file_input2 = [path_ns,subfolder,'data_theta_seq_info',subfix2,'_lap',num2str(nl),ncases,case3,midmod,lockat{L},'.mat'];
            file_input3 = [path_ns,subfolder,'data_theta_seq_info',subfix4,'_lap',num2str(nl),ncases,case3,midmod,lockat{L},'.mat'];
            if ~exist(file_input2,'file') || ~exist(file_input3,'file')
                break
            end
            Exfgcell = load(file_input2, 'TSS_mean1','TSS_mean2');
            Exnonfgcell = load(file_input3, 'TSS_mean1','TSS_mean2');
            fprintf(1,'input file:\n%s\n%s\n%s\n',file_input1,file_input2,file_input3)
            
            ThetaSeqm1 = Allcell.TSS_mean1([sp1,sp2],[tp1,tp2]);
            ThetaSeqm2 = Allcell.TSS_mean2([sp1,sp2],[tp1,tp2]);
            WC1(nl,1) = weightcorr(ThetaSeqm1,TPsweep,SPsweep);
            WC2(nl,1) = weightcorr(ThetaSeqm2,TPsweep,SPsweep);
            
            ThetaSeqm1 = Exfgcell.TSS_mean1([sp1,sp2],[tp1,tp2]);
            ThetaSeqm2 = Exfgcell.TSS_mean2([sp1,sp2],[tp1,tp2]);
            WC1(nl,2) = weightcorr(ThetaSeqm1,TPsweep,SPsweep);
            WC2(nl,2) = weightcorr(ThetaSeqm2,TPsweep,SPsweep);
            
            ThetaSeqm1 = Exnonfgcell.TSS_mean1([sp1,sp2],[tp1,tp2]);
            ThetaSeqm2 = Exnonfgcell.TSS_mean2([sp1,sp2],[tp1,tp2]);
            WC1(nl,3) = weightcorr(ThetaSeqm1,TPsweep,SPsweep);
            WC2(nl,3) = weightcorr(ThetaSeqm2,TPsweep,SPsweep);
        end
        nsn = nsn+1;
        WC1_all(:,:,nsn) = WC1;
        WC2_all(:,:,nsn) = WC2;
        
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
        title('not control speed')
        set(gca,'FontSize',18)
        subplot(1,2,2)
        plot(WC2,'LineWidth',2)
        box off
        xlim([0.5 5.5])
        ylim([-0.1 0.3])
        xlabel('Lap Num')
        title('control speed')
        set(gca,'FontSize',18)
        legend({'allcell','exfgcell','exnonfgcell'},'Position',[0.85,0.75,0.13,0.20])
        file_output1 =['Weight_Correlation' ncases Qsize case3 midmod lockat{L} '.png'];
%         saveas(ff1,[subfolder file_output1])
    end
    WC1_allmean = mean(WC1_all,3,'omitnan');
    WC2_allmean = mean(WC2_all,3,'omitnan');
    WC1_allsem = std(WC1_all,0,3,'omitnan')./sqrt(size(WC1_all,3));
    WC2_allsem = std(WC2_all,0,3,'omitnan')./sqrt(size(WC2_all,3));
    
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
    save(['H:\neuralynx\gamma in sequence result\Weight Correlation_CCW_' lockat{L} '_0627'],...
        'WC1_all','WC2_all')
end
