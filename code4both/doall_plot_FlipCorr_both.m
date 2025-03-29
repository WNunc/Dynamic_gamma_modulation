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
lockat = {'firstlap','alllap'};


for L = 1:2
    nsn = 0;
    for ns = 1:isession
        
        path_ns = path{ns};
        cd(path_ns)
        disp(path_ns)
        goodphase_ns = Seq_cutPhase{ns,2};
        case3 = num2str(goodphase_ns);
        
        inFolder = [resuletFolder path_ns(13:end)];
        outFolder = inFolder;
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
        FC1 = [];FC2 = [];
        sd1 = [];sd2 = [];
        sem1 = [];sem2 = [];
        Q = {};
        FC1 = nan(5,3);FC2 = nan(5,3);
        for    nl = 1:5
            file_input1 = [inFolder,'data_theta_seq_info',subfix,'_lap',num2str(nl),ncases,case3,midmod,'_v5-both.mat'];
            Allcell = load(file_input1,'TSS_mean1','TSS_mean2');
            file_input2 = [inFolder,'data_theta_seq_info',subfix2,'_lap',num2str(nl),ncases,case3,midmod,lockat{L},'_v5-both.mat'];
            file_input3 = [inFolder,'data_theta_seq_info',subfix4,'_lap',num2str(nl),ncases,case3,midmod,lockat{L},'_v5-both.mat'];
            if ~exist(file_input2,'file') || ~exist(file_input3,'file')
                break
            end
            Exfgcell = load(file_input2, 'TSS_mean1','TSS_mean2');
            Exnonfgcell = load(file_input3, 'TSS_mean1','TSS_mean2');
            fprintf(1,'input file:\n%s\n%s\n%s\n',file_input1,file_input2,file_input3)
            
            ThetaSeqm1 = Allcell.TSS_mean1([sp1,sp2],[tp1,tp2]);
            ThetaSeqm2 = Allcell.TSS_mean2([sp1,sp2],[tp1,tp2]);
            FC1(nl,1) = flipcorr(ThetaSeqm1,TPsweep,SPsweep);
            FC2(nl,1) = flipcorr(ThetaSeqm2,TPsweep,SPsweep);
            
            ThetaSeqm1 = Exfgcell.TSS_mean1([sp1,sp2],[tp1,tp2]);
            ThetaSeqm2 = Exfgcell.TSS_mean2([sp1,sp2],[tp1,tp2]);
            FC1(nl,2) = flipcorr(ThetaSeqm1,TPsweep,SPsweep);
            FC2(nl,2) = flipcorr(ThetaSeqm2,TPsweep,SPsweep);
            
            ThetaSeqm1 = Exnonfgcell.TSS_mean1([sp1,sp2],[tp1,tp2]);
            ThetaSeqm2 = Exnonfgcell.TSS_mean2([sp1,sp2],[tp1,tp2]);
            FC1(nl,3) = flipcorr(ThetaSeqm1,TPsweep,SPsweep);
            FC2(nl,3) = flipcorr(ThetaSeqm2,TPsweep,SPsweep);
        end
        nsn = nsn+1;
        FC1_all(:,:,nsn) = FC1;
        FC2_all(:,:,nsn) = FC2;
        
        nlap = 1:5;
        ff1 = figure('name',path_ns);
        set(ff1,'Position',[960 55 960 420]);
        subplot(1,2,1)
        plot(FC1,'LineWidth',2)
        box off
        xlim([0.5 5.5])
        ylim([-0.3 0.4])
        xlabel('Lap Num')
        ylabel('Flip Correlation')
        title('not control speed')
        set(gca,'FontSize',18)
        subplot(1,2,2)
        plot(FC2,'LineWidth',2)
        box off
        xlim([0.5 5.5])
        ylim([-0.3 0.4])
        xlabel('Lap Num')
        title('control speed')
        set(gca,'FontSize',18)
        legend({'allcell','exfgcell','exnonfgcell'},'Position',[0.85,0.75,0.13,0.20])
        file_output1 =['Flip Correlation' ncases Qsize case3 midmod lockat{L} '.png'];
        saveas(ff1,[outFolder file_output1])
        save([outFolder 'Flip Correlation' ncases Qsize case3 midmod lockat{L} '.mat'],'FC1','FC2')
    end
    FC1_allmean = mean(FC1_all,3,'omitnan');
    FC2_allmean = mean(FC2_all,3,'omitnan');
    FC1_allsem = std(FC1_all,0,3,'omitnan')./sqrt(size(FC1_all,3));
    FC2_allsem = std(FC2_all,0,3,'omitnan')./sqrt(size(FC2_all,3));
    
    ff2 = figure;
    set(ff2,'Position',[960 55 960 420]);
    subplot(1,2,1)
    errorbar(FC2_allmean,FC2_allsem,'LineWidth',2)
    box off
    xlim([0.5 5.5])
    ylim([-0.25 0.25])
    xlabel('Lap Num')
    ylabel('Flip Correlation')
    title('not control speed')
    set(gca,'FontSize',18)
    subplot(1,2,2)
    errorbar(FC2_allmean,FC2_allsem,'LineWidth',2)
    box off
    xlim([0.5 5.5])
    ylim([-0.25 0.25])
    xlabel('Lap Num')
    title('control speed')
    set(gca,'FontSize',18)
    legend({'allcell','exfgcell','exnonfgcell'},'Position',[0.85,0.75,0.13,0.20])
    saveas(ff2,['H:\neuralynx\gamma in sequence result\Flip Correlation_' lockat{L} '_both.png'])
    save(['H:\neuralynx\gamma in sequence result\Flip Correlation_' lockat{L} '_both.mat'],...
        'FC1_all','FC2_all')
end
