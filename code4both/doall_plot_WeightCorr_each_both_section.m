%%
clear
close all
directories_allData_v0_up
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


for L = 2 % 1:2
    nsn = 0;
    wc1_all = cell(3,5);wc2_all = cell(3,5);%3个部分5圈
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
        ff1 = figure('name',path_ns);
        wc1 = cell(3,5);wc2 = cell(3,5);% 3个部分×5lap
        for S = 1:3
            WC1 = nan(5,3);WC2 = nan(5,3);
            
            WC1_m = [];WC2_m = [];WC1_sem = [];WC2_sem = [];
            for    nl = 1:5
                file_input1 = [inFolder,'data_theta_seq_info',subfix,'_lap',num2str(nl),ncases,case3,midmod,'_v5-both-',num2str(S), '.mat'];
                if ~exist(file_input1,'file')
                    WC1_m(nl,1:3) = nan;WC1_sem(nl,1:3) = nan;
                    WC2_m(nl,1:3) = nan;WC2_sem(nl,1:3) = nan;
                    continue
                end
                Allcell = load(file_input1,'TSS','TSS_mean1','TSS_mean2');
                file_input2 = [inFolder,'data_theta_seq_info',subfix2,'_lap',num2str(nl),ncases,case3,midmod,lockat{L},'_v5-both-',num2str(S), '.mat'];
                file_input3 = [inFolder,'data_theta_seq_info',subfix4,'_lap',num2str(nl),ncases,case3,midmod,lockat{L},'_v5-both-',num2str(S), '.mat'];
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
                
                % check TSS size
                thetaSeq1 = Allcell.TSS; S1 = size(thetaSeq1);
                thetaSeq2 = Exfgcell.TSS; S2 = size(thetaSeq2);
                thetaSeq3 = Exnonfgcell.TSS; S3 = size(thetaSeq3);
                if  S1(1) ~= S2(1) || S1(1) ~= S3(1)
                    disp('something wrong')
                    continue
                end
                for tSi = 1:S1(1)
                    wc1{S,nl}(tSi,1) = nanweightcorr(thetaSeq1{tSi,1},TPsweep,SPsweep);
                    wc2{S,nl}(tSi,1) = nanweightcorr(thetaSeq1{tSi,2},TPsweep,SPsweep);
                    wc1{S,nl}(tSi,2) = nanweightcorr(thetaSeq2{tSi,1},TPsweep,SPsweep);
                    wc2{S,nl}(tSi,2) = nanweightcorr(thetaSeq2{tSi,2},TPsweep,SPsweep);
                    wc1{S,nl}(tSi,3) = nanweightcorr(thetaSeq3{tSi,1},TPsweep,SPsweep);
                    wc2{S,nl}(tSi,3) = nanweightcorr(thetaSeq3{tSi,2},TPsweep,SPsweep);
                end
                % 每个session的WC值
                WC1_m(nl,:) = mean(wc1{S,nl});WC1_sem(nl,:) = std(wc1{S,nl})./sqrt(S1(1));
                WC2_m(nl,:) = mean(wc2{S,nl});WC2_sem(nl,:) = std(wc2{S,nl})./sqrt(S1(1));
                % 所有session的WC值
                WC1_m_all{S}(nl,:,ns) = WC1_m(nl,:);
                WC2_m_all{S}(nl,:,ns) = WC2_m(nl,:);
                % 所有session的所有sequence的WC值
                wc1_all{S,nl} = [wc1_all{S,nl};wc1{S,nl}];
                wc2_all{S,nl} = [wc2_all{S,nl};wc2{S,nl}];
            end
            % 画单个sequence的WC
            figure(ff1);
            set(ff1,'Position',[463 172 1413 728]);
            subplot(2,3,S)
            errorbar(WC1_m,WC1_sem,'LineWidth',2)
            box off
            xlim([0.5 5.5])
            ylim([-0.1 0.3])
            
            ylabel('Weight Correlation')
            title('not control speed')
            set(gca,'FontSize',18)
            subplot(2,3,3+S)
            errorbar(WC2_m,WC2_sem,'LineWidth',2)
            box off
            xlim([0.5 5.5])
            ylim([-0.1 0.3])
            xlabel('Lap Num')
            ylabel('Weight Correlation')
            title('control speed')
            set(gca,'FontSize',18)
        end
%         saveas(ff1,[outFolder 'weight_corr_eachseq_part.png']);
%         save([outFolder 'weight_corr_eachseq_part.mat'],'wc1','wc2');
    end
    ff2 = figure;
    for S = 1:3
        for nl = 1:5
            seq_num = size(wc1_all{nl},1);
            wc1_mean{S}(nl,:) = mean(wc1_all{S,nl});wc1_sem{S}(nl,:) = std(wc1_all{S,nl})./sqrt(seq_num);
            wc2_mean{S}(nl,:) = mean(wc2_all{S,nl});wc2_sem{S}(nl,:) = std(wc2_all{S,nl})./sqrt(seq_num);
        end
        % 画全部单个sequence的wc
        set(ff2,'Position',[463 172 1413 728]);
        subplot(2,3,S)
        errorbar(wc1_mean{S},wc1_sem{S},'LineWidth',2)
        box off
        xlim([0.5 5.5])
        ylim([-0.1 0.3])
        ylabel('Weight Correlation')
        title('not control speed')
        set(gca,'FontSize',18)
        subplot(2,3,S+3)
        errorbar(wc2_mean{S},wc2_sem{S},'LineWidth',2)
        box off
        xlim([0.5 5.5])
        ylim([-0.1 0.3])
        xlabel('Lap Num')
        ylabel('Weight Correlation')
        title('control speed')
        set(gca,'FontSize',18)
        legend({'allcell','exfgcell','exnonfgcell'},'Position',[0.875,0.84,0.12,0.12])
%         saveas(ff2,[outFolder 'weight_corr_eachseq_all_part.png']);
    end
end