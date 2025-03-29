%%
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
    for ns = [1,2,4:6,8,10:13,15,16]%1:isession %[1;2;4;5;6;8;10;11;16;20;22;25]'% [1:11,16,20:22,25,26,28,29,30]%[1:4,6,7]% [1;2;5;6;8;10;11;13;20;22;23;24;27;28]'%[1:8,10:15,22:31]%[1:9,11,13:16]%
        
        path_ns = path{ns};
        cd(path_ns)
        disp(path_ns)
        goodphase_ns = Seq_cutPhase{ns,2};
        case3 = num2str(goodphase_ns);
        
        inFolder = [resuletFolder path_ns(13:end)];
        outFolder = inFolder;
        %         name = path_ns(14:end);
        %         name(name == '\') = '-';
        %         outFolder = ['I:\sequence\each\' name ];
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
        wc1_m = [];wc2_m = [];wc1_sem = [];wc2_sem = [];
        
        a = [];b = [];
        nsn = nsn +1;
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
            
            % check TSS size
            thetaSeq1 = Allcell.TSS; S1 = size(thetaSeq1);
            thetaSeq2 = Exfgcell.TSS; S2 = size(thetaSeq2);
            thetaSeq3 = Exnonfgcell.TSS; S3 = size(thetaSeq3);
            if  S1(1) ~= S2(1) || S1(1) ~= S3(1)
                disp('something wrong')
                continue
            end
            for tSi = 1:S1(1)
                wc1{nl}(tSi,1) = nanweightcorr(thetaSeq1{tSi,1},TPsweep,SPsweep);
                wc2{nl}(tSi,1) = nanweightcorr(thetaSeq1{tSi,2},TPsweep,SPsweep);
                wc1{nl}(tSi,2) = nanweightcorr(thetaSeq2{tSi,1},TPsweep,SPsweep);
                wc2{nl}(tSi,2) = nanweightcorr(thetaSeq2{tSi,2},TPsweep,SPsweep);
                wc1{nl}(tSi,3) = nanweightcorr(thetaSeq3{tSi,1},TPsweep,SPsweep);
                wc2{nl}(tSi,3) = nanweightcorr(thetaSeq3{tSi,2},TPsweep,SPsweep);
            end
            
            
            a = [a;[wc1{nl}(:,1),ones(S1(1),1)*nl]];
            b = [b;[wc2{nl}(:,1),ones(S1(1),1)*nl]];
            % 每个session中，每圈全部sequence的平均WC值
            wc1_m(nl,:) = mean(wc1{nl});wc1_sem(nl,:) = std(wc1{nl})./sqrt(S1(1));
            wc2_m(nl,:) = mean(wc2{nl});wc2_sem(nl,:) = std(wc2{nl})./sqrt(S1(1));
            % 所有session，每圈全部sequence的平均WC值
            wc1_m_all(nl,:,nsn) = wc1_m(nl,:);
            wc2_m_all(nl,:,nsn) = wc2_m(nl,:);
            % 所有session，所有sequence的WC值
            wc1_all{nl} = [wc1_all{nl};wc1{nl}];
            wc2_all{nl} = [wc2_all{nl};wc2{nl}];
            % 所有session的平均sequence的WC值
            WC1_m_all(nl,:,nsn) = WC1(nl,:);
            WC2_m_all(nl,:,nsn) = WC2(nl,:);
        end
        
        % 计算相关性
        nlap = 1:5;
        A(:,1) = wc1_m(:,1);
        A(:,2) = nlap';
        [R_A,P_A] = corr(A(:,1),A(:,2),'type','Pearson')
        [R_a,P_a] = corr(a(:,1),a(:,2),'type','Pearson')
        B(:,1) = wc2_m(:,1);
        B(:,2) = nlap';
        [R_B,P_B] = corr(B(:,1),B(:,2),'type','Pearson')
        [R_b,P_b] = corr(b(:,1),b(:,2),'type','Pearson')
        
        R(nsn,1) = R_A;R(nsn,2) = R_B;
        r(nsn,1) = R_a;r(nsn,2) = R_b;
        P(nsn,1) = P_A;P(nsn,2) = P_B;
        p(nsn,1) = P_a;p(nsn,2) = P_b;
        
        % 画单个sequence的WC
        ff1 = figure('name',path_ns);
        set(ff1,'Position',[960 55 960 420]);
        subplot(1,2,1)
        errorbar(wc1_m,wc1_sem,'LineWidth',2)
        box off
        xlim([0.5 5.5])
        ylim([-0.1 0.3])
        xlabel('Lap Num')
        ylabel('Weight Correlation')
        set(gca,'FontSize',18)
        title('not control speed',...
            ['R = ' num2str(R_A) ',P =' num2str(P_A),'  r = ' num2str(R_a) ',p =' num2str(P_a)],'FontSize',13)
        
        subplot(1,2,2)
        errorbar(wc2_m,wc2_sem,'LineWidth',2)
        box off
        xlim([0.5 5.5])
        ylim([-0.1 0.3])
        xlabel('Lap Num')
        set(gca,'FontSize',18)
        title('control speed',...
            ['R = ' num2str(R_B) ',P =' num2str(P_B),'  r = ' num2str(R_b) ',p =' num2str(P_b)],'FontSize',13)
        %saveas(ff1,[outFolder 'weight_corr_eachseq.png']);
        %save([outFolder 'weight_corr_eachseq.mat'],...
        %    'wc1','wc2','wc1_m','wc2_m','wc1_sem','wc2_sem','WC1','WC2');
    end
    
    for nl = 1:5
        seq_num = size(wc1_all{nl},1);
        wc1_mean(nl,:) = mean(wc1_all{nl});wc1_sem(nl,:) = std(wc1_all{nl})./sqrt(seq_num);
        wc2_mean(nl,:) = mean(wc2_all{nl});wc2_sem(nl,:) = std(wc2_all{nl})./sqrt(seq_num);
    end
    % 画全部单个sequence的wc
    ff2 = figure('name','groupdata_each_sequence');
    set(ff2,'Position',[960 55 960 420]);
    subplot(1,2,1)
    errorbar(wc1_mean,wc1_sem,'LineWidth',2)
    box off
    xlim([0.5 5.5])
    ylim([0 0.3])
    xlabel('Lap Num')
    ylabel('Weight Correlation')
    title('not control speed')
    set(gca,'FontSize',18)
    subplot(1,2,2)
    errorbar(wc2_mean,wc2_sem,'LineWidth',2)
    box off
    xlim([0.5 5.5])
    ylim([0 0.3])
    xlabel('Lap Num')
    title('control speed')
    set(gca,'FontSize',18)
end

%% 整理成SPSS统计的格式

% % cd('H:\neuralynx\gamma in sequence result\theta sequence indicator\230214')
% % 全部圈
% sessionIND = 1:isession;%[1:11,16,20:22,25,26,28,29,30];
% group = [ones(length(sessionIND),1);ones(length(sessionIND),1)*2];
% GPdataA1_fgcell = squeeze(wc1_m_all(:,2,sessionIND));
% GPdataA1_nfgcell = squeeze(wc1_m_all(:,3,sessionIND));
% GPdataA1 = [group,[GPdataA1_fgcell';GPdataA1_nfgcell']];
% GPdataA2_fgcell = squeeze(wc2_m_all(:,2,sessionIND));
% GPdataA2_nfgcell = squeeze(wc2_m_all(:,3,sessionIND));
% GPdataA2 = [group,[GPdataA2_fgcell';GPdataA2_nfgcell']];
% % 画图
% plotWC(GPdataA1_fgcell,GPdataA1_nfgcell,GPdataA2_fgcell,GPdataA2_nfgcell,'全部session')
% save('WC-allsession.mat','group','GPdataA1_fgcell','GPdataA1_nfgcell','GPdataA2_fgcell','GPdataA2_nfgcell')
%%%%%%%%%%%%%%%%%%%%%%%%%%重要%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 前两圈向上**********************************************
sessionIND = 1:12;%[1,2,4:6,8,10:13,15,16]%[2 4 5 10 11 15 16 18 19 20]%;
%group: 1=exfgcell;2=exnfgcell
group = [ones(length(sessionIND),1);ones(length(sessionIND),1)*2];
GPdataB1_fgcell = squeeze(wc1_m_all(:,2,sessionIND));
GPdataB1_nfgcell = squeeze(wc1_m_all(:,3,sessionIND));
GPdataB1 = [group,[GPdataB1_fgcell';GPdataB1_nfgcell']];
GPdataB2_fgcell = squeeze(wc2_m_all(:,2,sessionIND));
GPdataB2_nfgcell = squeeze(wc2_m_all(:,3,sessionIND));
GPdataB2 = [group,[GPdataB2_fgcell';GPdataB2_nfgcell']];

GPdataB1_fgcell0 = squeeze(WC1_m_all(:,2,sessionIND));
GPdataB1_nfgcell0 = squeeze(WC1_m_all(:,3,sessionIND));
GPdataB2_fgcell0 = squeeze(WC2_m_all(:,2,sessionIND));
GPdataB2_nfgcell0 = squeeze(WC2_m_all(:,3,sessionIND));

figure
plot_sem(GPdataB1_fgcell','Laps','Weight correlation','#F24444');
hold on
plot_sem(GPdataB1_nfgcell','Laps','Weight correlation','#F2CA50')
hold off
box off
xlim([0.5 5.5])
ylim([0 0.3])
set(gca,'FontSize',18)
legend({'exfgcell','exnonfgcell'},'Position',[0.85,0.75,0.13,0.20],'Box','off')

% 画图
plotWC(GPdataB1_fgcell,GPdataB1_nfgcell,GPdataB2_fgcell,GPdataB2_nfgcell,'前两圈up')
save('WC-f2lapup.mat','group',...
    'GPdataB1_fgcell','GPdataB1_nfgcell','GPdataB2_fgcell','GPdataB2_nfgcell',...
    'GPdataB1_fgcell0','GPdataB1_nfgcell0','GPdataB2_fgcell0','GPdataB2_nfgcell0')

%%%%%%%%%%%%%%%%%%%%%%%%%重要%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 全部圈向上
% sessionIND = [1,2,4,10,14];
% group = [ones(length(sessionIND),1);ones(length(sessionIND),1)*2];
% GPdataC1_fgcell = squeeze(wc1_m_all(:,2,sessionIND));
% GPdataC1_nfgcell = squeeze(wc1_m_all(:,3,sessionIND));
% GPdataC1 = [group,[GPdataC1_fgcell';GPdataC1_nfgcell']];
% GPdataC2_fgcell = squeeze(wc2_m_all(:,2,sessionIND));
% GPdataC2_nfgcell = squeeze(wc2_m_all(:,3,sessionIND));
% GPdataC2 = [group,[GPdataC2_fgcell';GPdataC2_nfgcell']];
% % 画图
% plotWC(GPdataC1_fgcell,GPdataC1_nfgcell,GPdataC2_fgcell,GPdataC2_nfgcell,'所有圈up')
% save('WC-alllapup.mat','group','GPdataC1_fgcell','GPdataC1_nfgcell','GPdataC2_fgcell','GPdataC2_nfgcell')
% %% 非参数检验
% lapnum = [1,2,3,4,5];
% 
% [r,c] = size(GPdataA1_fgcell);
% LAP = repmat(lapnum,1,2*c);
% TYPE = [ones(1,r*c),ones(1,r*c)*2];
% SRH_GPA1 = [reshape(GPdataA1_fgcell,1,r*c),reshape(GPdataA1_nfgcell,1,r*c)];
% disp('全部圈-not control speed')
% DATA_A1 = [SRH_GPA1',TYPE',LAP'];
% SRH_test(DATA_A1,'TPYE','LAP')
% SRH_GPA2 = [reshape(GPdataA2_fgcell,1,r*c),reshape(GPdataA2_nfgcell,1,r*c)];
% disp('全部圈-control speed')
% DATA_A2 = [SRH_GPA2',TYPE',LAP'];
% SRH_test(DATA_A2,'TPYE','LAP')
% 
% [r,c] = size(GPdataB1_fgcell);
% LAP = repmat(lapnum,1,2*c);
% TYPE = [ones(1,r*c),ones(1,r*c)*2];
% SRH_GPB1 = [reshape(GPdataB1_fgcell,1,r*c),reshape(GPdataB1_nfgcell,1,r*c)];
% disp('前2圈-not control speed')
% DATA_B1 = [SRH_GPB1',TYPE',LAP'];
% SRH_test(DATA_B1,'TPYE','LAP')
% SRH_GPB2 = [reshape(GPdataB2_fgcell,1,r*c),reshape(GPdataB2_nfgcell,1,r*c)];
% disp('前2圈-control speed')
% DATA_B2 = [SRH_GPB2',TYPE',LAP'];
% SRH_test(DATA_B2,'TPYE','LAP')
% 
% [r,c] = size(GPdataC1_fgcell);
% LAP = repmat(lapnum,1,2*c);
% TYPE = [ones(1,r*c),ones(1,r*c)*2];
% SRH_GPC1 = [reshape(GPdataC1_fgcell,1,r*c),reshape(GPdataC1_nfgcell,1,r*c)];
% disp('全up圈-not control speed')
% DATA_C1 = [SRH_GPC1',TYPE',LAP'];
% SRH_test(DATA_C1,'TPYE','LAP')
% SRH_GPC2 = [reshape(GPdataC2_fgcell,1,r*c),reshape(GPdataC2_nfgcell,1,r*c)];
% disp('全up圈-control speed')
% DATA_C2 = [SRH_GPC2',TYPE',LAP'];
% SRH_test(DATA_C2,'TPYE','LAP')



%%
function plotWC(GPdata1_fgcell,GPdata1_nfgcell,GPdata2_fgcell,GPdata2_nfgcell,Name)
[r,c] = size(GPdata1_fgcell);
mean1(:,1) = mean(GPdata1_fgcell,2);sem1(:,1) = std(GPdata1_fgcell,[],2)./sqrt(c);
mean1(:,2) = mean(GPdata1_nfgcell,2);sem1(:,2) = std(GPdata1_nfgcell,[],2)./sqrt(c);
mean2(:,1) = mean(GPdata2_fgcell,2);sem2(:,1)= std(GPdata2_fgcell,[],2)./sqrt(c);
mean2(:,2) = mean(GPdata2_nfgcell,2);sem2(:,2) = std(GPdata2_nfgcell,[],2)./sqrt(c);
ff = figure('name',Name);%'groupdata_each_sequence'
set(ff,'Position',[960 55 960 420]);
subplot(1,2,1)
h1 = errorbar(mean1,sem1,'LineWidth',2);
h1(1).Color = '#F24444';
h1(2).Color = '#F2CA50';
box off
xlim([0.5 5.5])
ylim([0 0.3])
xlabel('Lap Num')
ylabel('Weight Correlation')
title('not control speed')
set(gca,'FontSize',18)
subplot(1,2,2)
h2 = errorbar(mean2,sem2,'LineWidth',2);
h2(1).Color = '#F24444';
h2(2).Color = '#F2CA50';
box off
xlim([0.5 5.5])
ylim([0 0.3])
xlabel('Lap Num')
title('control speed')
set(gca,'FontSize',18)
legend({'exfgcell','exnonfgcell'},'Position',[0.85,0.75,0.13,0.20])
end






