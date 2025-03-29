clear
close all

result_folder = 'H:\neuralynx\gamma in sequence result\';
lockat = {'firstlap','alllap'};
indexname = {'Probability Difference',...
    'Weight Correlation',...
    'Flip Correlation'};
direction = {'_CW_','_CCW_'};
for L = 1:2
    for I = 1:3
        result_CW = [result_folder indexname{I} '_CW_' lockat{L}];
        result_CCW = [result_folder indexname{I} '_CCW_' lockat{L}];
        switch I
            case 1
                load(result_CW);
                all1_CW = PD1_all;
                all2_CW = PD2_all;
                clear PD1_all PD2_all
                load(result_CCW);
                all1_CCW = PD1_all;
                all2_CCW = PD2_all;
                clear PD1_all PD2_all
            case 2
                load(result_CW);
                all1_CW = WC1_all;
                all2_CW = WC2_all;
                Cclear WC1_all WC2_all
                load(result_CCW);
                all1_CCW = WC1_all;
                all2_CCW = WC2_all;
                clear WC1_all WC2_all
            case 3
                load(result_CW);
                all1_CW = FC1_all;
                all2_CW = FC2_all;
                clear FC1_all FC2_all
                load(result_CCW);
                all1_CCW = FC1_all;
                all2_CCW = FC2_all;
                clear FC1_all FC2_all
        end
        
        
        all1 = cat(3,all1_CW,all1_CCW);
        all2 = cat(3,all2_CW,all2_CCW);
        
        SPSS_all1 = [squeeze(all1(:,1,:)) squeeze(all1(:,2,:)) squeeze(all1(:,3,:))]';
        SPSS_all2 = [squeeze(all2(:,1,:)) squeeze(all2(:,2,:)) squeeze(all2(:,3,:))]';
        
        all1mean = mean(all1,3,'omitnan');
        all2mean = mean(all2,3,'omitnan');
        all1sem = std(all1,0,3,'omitnan')./sqrt(size(all1,3));
        all2sem = std(all2,0,3,'omitnan')./sqrt(size(all2,3));
        
        ff2 = figure('name',lockat{L});
        set(ff2,'Position',[960 55 960 420]);
        subplot(1,2,1)
        errorbar(all1mean,all1sem,'LineWidth',2)
        box off
        xlim([0.5 5.5])
        ylim([0 0.3])
        if I==3
            ylim([-0.25 0.25])
        end
        xlabel('Lap Num')
        ylabel(indexname{I})
        title('not control speed')
        set(gca,'FontSize',18)
        subplot(1,2,2)
        errorbar(all2mean,all2sem,'LineWidth',2)
        box off
        xlim([0.5 5.5])
        ylim([0 0.3])
        if I==3
            ylim([-0.25 0.25])
        end
        xlabel('Lap Num')
        title('control speed')
        set(gca,'FontSize',18)
        legend({'allcell','exfgcell','exnonfgcell'},'Position',[0.85,0.75,0.13,0.20])
        
        AA = all1(:,2,:)-all1(:,1,:);
        AA = squeeze(AA);
        AAm = mean(AA,2,'omitnan');
        
        BB = all1(:,3,:)-all1(:,1,:);
        BB = squeeze(BB);
        BBm = mean(BB,2,'omitnan');
        figure
        subplot(1,2,1)
        plot([zeros(5,1),AAm,BBm],'LineWidth',2);
        
        legend({'allcell','exfgcell','exnonfgcell'},'Position',[0.85,0.75,0.13,0.20])
        ylim([-0.1 0.05])
        
        AA = all2(:,2,:)-all2(:,1,:);
        AA = squeeze(AA);
        AAm = mean(AA,2,'omitnan');
        
        BB = all2(:,3,:)-all2(:,1,:);
        BB = squeeze(BB);
        BBm = mean(BB,2,'omitnan');
        subplot(1,2,2)
        plot([zeros(5,1),AAm,BBm],'LineWidth',2);
        
        legend({'allcell','exfgcell','exnonfgcell'},'Position',[0.85,0.75,0.13,0.20])
        ylim([-0.1 0.05])
        
    end
end




 
