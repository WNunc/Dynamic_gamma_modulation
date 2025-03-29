% 直接读取
% theta_sequence_averaging_traverse_both.m
% 的输出，统计每只鼠的最好切分相位


clear
close all
directories_allData_v0

Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
case1 = '-ontrack';
ncells = 3;nspkikes = 5;
case2 = sprintf('_%ucell_%uspk',ncells, nspkikes);
midmod = '_cycmid';
TPsweep = 7;
SPsweep = 10;


for rats = [2,3,10,11,16]%鼠的编号
    isession = find(Ind_Rat == rats);
    goodphase1all = [];goodphase2all = [];goodphase3all = [];
    for ns = isession'
        path_ns = path{ns};
        disp(path_ns)
        cd(path_ns)
        
        trackdata_ns = trackdata{ns};
        traverseIND = 0;
        phasestep = 10;
        ProbDiff_stat = [];WeightCorr_stat = [];FlipCorr_stat = [];
        
        while traverseIND < 36
            cut_phase = traverseIND * phasestep;
            traverseIND = traverseIND+1;
            case3 = num2str(cut_phase);
            TSS_all1 = {};
            TSS_all2 = {};
            
            for D = 1:2 %1 =CW, 2 = CCW
                subfolder1 = Directionfolder{D};
                subfolder2 = Phasecutfolder{D};
                for nl = 1:5
                    fileinput = [path_ns,subfolder2,'phase_', case3, '\data_theta_seq_info',...
                        case1,'_lap',num2str(nl),case2,case3,midmod,'_v5.mat'];
                    if ~exist(fileinput,'file')
                        disp([ subfolder2 ': lap_' num2str(nl) ' has no enough sequence'])
                        continue
                    end
                    if D == 1
                        seqCW = load(fileinput);
                        TSS = seqCW.TSS;
                    else
                        seqCCW = load(fileinput);
                        TSS = seqCCW.TSS;
                    end
                    TSS_all1 = cat(1,TSS_all1,TSS{:,1});
                    TSS_all2 = cat(1,TSS_all2,TSS{:,2});
                end
            end
            TSS_temp1 = []; TSS_temp2 = [];
            
            for i = 1:length(TSS_all1)
                TSS_temp1(:,:,i) =  TSS_all1{i};
                TSS_temp2(:,:,i) =  TSS_all2{i};
                %             [WeigCorr1(i),~] = nanweightcorr(TSS_all1{i},TPsweep,SPsweep);
                %             [WeigCorr2(i),~] = nanweightcorr(TSS_all2{i},TPsweep,SPsweep);
            end
            TSS_mean1 = mean(TSS_temp1,3,'omitnan');
            TSS_mean2 = mean(TSS_temp2,3,'omitnan');
            [ProbDiff,~] = probdiff(TSS_mean1,TPsweep,SPsweep);
            [WeigCorr,~] = nanweightcorr(TSS_mean1,TPsweep,SPsweep);
            [FlipCorr,~] = flipcorr(TSS_mean1,TPsweep,SPsweep);
            ProbDiff_stat(traverseIND,:) = ProbDiff;
            WeightCorr_stat(traverseIND,:) = WeigCorr;
            FlipCorr_stat(traverseIND,:) = FlipCorr;
        end
        [~,goodphase1]= maxk(ProbDiff_stat,10);
        goodphase1 = (goodphase1-1).*phasestep;
        [~,goodphase2]= maxk(WeightCorr_stat,10);
        goodphase2 = (goodphase2-1).*phasestep;
        [~,goodphase3]= maxk(FlipCorr_stat,10);
        goodphase3 = (goodphase3-1).*phasestep;
        
        goodphase1all = [goodphase1all; goodphase1];
        goodphase2all = [goodphase2all; goodphase2];
        goodphase3all = [goodphase3all; goodphase3];
    end
    binsize = 10;
    xph = 0:  binsize  :360; %this is set up for the hist which centers the bars
    figure
    set(gcf,'Position',[587 228 678 366])
    H = histogram([goodphase1all,goodphase2all,goodphase3all],xph);
    xlim([0 360])
    xlabel('phase of theta')
    ylabel('number of max indicator')
    goodphase = find(H.Values == max(H.Values));
    goodphase = round(mean(goodphase));
    goodphase = (goodphase-1).*phasestep;
    title(['Rat ' num2str(rats) ' goodphase: ' num2str(goodphase)])
end