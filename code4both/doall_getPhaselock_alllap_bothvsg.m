% 观察sg相锁，修改自
% doall_getPhaselock_alllap_bothv2.m



clear
close all
% directories_allData_v2
directories_allData_v0_allgood
peak0 = 1; % firing rate threshold to remove place cells
nsn = 0;
TTList0 = 'TTList_dCA1_pyr.txt';
resuletFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};

ffs = figure;
for ns = [1,2,4:6,8,10:13,15,16]%[1:4,6,7]%[1;2;4;5;6;8;10;11;16;20;22;25]'%[1:11,16,20:22,26,27] %1:isession%[1;2;5;6;8;10;11;13;20;22;23;24;27;28]'%
    nsn = nsn+1;
    path_ns = path{ns};
    outFolder = [resuletFolder path_ns(13:end)];
    cd(path_ns);
    Ncell = getnumberofcells_cz_v1(TTList0);
    trackdata_ns = trackdata{ns};load(trackdata_ns,'Ang_RewardLoc_ontrack')
    sgcell = 0;sgcellinfo = {};
    for D = 1:2
        
        subfolder1 = Directionfolder{D};
        subfolder2 = Phasecutfolder{D};
        
        file_input1 = [subfolder1 'Cells_allsegment_v1_vel_0.mat'];  % use to get all spikes
        file_input2 = [subfolder1 'Cells_allsegment_v1_vel_5.mat'];  % use to get ratemap
        file_input3 = [subfolder1 'Cells_singleLap_v2_vel_5.mat'];  % use itself as the decoder
        file_input4 = [subfolder1 'data_phaselocking_TSlap_vel0_new.mat'];  % ** phase locking feature trial TS 用vel0 spike 的结果,不带new也可以
        file_input5 = [subfolder1 'data_phaselocking_TSlap_vel0_f2lap.mat'];  % ** phase locking feature trial TS 用vel0 spike 的前两圈相锁的结果
        % file_input5 = [subfolder1 'data_phaselocking_spkmin.mat'];  % phase locking feature running TS
        fprintf(1,'all spike:\t%s\n%s\ndecoder:\t%s\nphase lock:\t%s\n%s\n',file_input1,file_input2,file_input3,file_input4,file_input5)
        load(file_input4)
        load(file_input5)
        if Ncell < 1
            continue
        end
        disp(pwd)
        for nt = 1% this file for pre_running only
            switch nt
                case 1
                    laps=5;
                case 2
                    laps=8;
                case 3
                    laps=8;
                case 4
                    laps=6;
            end
            % Load Spikes data
            S1 = load(file_input1);  % used to get all spikes
            spikes = S1.spikes;
            S2 = load(file_input2);  % used to get ratemap
            Ratemap_seg = S2.Ratemap_seg{nt};
            % Remove the place cells whose peak firing rate<1Hz in all laps
            peak_all = max(Ratemap_seg); % use S2 for detecting
            ind = find(peak_all >= peak0);
            spikes = spikes(ind,:);
            % find cell on track as detecting
            X = cell2struct(S2.fieldProp_seg{nt}(ind),'placefield',1);
            TT0t64 = ReadFileList(S2.TTList0);
            TT0t64 = TT0t64(ind);
            COM = [];
            for ncell  = 1:length(ind)
                COM(ncell) = X(ncell).placefield(1).x_COM;
            end
            
            cind_ot = find(COM>Ang_RewardLoc_ontrack(1) & COM<Ang_RewardLoc_ontrack(18));
            cellnumontrack = length(cind_ot);
            sgamma_rayleighP = [plFeature_alllap{1,2}.rayleighP];% all lap 算出来的相锁
%             sgamma_rayleighP = [plFeature_f2lap{1,2}.rayleighP];% 前两圈 lap 算出来的相锁
            sgamma_rayleighP_ontrack = sgamma_rayleighP(ind);
            sgamma_rayleighP_ontrack = sgamma_rayleighP_ontrack(cind_ot);
            
            sgPhase_singlap_ot = spikePhase_singlap(ind,2,:);
            sgPhase_singlap_ot = sgPhase_singlap_ot(cind_ot,1,:);
            sgPhase_singlap_ot = squeeze(sgPhase_singlap_ot);
            [spkNum_singlap,spkNum_ind] = spknumfinder(sgPhase_singlap_ot,10);
            spkNum_alllap = sum(spkNum_singlap,2);
            cind_ot_sg = find(sgamma_rayleighP_ontrack<0.05 & spkNum_ind'==1)%find(sgamma_rayleighP_ontrack<0.05 & ) % ontrack的ind
            cellnum_sg(nsn,D) = length(cind_ot_sg);
        end
        if cellnum_sg(nsn,D)==0
            disp([path_ns '没有slowgamma相锁的神经元'])
            continue
        end
        
        
        
        
        %% 组合
        icell = cellnum_sg(nsn,D);
        phasealllap = spikePhase_alllap(ind, 2);
        phasealllap = phasealllap(cind_ot,:);
        anglealllap = [plFeature_alllap{2}.meanAngle];
        anglealllap = anglealllap(ind);
        anglealllap = anglealllap(cind_ot);
        vlengthalllap = [plFeature_alllap{2}.vectorLength];
        vlengthalllap = vlengthalllap(ind);
        vlengthalllap = vlengthalllap(cind_ot);
        fig = figure;
        for isgcell = 1:icell
            %             set(ff, 'Position',[0 650 350*icell 350])
            cellid = cind_ot_sg(isgcell);
            ffp = figure('Position',[500 500 300 450]);
            subplot(2,1,1)
            ph = polarhistogram(phasealllap{cellid}*pi/180,(0:15:360)*pi/180,'Normalization','pdf');
            Ua = anglealllap(cellid);Va = vlengthalllap(cellid);
            title(TT0t64{cellid})
            subplot(2,1,2)
            h = histogram([phasealllap{cellid},phasealllap{cellid}+360],0:15:720);
            xlim([0,720])
            xticks([0,180,360,540,720])
            title(['mAngle = ' num2str(Ua) ', vLength = ' num2str(Va)])
            saveas(ffp,[outFolder 'SG ' TT0t64{cellid} '-' num2str(D) '.png']);
            saveas(ffp,[outFolder 'SG ' TT0t64{cellid} '-' num2str(D) '.eps'],'epsc');
            for nl = 1:5
                sgamma_mAngle = [plFeature_singlap{nl,3}.meanAngle];
                sgamma_mAngle(sgamma_mAngle == 0) = nan;
                sgamma_mAngle_ontrack = sgamma_mAngle(ind);
                sgamma_mAngle_ontrack = sgamma_mAngle_ontrack(cind_ot);
                
                sgamma_vLength = [plFeature_singlap{nl,3}.vectorLength];
                sgamma_vLength(sgamma_vLength == 0|sgamma_vLength == 1) = nan;
                sgamma_vLength_ontrack = sgamma_vLength(ind);
                sgamma_vLength_ontrack = sgamma_vLength_ontrack(cind_ot);
                
                sgamma_rP = [plFeature_singlap{nl,3}.rayleighP];
                sgamma_rP(sgamma_rP == 0) = nan;
                sgamma_rP_ontrack = sgamma_rP(ind);
                sgamma_rP_ontrack = sgamma_rP_ontrack(cind_ot);
                lineWid(nl) = sgamma_rP_ontrack(cellid);% P值作为线宽
                
                U(nl) = sgamma_mAngle_ontrack(cellid);
                V(nl) = sgamma_vLength_ontrack(cellid);
                
            end
            sgcell = sgcell+1;
            [x,y] = pol2cart(U/180*pi,V);
            sgcellinfo{sgcell,1} = U;
            sgcellinfo{sgcell,2} = V;
            sgcellinfo{sgcell,3} = lineWid;
            lineWid(lineWid>0.05 | isnan(lineWid)) = 0.5;% P值作为线宽 不显著的用细线
            lineWid(lineWid<0.05) = 2;% P值作为线宽 显著的用粗线
            figure(fig)
            subplot(1,icell,isgcell)
            c = compass(x,y);
            c(1).LineWidth = lineWid(1);c(2).LineWidth = lineWid(2);c(3).LineWidth = lineWid(3);
            c(4).LineWidth = lineWid(4);c(5).LineWidth = lineWid(5);
            c(1).Color =  [1,0,0];c(2).Color =  [1,0,1];c(3).Color =  [0,1,0];
            c(4).Color =  [0,1,1];c(5).Color =  [0,0,1];
            
        end
        suptitle(path_ns)
    
    legend({'lap1','lap2','lap3','lap4','lap5'},'Position',[0.93,0.7,0.05,0.2],'FontSize',12)
    end
    %         nsn = nsn+1;
    if isempty(sgcellinfo)
        continue
    end
    allcellvL_sg = cell2mat(sgcellinfo(:,2));
    vecLength_slap(nsn,:) = mean(allcellvL_sg,'omitnan');% singlelap
    vecLength_flap(nsn,1) = mean(vecLength_slap(nsn,1:2),'omitnan');% first lap to other
%     vecLength_flap(nsn,2) = vecLength_slap(nsn,2);% first lap to other
    vecLength_flap(nsn,2) = mean(vecLength_slap(nsn,3:5),'omitnan');% first lap to other
    
    figure(ffs)
    set(gcf, 'Position',[0 150 800 400])
    subplot(1,2,1)
    plot(vecLength_slap','k','LineWidth',1.5)
    box off
    ylim([0,1])
    xlim([0.5,5.5])
    ylabel('vector length');xlabel('Lap Num')
    set(gca,'FontSize',14)
    subplot(1,2,2)
    plot(vecLength_flap','k','LineWidth',1.5)
    box off
    ylim([0,1])
    xlim([0.5,2.5])
    ylabel('vector length');xlabel('Lap Num')
    xticks([1,2])
    xticklabels({'1','other'})
    set(gca,'FontSize',14)
%     close all
    
end
sum(sum(cellnum_sg))



function [N,L] = spknumfinder(spike, threshold)
% 找出来每一圈数量不超过threshold的spike
% 输出每一圈spike数量和每一圈满足要求的cell的逻辑值
[i,j] = size(spike);
for I = 1:i
    for J = 1:j
        N(I,J) = length(spike{I,J});
        
        if N(I,J)>=threshold
            L(I,J) = 1;
        else
            L(I,J) = 0;
        end
    end
end
L = sum(L,2);
L = L == j;
end