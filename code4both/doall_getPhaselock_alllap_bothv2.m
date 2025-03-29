% 观察挑选pre前5圈计算得到的相锁神经元在每一圈的相锁变化



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
PreferPhase_allcell = [];
PreferPhase_fgcell = [];

ffs = figure;
for ns = [1,2,4:6,8,10:13,15,16]%[1:4,6,7]%[1;2;4;5;6;8;10;11;16;20;22;25]'%[1:11,16,20:22,26,27] %1:isession%[1;2;5;6;8;10;11;13;20;22;23;24;27;28]'%
    nsn = nsn+1;
    path_ns = path{ns};
    outFolder = [resuletFolder path_ns(13:end)];
    cd(path_ns);
    Ncell = getnumberofcells_cz_v1(TTList0);
    trackdata_ns = trackdata{ns};load(trackdata_ns,'Ang_RewardLoc_ontrack')
    fgcell = 0;fgcellinfo = {};
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
            fgamma_rayleighP = [plFeature_alllap{1,3}.rayleighP];% all lap 算出来的相锁
%             fgamma_rayleighP = [plFeature_f2lap{1,3}.rayleighP];% 前两圈 lap 算出来的相锁
            fgamma_rayleighP_ontrack = fgamma_rayleighP(ind);
            fgamma_rayleighP_ontrack = fgamma_rayleighP_ontrack(cind_ot);
            
            fgPhase_singlap_ot = spikePhase_singlap(ind,3,:);
            fgPhase_singlap_ot = fgPhase_singlap_ot(cind_ot,1,:);
            fgPhase_singlap_ot = squeeze(fgPhase_singlap_ot);
            [spkNum_singlap,spkNum_ind] = spknumfinder(fgPhase_singlap_ot,10);
            spkNum_alllap = sum(spkNum_singlap,2);
            cind_ot_fg = find(fgamma_rayleighP_ontrack<0.05 & spkNum_ind'==1)%find(fgamma_rayleighP_ontrack<0.05 & ) % ontrack的ind
            cellnum_fg(nsn,D) = length(cind_ot_fg);
        end
        if cellnum_fg(nsn,D)==0
            disp([path_ns '没有fastgamma相锁的神经元'])
            continue
        end
        
        
        
        
        %% 组合
        icell = cellnum_fg(nsn,D);
        phasealllap = spikePhase_alllap(ind, 3);
        phasealllap = phasealllap(cind_ot,:);
        anglealllap = [plFeature_alllap{3}.meanAngle];
        anglealllap = anglealllap(ind);
        anglealllap = anglealllap(cind_ot);
        PreferPhase_allcell = [PreferPhase_allcell,anglealllap];
        vlengthalllap = [plFeature_alllap{3}.vectorLength];
        vlengthalllap = vlengthalllap(ind);
        vlengthalllap = vlengthalllap(cind_ot);
        TT0t64 = TT0t64(cind_ot);
        fig = figure;
        for ifgcell = 1:icell
            %             set(ff, 'Position',[0 650 350*icell 350])
            cellid = cind_ot_fg(ifgcell);
            ffp = figure('Position',[500 500 300 450]);
            subplot(2,1,1)
            ph = polarhistogram(phasealllap{cellid}*pi/180,(0:15:360)*pi/180,'Normalization','pdf');
            Ua = anglealllap(cellid);Va = vlengthalllap(cellid);
            PreferPhase_fgcell = [PreferPhase_fgcell,Ua];
            title(TT0t64{cellid})
            subplot(2,1,2)
            h = histogram([phasealllap{cellid},phasealllap{cellid}+360],0:15:720);
            xlim([0,720])
            xticks([0,180,360,540,720])
            title(['mAngle = ' num2str(Ua) ', vLength = ' num2str(Va)])
            set(ffp,'Renderer','painters');
%             saveas(ffp,[outFolder 'FG ' TT0t64{cellid} '-' num2str(D) '.png']);
%             saveas(ffp,[outFolder 'FG ' TT0t64{cellid} '-' num2str(D) '.eps'],'epsc');
            for nl = 1:5
                fgamma_mAngle = [plFeature_singlap{nl,3}.meanAngle];
                fgamma_mAngle(fgamma_mAngle == 0) = nan;
                fgamma_mAngle_ontrack = fgamma_mAngle(ind);
                fgamma_mAngle_ontrack = fgamma_mAngle_ontrack(cind_ot);
                
                fgamma_vLength = [plFeature_singlap{nl,3}.vectorLength];
                fgamma_vLength(fgamma_vLength == 0|fgamma_vLength == 1) = nan;
                fgamma_vLength_ontrack = fgamma_vLength(ind);
                fgamma_vLength_ontrack = fgamma_vLength_ontrack(cind_ot);
                
                fgamma_rP = [plFeature_singlap{nl,3}.rayleighP];
                fgamma_rP(fgamma_rP == 0) = nan;
                fgamma_rP_ontrack = fgamma_rP(ind);
                fgamma_rP_ontrack = fgamma_rP_ontrack(cind_ot);
                lineWid(nl) = fgamma_rP_ontrack(cellid);% P值作为线宽
                
                U(nl) = fgamma_mAngle_ontrack(cellid);
                V(nl) = fgamma_vLength_ontrack(cellid);
                
            end
            fgcell = fgcell+1;
            [x,y] = pol2cart(U/180*pi,V);
            fgcellinfo{fgcell,1} = U;
            fgcellinfo{fgcell,2} = V;
            fgcellinfo{fgcell,3} = lineWid;
            lineWid(lineWid>0.05 | isnan(lineWid)) = 0.5;% P值作为线宽 不显著的用细线
            lineWid(lineWid<0.05) = 2;% P值作为线宽 显著的用粗线
            figure(fig)
            subplot(1,icell,ifgcell)
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
    allcellvL_fg = cell2mat(fgcellinfo(:,2));
    vecLength_slap(nsn,:) = mean(allcellvL_fg,'omitnan');% singlelap
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
sum(sum(cellnum_fg))



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