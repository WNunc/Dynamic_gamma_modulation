% 后一个slow gamma周期中的spike相位减去前一个周期的相位，统计这一差值的分布
% 在doall
clear
close all
n_std = 1.5;
directories_allData_v0_allgood
resuletFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
Dx = {'-cw','-ccw'};
case1 = '-ontrack';
ncells = 3;nspkikes = 5;
case2 = sprintf('_%ucell_%uspk',ncells, nspkikes);
midmod = '_cycmid';
nsn = 0;
tcolor = {'black','red'};
peak0 = 1; % firing rate threshold to remove place cells

fileinput1 = 'sgamma_phase.mat';% sgamma phase
fileinput2 = ['sgamma_dominant_thetacyc_IND_std' num2str(n_std) '_v2.mat'];% sg dominant theta cycle ind
fileinput3 = [];% good theta cycle in each lap
fileinput4 = [];% all theta cycle in each lap
inputFolder = ['H:\neuralynx\sgamma result v3 std-' num2str(n_std) '\'];
fileoutput = 'sgphase_difference_v3_3spk.mat';
nlap = 5;
PhaseDiff1 = [];PhaseDiff2 = [];% 汇总
for ns = [1,2,4:6,8,10:13,15,16]%1:isession%
    path_ns = path{ns};
    disp(path_ns)
    cd(path_ns)
    goodphase_ns = Seq_cutPhase{ns,1};%***
    case3 = num2str(goodphase_ns);% directories_allData_v0_allgood 俩方向相位一样
    inFolder1 = [resuletFolder path_ns(13:end)];
    inFolder2 = [inputFolder path_ns(13:end)];
    outFolder = inFolder2;
    % 重要数据 phasediff1 按theta cyc 统计
    % column1 = phase diff
    % column2 = lap ID
    % column3 = theta cycle ID
    % column4 = sequence or not
    phasediff1 = [];
    % 重要数据 phasediff2 按cell统计
    % column1 = phase diff
    % column2 = cell ID
    % column3 = lap ID
    % column4 = n-spk pair
    % column5 = theta cycle ID
    % column6 = sequence or not
    phasediff2 = [];
    
    phsdf_n1 = 0;phsdf_n2 = 0;
    % 导入 sgamma phase 和 sg dominant theta ind
    load([inFolder2,fileinput1])
    load([inFolder1,fileinput2])
    icell = length(SPK{1});
    % 导入theta cycle
    fileinput4 = ['data_theta_seq_info_AllLap',case1,case2,case3,midmod,'_v5-both.mat'];
    load([inFolder1,fileinput4])
    for nl = 1:nlap
        fileinput3 = ['data_theta_seq_info',case1,'_lap',num2str(nl),case2,case3,midmod,'_v5-both.mat'];
        load([inFolder1,fileinput3],'ThetaGood')
        thcyc_sl = ind_th_sl{nl};
        nthc_sl = length(thcyc_sl);
        thcyc_al = ind_th_al{nl};
        nthc_al = length(thcyc_al);
        %% 单圈
        if ~isempty(thcyc_sl)
            phsd_acorssthcyc = [];
            for nth = 1:nthc_sl
                phsd_t = [];
                indth = thcyc_sl(nth);
                if indth > size(theta_INFO{nl},1)
                    continue
                end
                thcyc_onset = theta_INFO{nl}{indth,3}(1);
                thcyc_offset = theta_INFO{nl}{indth,3}(2);
                thseq = theta_INFO{nl}{indth,8};
                thseq_t = thcyc_onset:0.005:thcyc_offset;
                thseq_s = 0.5*2*pi/90:2*pi/90:2*pi-0.5*2*pi/90;
                for nc = 1:icell
                    spkind = find(SPK{nl}{nc}>thcyc_onset & SPK{nl}{nc}<thcyc_offset);
                    if isempty(spkind)
                        continue
                    else
                        if ~ismember(theta_INFO{nl}{indth,1},[ThetaGood{:,1}])
                            c = 0; % 不在满足条件的sequence
                            figname = sprintf('B Lap-%u thetaCyc-%u cell-%u',nl, indth, nc);
                        else
                            c = 1; % 在满足条件的sequence
                            figname = sprintf('A Lap-%u thetaCyc-%u cell-%u',nl, indth, nc);
                        end
                        t = SPK{nl}{nc}(spkind);
                        p = SPKPhase{nl}{nc}(spkind);
                        itt = LFP{nl,1}>thcyc_onset&LFP{nl,1}<thcyc_offset;
                        SGTT = LFP{nl,1}(itt);
                        SGwave = LFP{nl,2}(CSCnum(nc),itt)+LFP{nl,3}(CSCnum(nc),itt);
                        SGphase = LFP{nl,4}(CSCnum(nc),itt);
                        t_sgcyc = findpeaks(SGphase);
                        t_sgcyc = t_sgcyc.loc;
                        t_sgcyc = SGTT(t_sgcyc(1:end-1));
                        [sgcycind] = get_sgphase(SGTT, SGphase,[thcyc_onset,thcyc_offset],SPK{nl}{nc}(spkind));
                        sgphs_cyc0 = [p*pi/180,sgcycind];
                        sgphs_cyc = [];
                        trigger = length(unique(sgcycind));% spike 分布在sg的周期数
                        phsd = [];
                        if trigger>=4 % spike个数限制2or3
                            sgcIND = unique(sgcycind);
                            % 将一个sg周期中的spike相位平均值作为这个周期的相位
                            [sgphs_cyc] = meanSameCyclePhase(sgphs_cyc0);
                            for ndf = 1:trigger-1
                                %                                 if (sgphs_cyc(ndf+1,1) - sgphs_cyc(ndf,1))>pi
                                %                                     phsd(ndf) = (sgphs_cyc(ndf+1,1) - sgphs_cyc(ndf,1))-2*pi;
                                %                                 elseif (sgphs_cyc(ndf+1,1) - sgphs_cyc(ndf,1))<-pi
                                %                                     phsd(ndf) = (sgphs_cyc(ndf+1,1) - sgphs_cyc(ndf,1))+2*pi;
                                %                                 else
                                %                                     phsd(ndf) = (sgphs_cyc(ndf+1,1) - sgphs_cyc(ndf,1));
                                %                                 end
                                phsd(ndf) = sgphs_cyc(ndf+1,1) - sgphs_cyc(ndf,1);
                            end
                            if trigger>=3
                                % pause(0.1)
                            end
                            phsdf_n2 = phsdf_n2 + 1;
                            phasediff2(phsdf_n2,1) = mean(phsd);%circ_mean(phsd,[],2);%
                            phasediff2(phsdf_n2,2) = nc;
                            phasediff2(phsdf_n2,3) = nl;
                            phasediff2(phsdf_n2,4) = trigger-1;
                            phasediff2(phsdf_n2,5) = indth;
                            phasediff2(phsdf_n2,6) = c;
                        end
                        aa = mean(phsd);
                        % aa = circ_mean(phsd,[],2);
                    end
                    
                    phsd_t = [phsd_t,aa];%,phsd
                    phsd_t = phsd_t(~isnan(phsd_t));
                end
                if isempty(phsd_t)
                    continue
                end
                phsdf_n1 = phsdf_n1 + 1;
                phasediff1(phsdf_n1,1) = circ_mean(phsd_t,[],2);% (sgphs_cyc(ndf+1,1) - sgphs_cyc(ndf,1))-sgphs_cyc(ndf,1)+pi;%
                phasediff1(phsdf_n1,2) = nl;
                phasediff1(phsdf_n1,3) = indth;
                phasediff1(phsdf_n1,4) = c;
            end
        end
        %
        %         if ~isempty(thcyc_al)
        %             for nth = 1:nthc_al
        %                 indth = thcyc_al(nth);
        %                 thcyc_onset = ThetaGood{indth,3}(1);
        %                 thcyc_offset = ThetaGood{indth,3}(2);
        %                 for nc = 1:icell
        
        %                 end
        %             end
        %         end
    end

    PhaseDiff1 = [PhaseDiff1;phasediff1];
    PhaseDiff2 = [PhaseDiff2;phasediff2];
%     save([outFolder fileoutput],'phasediff1','phasediff2');
end

figure;histogram(PhaseDiff1(:,1),-2*pi:pi/8:2*pi,'Normalization','probability');
[h,p_t,ci,stats] = ttest(PhaseDiff1(:,1));
[p_r,z] = circ_rtest(PhaseDiff1(:,1));
xlabel PhaseDiff
ylabel('Relative probability')
title(['all laps thetacyc'],...
    sprintf(['t-test p = ' num2str(p_t) '\nr-test p = ' num2str(p_r)]))

xticks([-pi,0,pi])
xticklabels({'-π','0','π'})
set(gca,'FontSize',12)
axis square

figure;histogram(PhaseDiff2(:,1),-2*pi:pi/8:2*pi,'Normalization','probability');
[h,p_t,ci,stats] = ttest(PhaseDiff2(:,1));
[p_r,z] = circ_rtest(PhaseDiff2(:,1));
xlabel PhaseDiff
ylabel('Relative probability')
title(['all laps cell'],...
    sprintf(['t-test p = ' num2str(p_t) '\nr-test p = ' num2str(p_r)]))
xticks([-pi,0,pi])
xticklabels({'-π','0','π'})
set(gca,'FontSize',12)
axis square
% temp = PhaseDiff(:,1);
% temp(temp>pi) = temp(temp>pi)-2*pi;
% temp(temp<-pi) = temp(temp<-pi)+2*pi;
% figure;histogram(temp,-2*pi:pi/8:2*pi,'Normalization','pdf');
% xlabel PhaseDiff
% ylabel count
% [h,p,ci,stats] = ttest(temp)
% title(['all laps '],['t-test p = ' num2str(p)])
figure
for i = 1:5
    subplot(1,5,i)
    histogram(PhaseDiff1(PhaseDiff1(:,2)==i,1),-2*pi:pi/8:2*pi,'Normalization','probability');
    xlabel PhaseDiff
    ylabel('Relative probability')
    xticks([-pi,0,pi])
    xticklabels({'-π','0','π'})
    
    axis square
    [h,p_t,ci,stats] = ttest(PhaseDiff1(PhaseDiff1(:,2)==i,1));
    [p_r,z] = circ_rtest(PhaseDiff1(PhaseDiff1(:,2)==i,1));
    title(['Lap-' num2str(i)],...
    sprintf(['t-test p = ' num2str(p_t) '\nr-test p = ' num2str(p_r)]))
    set(gca,'FontSize',12)
    ylim([0,0.32])
end

figure
for i = 1:5
    subplot(1,5,i)
    histogram(PhaseDiff2(PhaseDiff2(:,3)==i,1),-2*pi:pi/8:2*pi,'Normalization','probability');
    xlabel PhaseDiff
    ylabel('Relative probability')
    xticks([-pi,0,pi])
    xticklabels({'-π','0','π'})
    
    axis square
    [h,p_t,ci,stats] = ttest(PhaseDiff2(PhaseDiff2(:,3)==i,1));
    [p_r,z] = circ_rtest(PhaseDiff2(PhaseDiff2(:,3)==i,1));
    title(['Lap-' num2str(i)],...
    sprintf(['t-test p = ' num2str(p_t) '\nr-test p = ' num2str(p_r)]))
    set(gca,'FontSize',12)
    ylim([0,0.32])
end

% 将一个sg周期中的spike相位平均值作为这个周期的相位
function [sgphs_cyc] = meanSameCyclePhase(sgphs_cyc0)
cycind = unique(sgphs_cyc0(:,2));
for c = 1:length(cycind)
    r = sgphs_cyc0(:,2) == cycind(c);
    % sgphs_cyc(c,1) = circ_mean(sgphs_cyc0(r,1),[],1);
    sgphs_cyc(c,1) = mean(sgphs_cyc0(r,1));
    sgphs_cyc(c,2) = cycind(c);
end
end