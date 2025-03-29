% investigating the relationship between the sequence slope and the phase shift
% step 2 linear regression between the phasediff and slope or
% weightcorrelation
% load all sequence slope
% load all sequence weightcorr
% load all phase diff
% phasediff1 按thetacyc统计
clear
close all
n_std = 1.5;
directories_allData_v0_allgood
resuletFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
Dx = {'-cw','-ccw'};
case1 = {'-ontrack','-ontrack_exfg','-ontrack_dsfg'};
ncells = 3;nspkikes = 5;
case2 = sprintf('_%ucell_%uspk',ncells, nspkikes);
midmod = '_cycmid';

fileinput2 = [];% good theta cycle in each lap
fileinput3 = [];% all theta cycle in each lap

lockat = {'firstlap','alllap','f2lap'};
L = 2;
PHSDvs = [];
for ns = [1,2,4:6,8,10:13,15,16]%1:isession%
    path_ns = path{ns};
    disp(path_ns);
    outFolder = [resuletFolder path_ns(13:end)];
    cd(outFolder);
    mkdir('eachSequence')
    pause(0.1)
    goodphase_ns = Seq_cutPhase{ns,1};%***CW和CCW分割相位一致，用1和2都一样
    case3 = num2str(goodphase_ns);
    % all theta cycle in each lap both Direction
    TGa = load([outFolder,'data_theta_seq_info','_AllLap',case1{1},case2,case3,midmod,'_v5-both.mat'],...
        'theta_INFO');
    TGa = TGa.theta_INFO;
    
    % phase difference
    pd_Folder = ['H:\neuralynx\sgamma result v3 std-' num2str(n_std) '\' path_ns(13:end) '\'];
    pd_file = 'sgphase_difference_v3_3spk.mat';
    load([pd_Folder pd_file]);
    % weight correlation
    WC = load([outFolder 'weight_corr_eachseq.mat']);
    WC = WC.wc1;
    % all sequence slope
    load([outFolder 'Slope_eachgoodseq.mat'])
    
    sn = 0;
    phsvs2 = [];
    for nl = 1:5
        Da = [];Dg = [];
        disp(['======== plot lap ' num2str(nl) ' ========'])
        % 把thetagood和thetainfo还有WC对应一下
        % load good theta cycle in each lap both Direction
        TG = load([outFolder,'data_theta_seq_info',case1{1},'_lap',num2str(nl),case2,case3,midmod,'_v5-both.mat'],'ThetaGood');
        TG = TG.ThetaGood;
        indgood = [TG{:,1}];
        Lg = length(indgood);
        Sg = find(diff(indgood)<0);% theta good 中CW和CCW的分界
        if ~isempty(Sg)
            Dg(1:Sg) = 1;
            Dg(Sg+1:Lg) = 2;
        end
        
        indall = [TGa{nl}{:,1}];
        La = length(indall);
        Sa = find(diff(indall)<0);% theta info 中CW和CCW的分界
        Da(1:Sa) = 1;
        Da(Sa+1:La) = 2;
        
        Num_cyc = size(TGa{nl},1);
        for t = 1:Num_cyc
            ffname = ['all lap ' num2str(nl) ' cyc-' num2str(TGa{nl}{t,1}) '-' num2str(t)];
            inds = find(phasediff1(:,2)==nl & phasediff1(:,3)==t);
            
            % 修正了phasediff1，把存在phasediff是负的cell的
            % thetacycle，修正为cell中phasediff的最小值
            Z = phasediff2((phasediff2(:,5) == t),1);
            X = sum(sum(Z<0));
            if X>0 && length(Z) >1
                a = min(Z);
                phasediff1(phasediff1(:,3)==t,1) = a;
            end
            if ~isempty(inds)
                sn = sn+1;
                phsvs2(sn,1) = phasediff1(inds,1);
                phsvs2(sn,2) = slope_all{nl}{t,5};
                if ~isempty(Dg) % 如果不为空，说明thetagood可能有重复的cycid
                    wcid = find(indgood==indall(t) & Dg==Da(t));
                else % 否则是没有重复的
                    wcid = find(indgood==indall(t));
                end
                
                if ~isempty(wcid) % 如果不为空则是goodseq
                    phsvs2(sn,3) = WC{nl}(wcid);
                else %否则是个不满足要求的seq
                    % 计算这个sequence的WC
                    all = TGa{nl}{t,8};% seq
                    loc = TGa{nl}{t,7};% real loc bin
                    tbins = size(all,2);
                    TPsweep = 7;
                    SPsweep = 10;
                    TPmid = round(tbins*0.5);
                    SPmid = loc(TPmid);
                    tp1 = [TPmid-TPsweep:TPmid];
                    tp2 = [TPmid+1:TPmid+1+TPsweep];
                    sp1 = [SPmid-SPsweep:SPmid];
                    sp2 = [SPmid+1:SPmid+1+SPsweep];
                    [wca,~]= nanweightcorr(all([sp1,sp2],[tp1,tp2]),TPsweep,SPsweep);
                    phsvs2(sn,3) = wca;
                end
                phsvs2(sn,4) = nl;
                phsvs2(sn,5) = t;
                phsvs2(sn,6) = phasediff1(inds,4);%good seq or not
            end
        end
    end
    %     save([outFolder 'Slope_eachgoodseq.mat'],'slope_good','slope_all','phsvsslope')
    PHSDvs = [PHSDvs;phsvs2];
end

%%
load('H:\neuralynx\sgamma result v3 std-1.5\matlab.mat', 'PhaseDiff1_sh')
close all
phbin = -pi:pi/10:pi;
wcbin = -1:1/10:1;
slbin = -40:2:40;
nphsDiff = size(PHSDvs,1);
% weight corr
x = PHSDvs(:,1);
y = PHSDvs(:,3);
% x = PHSDvs(PHSDvs(:,6)==1,1);
% y = PHSDvs(PHSDvs(:,6)==1,3);
figure('Units','centimeters','Position',[8 8 38 6],'Renderer','painters')
[p1_real,P1s_real,b,bint,r,rint,stats] = scatimg_XY(x,y,phbin,wcbin,'PhaseDiff','WeightCorr',1);
caxis([0,0.01])
Qratio.Q1 = length(find(x<0&y>=0))/nphsDiff;
Qratio.Q2 = length(find(x>0&y>=0))/nphsDiff;
Qratio.Q3 = length(find(x<0&y<0))/nphsDiff;
Qratio.Q4 = length(find(x>0&y<0))/nphsDiff;


Q1sh = [];Q2sh = [];Q3sh = [];Q4sh = [];
Psh1 = zeros(20,20);P1ssh1 = zeros(20,20); %P1_smooth_shuffle
for nsh = 1:1000
x_sh = PhaseDiff1_sh(:,nsh);
% x_sh = PhaseDiff1_sh(PHSDvs(:,6)==1,1);
[p1,P1s,b,bint,r,rint,stats] = scatimg_XY(x_sh,y,phbin,wcbin,'PhaseDiff','WeightCorr',0);
Q1sh(nsh) = length(find(x_sh<0&y>=0))/nphsDiff;
Q2sh(nsh) = length(find(x_sh>0&y>=0))/nphsDiff;
Q3sh(nsh) = length(find(x_sh<0&y<0))/nphsDiff;
Q4sh(nsh) = length(find(x_sh>0&y<0))/nphsDiff;
p1 = p1/nphsDiff;
Psh1 = Psh1+p1;
P1ssh1 = P1ssh1+P1s;
end
Psh1 = Psh1/1000;
P1ssh1 = P1ssh1/1000;
subplot(1,4,3)
imagesc(phbin,wcbin,P1ssh1)
axis xy; axis square;
xlabel('PhaseDiff','Fontsize',11)
set(gca,'XTick',-pi:pi:pi);
set(gca,'XTickLabel',{'-π','0','π'},'Fontsize',11)
ylabel('WeightCorr','Fontsize',11)
colormap(jet(128));
caxis([0,0.01])
c = colorbar;
c.Label.String = 'probability';

% figure('Units','centimeters','Position',[38 8 8 6])
subplot(1,4,4)
plot_hist_shandreal(Q1sh,Qratio.Q1,0.1:0.04:0.9)
p = (numel(find(Q1sh>Qratio.Q1))+1)/1001
axis square
title('alllap')
xlabel('Probability in Q1','FontSize',11)
ylabel('Count','FontSize',11)
set(gca,'FontSize',11)
% save('H:\neuralynx\sgamma result v3 std-1.5\alllap_phsvswc.mat','p1_real','P1s_real','Qratio',...
%     'Q1sh','Q2sh','Q3sh','Q4sh','Psh1','P1ssh1')
% saveas(gcf,'H:\neuralynx\sgamma result v3 std-1.5\alllap_phsvswc.png');
% saveas(gcf,'H:\neuralynx\sgamma result v3 std-1.5\alllap_phsvswc','epsc');

% 前期
indearly = (PHSDvs(:,4)==1 | PHSDvs(:,4)==2);%&PHSDvs(:,6)==1;
x = PHSDvs(indearly,1);
y = PHSDvs(indearly,3);
nphsDiff_e = length(x);
figure('Units','centimeters','Position',[8 8 38 6],'Renderer','painters')
[p1el_real,P1sel_real,b,bint,r,rint,stats] = scatimg_XY(x,y,phbin,wcbin,'PhaseDiff','WeightCorr',1);
caxis([0,0.01])
Qratio.Q1 = length(find(x<0&y>=0))/nphsDiff_e;
Qratio.Q2 = length(find(x>0&y>=0))/nphsDiff_e;
Qratio.Q3 = length(find(x<0&y<0))/nphsDiff_e;
Qratio.Q4 = length(find(x>0&y<0))/nphsDiff_e;


Q1sh = [];Q2sh = [];Q3sh = [];Q4sh = [];
Psh1 = zeros(20,20);P1ssh1 = zeros(20,20); %P1_smooth_shuffle
for nsh = 1:1000
x_sh = PhaseDiff1_sh(indearly,nsh);
[p1,P1s,b,bint,r,rint,stats] = scatimg_XY(x_sh,y,phbin,wcbin,'PhaseDiff','WeightCorr',0);
Q1sh(nsh) = length(find(x_sh<0&y>=0))/nphsDiff_e;
Q2sh(nsh) = length(find(x_sh>0&y>=0))/nphsDiff_e;
Q3sh(nsh) = length(find(x_sh<0&y<0))/nphsDiff_e;
Q4sh(nsh) = length(find(x_sh>0&y<0))/nphsDiff_e;
p1 = p1/nphsDiff_e;
Psh1 = Psh1+p1;
P1ssh1 = P1ssh1+P1s;
end
Psh1 = Psh1/1000;
P1ssh1 = P1ssh1/1000;
subplot(1,4,3)
imagesc(phbin,wcbin,P1ssh1)
axis xy; axis square;
xlabel('PhaseDiff','Fontsize',11)
set(gca,'XTick',-pi:pi:pi);
set(gca,'XTickLabel',{'-π','0','π'},'Fontsize',11)
ylabel('WeightCorr','Fontsize',11)
colormap(jet(128));
caxis([0,0.01])
c = colorbar;
c.Label.String = 'probability';

% figure('Units','centimeters','Position',[38 8 8 6])
subplot(1,4,4)
plot_hist_shandreal(Q1sh,Qratio.Q1,0.1:0.04:0.9)
axis square
title('early')
xlabel('Probability in Q1','FontSize',11)
ylabel('Count','FontSize',11)
set(gca,'FontSize',11)
% save('H:\neuralynx\sgamma result v3 std-1.5\earlystage_phsvswc.mat','p1el_real','P1sel_real','Qratio',...
%     'Q1sh','Q2sh','Q3sh','Q4sh','Psh1','P1ssh1')
% saveas(gcf,'H:\neuralynx\sgamma result v3 std-1.5\earlystage_phsvswc.png');
% saveas(gcf,'H:\neuralynx\sgamma result v3 std-1.5\earlystage_phsvswc','epsc');

% 后期
indlate = (PHSDvs(:,4)==4 | PHSDvs(:,4)==5);%&PHSDvs(:,6)==1;
x = PHSDvs(indlate,1);
y = PHSDvs(indlate,3);
nphsDiff_l = length(x);
figure('Units','centimeters','Position',[8 8 38 6],'Renderer','painters')
[p1lt_real,P1slt_real,b,bint,r,rint,stats] = scatimg_XY(x,y,phbin,wcbin,'PhaseDiff','WeightCorr',1);
caxis([0,0.01])
Qratio.Q1 = length(find(x<0&y>=0))/nphsDiff_l;
Qratio.Q2 = length(find(x>0&y>=0))/nphsDiff_l;
Qratio.Q3 = length(find(x<0&y<0))/nphsDiff_l;
Qratio.Q4 = length(find(x>0&y<0))/nphsDiff_l;


Q1sh = [];Q2sh = [];Q3sh = [];Q4sh = [];
Psh1 = zeros(20,20);P1ssh1 = zeros(20,20); %P1_smooth_shuffle
for nsh = 1:1000
x_sh = PhaseDiff1_sh(indlate,nsh);
[p1,P1s,b,bint,r,rint,stats] = scatimg_XY(x_sh,y,phbin,wcbin,'PhaseDiff','WeightCorr',0);
Q1sh(nsh) = length(find(x_sh<0&y>=0))/nphsDiff_l;
Q2sh(nsh) = length(find(x_sh>0&y>=0))/nphsDiff_l;
Q3sh(nsh) = length(find(x_sh<0&y<0))/nphsDiff_l;
Q4sh(nsh) = length(find(x_sh>0&y<0))/nphsDiff_l;
p1 = p1/nphsDiff_l;
Psh1 = Psh1+p1;
P1ssh1 = P1ssh1+P1s;
end
Psh1 = Psh1/1000;
P1ssh1 = P1ssh1/1000;
subplot(1,4,3)
imagesc(phbin,wcbin,P1ssh1)
axis xy; axis square;
xlabel('PhaseDiff','Fontsize',11)
set(gca,'XTick',-pi:pi:pi);
set(gca,'XTickLabel',{'-π','0','π'},'Fontsize',11)
ylabel('WeightCorr','Fontsize',11)
colormap(jet(128));
caxis([0,0.01])
c = colorbar;
c.Label.String = 'probability';

% figure('Units','centimeters','Position',[38 8 8 6])
subplot(1,4,4)
plot_hist_shandreal(Q1sh,Qratio.Q1,0.1:0.04:0.9)
axis square
title('late')
xlabel('Probability in Q1','FontSize',11)
ylabel('Count','FontSize',11)
set(gca,'FontSize',11)
% save('H:\neuralynx\sgamma result v3 std-1.5\latestage_phsvswc.mat','p1lt_real','P1slt_real','Qratio',...
%     'Q1sh','Q2sh','Q3sh','Q4sh','Psh1','P1ssh1')
% saveas(gcf,'H:\neuralynx\sgamma result v3 std-1.5\latestage_phsvswc.png');
% saveas(gcf,'H:\neuralynx\sgamma result v3 std-1.5\latestage_phsvswc','epsc');
% % 斜率
% x = PHSDvs(:,1);
% y = PHSDvs(:,2);
% [b,bint,r,rint,stats] = scatimg_XY(x,y,phbin,slbin,'PhaseDiff','slope',0);
% caxis([0,0.01])
% fgn = 'all lap phs vs wc all late 2';
% saveas(gcf,fgn,'epsc')
% saveas(gcf,[fgn,'.png'])
% function plot_hist_shandreal(Qsh,Qreal,bin)
% histogram(Qsh,bin)
% CI95 = prctile(Qsh,95);
% xline(CI95,'k-')
% xline(Qreal,'r-')
% end


% function [b,bint,r,rint,stats] = scatimg_XY(x,y,xbin,ybin,xl,yl)
% % 画图散点图和拟合直线以及概率热图
% % xbin = -pi:pi/20:pi;
% % ybin = -1:1/20:1;
% 
% X = [ones(length(x),1),x];
% 
% % linear regress
% [b,bint,r,rint,stats] = regress(y,X);
% figure
% scatter(x,y)
% % 画散点图和拟合直线
% xFit = -pi:0.1*pi:pi;
% yFit = b(1) + b(2)*xFit;
% hold on
% plot(xFit,yFit,'k--','LineWidth',1.5)
% hold off
% title(strcat(['R^2 = ',num2str(stats(1)),'    p = ',num2str(stats(3))]),'FontSize',10)
% axis xy; axis square;
% xlabel(xl,'Fontsize',18)
% set(gca,'XTick',-pi:pi:pi);
% set(gca,'XTickLabel',{'-π','0','π'},'Fontsize',15)
% ylabel(yl,'Fontsize',18)
% ylim([-1,1])
% % 画平滑后的概率热图
% p1=hist2d(x,y,xbin,ybin);
% P1=p1/sum(sum(p1));
% cyc_1 = smooth2a(normz(P1,0),5,5);
% figure
% imagesc(xbin,ybin,cyc_1)
% axis xy; axis square;
% xlabel(xl,'Fontsize',18)
% set(gca,'XTick',-pi:pi:pi);
% set(gca,'XTickLabel',{'-π','0','π'},'Fontsize',15)
% ylabel(yl,'Fontsize',18)
% colormap(jet(128));
% end

