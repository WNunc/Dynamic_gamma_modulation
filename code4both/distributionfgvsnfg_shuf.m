
% note PhaseDiff2 have 7 columns. if not, need re-run the 
% doall_sgPhsPreCell_vs_fgPhsLocCell.m
% 
load('H:\neuralynx\sgamma result v3 std-1.5\matlab.mat')
load('H:\neuralynx\sgamma result v3 std-1.5\Phasediff2.mat')

fgcphsdiff = PhaseDiff2(PhaseDiff2(:,7)==1,1);
nfgcphsdiff = PhaseDiff2(PhaseDiff2(:,7)==0,1);
fgPhaseDiff2_sh = PhaseDiff2_sh((PhaseDiff2(:,7)==1),:);
nfgPhaseDiff2_sh = PhaseDiff2_sh((PhaseDiff2(:,7)==0),:);
% shuffle
X1 = fgPhaseDiff2_sh;
Xmid1 = median(X1);
%Xmid1 = mean(X1);
Xm1_sort = sort(Xmid1);
Xm1_sort95 = [Xm1_sort(26),Xm1_sort(975)];
Xm1_sort99 = [Xm1_sort(6),Xm1_sort(995)];

X0 = nfgPhaseDiff2_sh;
Xmid0 = median(X0);
%Xmid0 = mean(X0);
Xm0_sort = sort(Xmid0);
Xm0_sort95 = [Xm0_sort(26),Xm0_sort(975)];
Xm0_sort99 = [Xm0_sort(6),Xm0_sort(995)];

ffa = figure;

subplot(2,6,1)
hfg = histogram(fgcphsdiff,-2*pi:pi/8:2*pi,'Normalization','probability','FaceColor','#F24444');
hold on
hnfg = histogram(nfgcphsdiff,-2*pi:pi/8:2*pi,'Normalization','probability','FaceColor','#F2CA50');
H1_shuf = histogram(PhaseDiff2_sh((PhaseDiff2(:,7)==1),:),-2*pi:pi/8:2*pi,...
    'Normalization','probability','FaceColor',[0.2 0.2 0.2]);
H0_shuf = histogram(PhaseDiff2_sh((PhaseDiff2(:,7)==0),:),-2*pi:pi/8:2*pi,...
    'Normalization','probability','FaceColor',[0.2 0.2 0.2]);
hold off

title(['all laps cell'])%,...
%sprintf(['t-test p = ' num2str(p_t) '\nr-test p = ' num2str(p_r)]))

subplot(2,6,7)

hold on
bar(-2*pi+pi/16:pi/8:2*pi-pi/16,H1_shuf.Values,1,'FaceColor',[0.2 0.2 0.2],'FaceAlpha',0.6)
bar(-2*pi+pi/16:pi/8:2*pi-pi/16,-H0_shuf.Values,1,'FaceColor',[0.2 0.2 0.2],'FaceAlpha',0.6)
bar(-2*pi+pi/16:pi/8:2*pi-pi/16,hfg.Values,1,'FaceColor','#F24444','FaceAlpha',0.85)
bar(-2*pi+pi/16:pi/8:2*pi-pi/16,-hnfg.Values,1,'FaceColor','#F2CA50','FaceAlpha',0.85)
stem(median(PhaseDiff2(( PhaseDiff2(:,7)==1),1)),0.35,'Marker','none','Color','#05F2DB','LineWidth',2)
stem(Xm1_sort95,[0.35,0.35],'Marker','none','Color','#000000')
stem(median(PhaseDiff2((PhaseDiff2(:,7)==0),1)),-0.35,'Marker','none','Color','#05F2DB','LineWidth',2)
stem(Xm0_sort95,[-0.35,-0.35],'Marker','none','Color','#000000')
hold off

% 统计fgcell
[h,p_t1,ci,stats] = ttest(fgcphsdiff);
[p_r1,z] = circ_rtest(fgcphsdiff);
% 统计nfgcell
[h,p_t0,ci,stats] = ttest(nfgcphsdiff);
[p_r0,z] = circ_rtest(nfgcphsdiff);


title(['all Lap cell'],...
    sprintf(['t-test p = ' num2str(p_t1) '\nr-test p = ' num2str(p_r1) ...
    '\nt-test p = ' num2str(p_t0) '\nr-test p = ' num2str(p_r0)]))

xticks([-pi,0,pi])
xticklabels({'-π','0','π'})
yticks([-0.3,0,0.3])
yticklabels({'0.3','0','0.3'})
set(gca,'FontSize',12)
ylim([-0.3,0.3])
axis square
xlabel PhaseDiff
ylabel('Relative probability')
box on
ffb = figure('Position', [643 324 1032 342],'Color', [1 1 1]);
%%
for i = 1:5
    figure(ffa)
    subplot(2,6,i+1)
    H1 = histogram(PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==1),1),-2*pi:pi/8:2*pi,...
        'Normalization','probability','FaceColor','#F24444');
    hold on
    H1_shuf = histogram(PhaseDiff2_sh((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==1),:),-2*pi:pi/8:2*pi,...
        'Normalization','probability','FaceColor',[0.2 0.2 0.2]);
    H0 = histogram(PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==0),1),-2*pi:pi/8:2*pi,...
        'Normalization','probability','FaceColor','#F2CA50');
    H0_shuf = histogram(PhaseDiff2_sh((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==0),:),-2*pi:pi/8:2*pi,...
        'Normalization','probability','FaceColor',[0.2 0.2 0.2]);
    hold off
    
    % 统计fgcell
    [h,p_t1,ci,stats] = ttest(PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==1),1));
    [p_r1,z] = circ_rtest(PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==1),1));
    % 统计nfgcell
    [h,p_t0,ci,stats] = ttest(PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==0),1));
    [p_r0,z] = circ_rtest(PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==0),1));
    
    X1 = PhaseDiff2_sh((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==1),:);
    Xmid1 = median(X1);
    %Xmid1 = mean(X1);
    Xm1_sort = sort(Xmid1);
    Xm1_sort95(i,:) = [Xm1_sort(26),Xm1_sort(975)];
    Xm1_sort99(i,:) = [Xm1_sort(6),Xm1_sort(995)];
    
    X0 = PhaseDiff2_sh((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==0),:);
    Xmid0 = median(X0);
    %Xmid0 = mean(X0);
    Xm0_sort = sort(Xmid0);
    Xm0_sort95(i,:) = [Xm0_sort(26),Xm0_sort(975)];
    Xm0_sort99(i,:) = [Xm0_sort(6),Xm0_sort(995)];
    
    subplot(2,6,i+7)
    bar(-2*pi+pi/16:pi/8:2*pi-pi/16,H1_shuf.Values,1,'FaceColor',[0.2 0.2 0.2],'FaceAlpha',0.6)
    hold on
    bar(-2*pi+pi/16:pi/8:2*pi-pi/16,H1.Values,1,'FaceColor','#F24444','FaceAlpha',0.8)
    %H1 = histogram(PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==1),1),-2*pi:pi/8:2*pi,...
    %'Normalization','probability','FaceColor','blue');
    %     H1_shuf = histogram(PhaseDiff2_sh((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==1),:),-2*pi:pi/8:2*pi,...
    %         'Normalization','probability','FaceColor','#8A88C8');
    % stem(median(PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==1),1)),0.35,'Marker','none','Color','black','LineWidth',1.5)
    fgmid = median(PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==1),1));
    p1 = (numel(find(Xmid1<fgmid))+1)/1001
    stem(fgmid,0.35,'Marker','none','Color','#05F2DB','LineWidth',2)
    stem(Xm1_sort95(i,:),[0.35,0.35],'Marker','none','Color','#000000')
    bar(-2*pi+pi/16:pi/8:2*pi-pi/16,-H0_shuf.Values,1,'FaceColor',[0.2 0.2 0.2],'FaceAlpha',0.6)
    bar(-2*pi+pi/16:pi/8:2*pi-pi/16,-H0.Values,1,'FaceColor','#F2CA50','FaceAlpha',0.8)
    %     H0 = histogram(PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==0),1),-2*pi:pi/8:2*pi,...
    %         'Normalization','probability','FaceColor','yellow');
    %     H0_shuf = histogram(PhaseDiff2((PhaseDiff2_sh(:,3)==i & PhaseDiff2(:,7)==0),:),-2*pi:pi/8:2*pi,...
    %         'Normalization','probability','FaceColor','#73734B');
    % stem(median(PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==0),1)),-0.35,'Marker','none','Color','black','LineWidth',1.5)
    nfgmid = median(PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==0),1));
    p0 = (numel(find(Xmid0<nfgmid))+1)/1001
    
    stem(nfgmid,-0.35,'Marker','none','Color','#05F2DB','LineWidth',2)
    stem(Xm0_sort95(i,:),[-0.35,-0.35],'Marker','none','Color','#000000')
    hold off
    
    
    xlabel PhaseDiff
    ylabel('Relative probability')
    xticks([-pi,0,pi])
    xticklabels({'-π','0','π'})
    yticks([-0.35,0,0.35])
    yticklabels({'0.35','0','0.35'})
    axis square
    title(['Lap-' num2str(i)],...
    sprintf(['t-test p = ' num2str(p_t1) '\nr-test p = ' num2str(p_r1) ...
    '\nt-test p = ' num2str(p_t0) '\nr-test p = ' num2str(p_r0)]))
    set(gca,'FontSize',12)
    ylim([-0.35,0.35])
    
    % fraction
    fgcell_sg = PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==1 & PhaseDiff2(:,1)<0),:);
    fgcell_nsg = PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==1 & PhaseDiff2(:,1)>=0),:);
    fgcell_all = PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==1),:);
    
    nfgcell_sg = PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==0 & PhaseDiff2(:,1)<0),:);
    nfgcell_nsg = PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==0 & PhaseDiff2(:,1)>=0),:);
    nfgcell_all = PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,7)==0),:);
    figure(ffb);
    subplot(2,5,i)
    p = pie([size(fgcell_sg,1),size(fgcell_nsg,1)]);
    p(1).FaceColor = '#0583F2';
    p(3).FaceColor = '#9ABBD9';
    title(['fgcell number ' num2str(size(fgcell_all,1))])
    set(gca,'FontSize',11)
    subplot(2,5,i+5)
    p = pie([size(nfgcell_sg,1),size(nfgcell_nsg,1)]);
    p(1).FaceColor = '#0583F2';
    p(3).FaceColor = '#9ABBD9';
    title(['nfgcell number ' num2str(size(nfgcell_all,1))])
    set(gca,'FontSize',11)
    legend1 = legend('sgcell','nsgcell');
    set(legend1,...
        'Position',[0.896815037753413 0.544809019028459 0.0939922467336172 0.124269002536584],...
        'FontSize',11);
    % %     stem(Xm_sort99(i,:),[0.28,0.28],'Marker','none','Color','#707070')
    %     stem(Xm_sort95(i,:),[0.28,0.28],'Marker','none','Color','#000000')
    %     stem(Diff_mid(i,:),[0.28],'Marker','none','Color','red','LineWidth',1.5)
    %
    %     hold off
    % %     xlabel PhaseDiff
    % %     ylabel('Relative probability')
    % %     xticks([-pi,0,pi])
    % %     xticklabels({'-π','0','π'})
    %
    % %     axis square
    % %     [h,p_t,ci,stats] = ttest(PhaseDiff2_sh(PhaseDiff2(:,2)==i,1));
    % %     [p_r,z] = circ_rtest(PhaseDiff2_sh(PhaseDiff2(:,2)==i,1));
    % %     title(['Lap-' num2str(i)],...
    % %     sprintf(['t-test p = ' num2str(p_t) '\nr-test p = ' num2str(p_r)]))
    % %     set(gca,'FontSize',12)
    % %     ylim([0,0.32])
    
end
%% 前期后期
ffc1 = figure;
ffc2 = figure;
II = [1,2;4,5];
gp = {'early','late'};

for i = 1:2%前期后期
    Xm1_sort95 = [];
    Xm1_sort99 = [];
    figure(ffc1)
    e = II(i,:);
    fgcase = (PhaseDiff2(:,3)==e(1)|PhaseDiff2(:,3)==e(2))& PhaseDiff2(:,7)==1;
    nfgcase = (PhaseDiff2(:,3)==e(1)|PhaseDiff2(:,3)==e(2))& PhaseDiff2(:,7)==0;
    subplot(1,2,i)
    H1 = histogram(PhaseDiff2((fgcase),1),-2*pi:pi/8:2*pi,...
        'Normalization','probability','FaceColor','#F24444');
    hold on
    H1_shuf = histogram(PhaseDiff2_sh((fgcase),:),-2*pi:pi/8:2*pi,...
        'Normalization','probability','FaceColor',[0.2 0.2 0.2]);
    H0 = histogram(PhaseDiff2((nfgcase),1),-2*pi:pi/8:2*pi,...
        'Normalization','probability','FaceColor','#F2CA50');
    H0_shuf = histogram(PhaseDiff2_sh((nfgcase),:),-2*pi:pi/8:2*pi,...
        'Normalization','probability','FaceColor',[0.2 0.2 0.2]);
    hold off
    
    % 统计fgcell
    [h,p_t1,ci,stats] = ttest(PhaseDiff2((fgcase),1));
    [p_r1,z] = circ_rtest(PhaseDiff2((fgcase),1));
    % 统计nfgcell
    [h,p_t0,ci,stats] = ttest(PhaseDiff2((nfgcase),1));
    [p_r0,z] = circ_rtest(PhaseDiff2((nfgcase),1));
    
    X1 = PhaseDiff2_sh((fgcase),:);
    Xmid1 = median(X1);
    %Xmid1 = mean(X1);
    Xm1_sort = sort(Xmid1);
    Xm1_sort95(i,:) = [Xm1_sort(26),Xm1_sort(975)];
    Xm1_sort99(i,:) = [Xm1_sort(6),Xm1_sort(995)];
    
    X0 = PhaseDiff2_sh((nfgcase),:);
    Xmid0 = median(X0);
    %Xmid0 = mean(X0);
    Xm0_sort = sort(Xmid0);
    Xm0_sort95(i,:) = [Xm0_sort(26),Xm0_sort(975)];
    Xm0_sort99(i,:) = [Xm0_sort(6),Xm0_sort(995)];
    figure(ffc2)
    subplot(1,2,i)
    bar(-2*pi+pi/16:pi/8:2*pi-pi/16,H1_shuf.Values,1,'FaceColor',[0.2 0.2 0.2],'FaceAlpha',0.6)
    hold on
    bar(-2*pi+pi/16:pi/8:2*pi-pi/16,H1.Values,1,'FaceColor','#F24444','FaceAlpha',0.8)
    %H1 = histogram(PhaseDiff2((fgcase),1),-2*pi:pi/8:2*pi,...
    %'Normalization','probability','FaceColor','blue');
    %     H1_shuf = histogram(PhaseDiff2_sh((fgcase),:),-2*pi:pi/8:2*pi,...
    %         'Normalization','probability','FaceColor','#8A88C8');
    % stem(median(PhaseDiff2((fgcase),1)),0.35,'Marker','none','Color','black','LineWidth',1.5)
    
    stem(median(PhaseDiff2((fgcase),1)),0.35,'Marker','none','Color','#05F2DB','LineWidth',2)
    stem(Xm1_sort95(i,:),[0.35,0.35],'Marker','none','Color','#000000')
    bar(-2*pi+pi/16:pi/8:2*pi-pi/16,-H0_shuf.Values,1,'FaceColor',[0.2 0.2 0.2],'FaceAlpha',0.6)
    bar(-2*pi+pi/16:pi/8:2*pi-pi/16,-H0.Values,1,'FaceColor','#F2CA50','FaceAlpha',0.8)
    %     H0 = histogram(PhaseDiff2((nfgcase),1),-2*pi:pi/8:2*pi,...
    %         'Normalization','probability','FaceColor','yellow');
    %     H0_shuf = histogram(PhaseDiff2((PhaseDiff2_sh(:,3)==i & PhaseDiff2(:,7)==0),:),-2*pi:pi/8:2*pi,...
    %         'Normalization','probability','FaceColor','#73734B');
    % stem(median(PhaseDiff2((nfgcase),1)),-0.35,'Marker','none','Color','black','LineWidth',1.5)
    p0 = numel(find(Xmid0<median(PhaseDiff2((nfgcase),1))))+1/1001
    stem(median(PhaseDiff2((nfgcase),1)),-0.35,'Marker','none','Color','#05F2DB','LineWidth',2)
    stem(Xm0_sort95(i,:),[-0.35,-0.35],'Marker','none','Color','#000000')
    hold off
    
    
    xlabel PhaseDiff
    ylabel('Relative probability')
    xticks([-pi,0,pi])
    xticklabels({'-π','0','π'})
    yticks([-0.35,0,0.35])
    yticklabels({'0.35','0','0.35'})
    axis square
    title([gp{i}],...
    sprintf(['t-test p = ' num2str(p_t1) '\nr-test p = ' num2str(p_r1) ...
    '\nt-test p = ' num2str(p_t0) '\nr-test p = ' num2str(p_r0)]))
    set(gca,'FontSize',12)
    ylim([-0.35,0.35])
    
%     % fraction
%     fgcell_sg = PhaseDiff2((fgcase & PhaseDiff2(:,1)<0),:);
%     fgcell_nsg = PhaseDiff2((fgcase & PhaseDiff2(:,1)>=0),:);
%     fgcell_all = PhaseDiff2((fgcase),:);
%     
%     nfgcell_sg = PhaseDiff2((nfgcase & PhaseDiff2(:,1)<0),:);
%     nfgcell_nsg = PhaseDiff2((nfgcase & PhaseDiff2(:,1)>=0),:);
%     nfgcell_all = PhaseDiff2((nfgcase),:);
%     figure(ffb);
%     subplot(2,5,i)
%     p = pie([size(fgcell_sg,1),size(fgcell_nsg,1)]);
%     p(1).FaceColor = '#0583F2';
%     p(3).FaceColor = '#9ABBD9';
%     
%     subplot(2,5,i+5)
%     p = pie([size(nfgcell_sg,1),size(nfgcell_nsg,1)]);
%     p(1).FaceColor = '#0583F2';
%     p(3).FaceColor = '#9ABBD9';
    % %     stem(Xm_sort99(i,:),[0.28,0.28],'Marker','none','Color','#707070')
    %     stem(Xm_sort95(i,:),[0.28,0.28],'Marker','none','Color','#000000')
    %     stem(Diff_mid(i,:),[0.28],'Marker','none','Color','red','LineWidth',1.5)
    %
    %     hold off
    % %     xlabel PhaseDiff
    % %     ylabel('Relative probability')
    % %     xticks([-pi,0,pi])
    % %     xticklabels({'-π','0','π'})
    %
    % %     axis square
    % %     [h,p_t,ci,stats] = ttest(PhaseDiff2_sh(PhaseDiff2(:,2)==i,1));
    % %     [p_r,z] = circ_rtest(PhaseDiff2_sh(PhaseDiff2(:,2)==i,1));
    % %     title(['Lap-' num2str(i)],...
    % %     sprintf(['t-test p = ' num2str(p_t) '\nr-test p = ' num2str(p_r)]))
    % %     set(gca,'FontSize',12)
    % %     ylim([0,0.32])
    
end
