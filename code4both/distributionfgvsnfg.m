


fgcphsdiff = PhaseDiff2(PhaseDiff2(:,end)==1,1);
nfgcphsdiff = PhaseDiff2(PhaseDiff2(:,end)==0,1);
figure
subplot(1,6,1)
histogram(fgcphsdiff,-2*pi:pi/8:2*pi,'Normalization','probability','FaceColor','blue');
hold on
histogram(nfgcphsdiff,-2*pi:pi/8:2*pi,'Normalization','probability','FaceColor','yellow');
hold off
title(['all laps thetacyc'])%,...
    %sprintf(['t-test p = ' num2str(p_t) '\nr-test p = ' num2str(p_r)]))

xticks([-pi,0,pi])
xticklabels({'-π','0','π'})
set(gca,'FontSize',12)
axis square
    xlabel PhaseDiff
    ylabel('Relative probability')
%%
for i = 1:5
    subplot(1,6,i+1)

    histogram(PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,end)==1),1),-2*pi:pi/8:2*pi,...
        'Normalization','probability','FaceColor','blue');
    hold on
    histogram(PhaseDiff2((PhaseDiff2(:,3)==i & PhaseDiff2(:,end)==0),1),-2*pi:pi/8:2*pi,...
        'Normalization','probability','FaceColor','yellow');
    hold off
    
    
    xlabel PhaseDiff
    ylabel('Relative probability')
    xticks([-pi,0,pi])
    xticklabels({'-π','0','π'})
    axis square
    title(['Lap-' num2str(i)])%,...
    %sprintf(['t-test p = ' num2str(p_t) '\nr-test p = ' num2str(p_r)]))
    set(gca,'FontSize',12)
    ylim([0,0.28])
    
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
