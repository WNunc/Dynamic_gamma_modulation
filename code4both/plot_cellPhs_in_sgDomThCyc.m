% doall_getInfo_sgamma_v3_1 画图子程序
% 单个thetacycle中的单个cell画图

FFA = figure('Visible','off');
set(FFA,'Position',[771 180 395 629])
subplot(5,1,1:2)
imagesc(thseq_t,thseq_s,thseq);axis xy
ylabel('position(rad)')
colormap(jet)
xticklabels([])
xlim([thcyc_onset,thcyc_offset])
set(gca,'FontSize',9);
title(figname,'FontSize',14,'color',tcolor{c})
subplot(5,1,3:4)
scatter(t,p,sz, 'k|','LineWidth',1.5);%SPIKE
hold on
stem(t_sgcyc,720*ones(1,length(t_sgcyc)),...
    'r','LineWidth',1,'Marker','none','ShowBaseLine','off')
hold off
xlim([thcyc_onset,thcyc_offset])
ylim([0 720])
yticks([0,360,720])
yticklabels([0,360,720])
xticklabels([])
ylabel('gamma_s phase(deg)')
set(gca,'FontSize',9);
title(num2str(sgcycind'))
subplot(5,1,5)
plot(SGTT,SGwave)
maxy = get(gca,'YLim');
% set(gca,'ytick',[])
hold on
stem(t_sgcyc,maxy(2)*ones(1,length(t_sgcyc)),...
    'r','LineWidth',1,'Marker','none','ShowBaseLine','off')
stem(t_sgcyc,maxy(1)*ones(1,length(t_sgcyc)),...
    'r','LineWidth',1,'Marker','none','ShowBaseLine','off')
hold off
xlim([thcyc_onset,thcyc_offset])
xlabel('time(s)')
axis tight
ylabel('amp(μV)')
set(gca,'FontSize',9);

% pause
saveas(FFA,[outputFolder path_ns(13:end) figname '.png'])
saveas(FFA,[outputFolder path_ns(13:end) figname],'epsc')
clear FFA

