% 调试到thetawindows_all_CT_v4_EEG_wn.m 到16行
% 画相位直方图

figure
set(gcf,'Position',[587 228 678 366])
subplot(121)
hist(spkphases,xph)
xlabel('phase of theta')
ylabel('number of spikes')
title('spikes in session')
subplot(122)
hist(spkotphases,xph)
xlabel('phase of theta')
ylabel('number of spikes')
title('spikes on track')
saveas(gcf,[path '\Tseq\hist_phase_spkmin.png'])