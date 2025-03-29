function plot_decodeingandtheta_CW(scores,thetacycle,nl,TS_running,mth)
lap_ts = scores{nl,6};%行为学时间
lap_ang = scores{nl,5};%位置（角度）
lap_pxn = scores{nl,3};%解码位置（Pxn）
% lap_pos = scores{nl,5}(scores{nl,4});%实际位置
lap_pos = scores{nl,9};%实际位置
lap_EEG_ts = scores{nl,10};%EEG时间
lap_theta = scores{nl,18};%scores中的theta
IND1 = thetacycle.ind{nl};
limit = TS_running(nl,:);
set(gcf,'Position',[207 353 1438 416]);
set(axes,'Position', [0.05,0.45,0.88,0.42],'FontSize',14)
imagesc(lap_ts,lap_ang,lap_pxn)
hold on
plot(lap_ts,lap_pos,'--w','LineWidth',1.5) % 画运动轨迹
stem(lap_ts(IND1),repmat( [lap_ang(90)],1,length(IND1)),'w','LineWidth',1,'Marker','none')
hold off
colormap jet
colorbar('Position',[0.94,0.45,0.015,0.42])
caxis([0,0.25])
axis xy
xlim(limit)
title('decoding')
set(gca,'FontSize',11)
set(axes,'Position', [0.05,0.10,0.88,0.35],'FontSize',14)
plot(lap_EEG_ts,lap_theta(mth,:)','k')
xlim(limit)
set(gca,'FontSize',11)
maxy = get(gca,'YLim');
set(gca,'ytick',[]);
hold on
stem(lap_ts(IND1),repmat( [maxy(2)],1,length(IND1)),'r','LineWidth',1,'Marker','none','ShowBaseLine','off')
stem(lap_ts(IND1),repmat( [maxy(1)],1,length(IND1)),'r','LineWidth',1,'Marker','none','ShowBaseLine','off')
hold off


