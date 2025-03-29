%% 1. plotdata

phaseCOM_total = [];
phasePEAK_total = [];
th = 0.3;
close all
% load(trackdata)
% load(Data_angle_ontrack.mat)
for nl = 1:6
plotdata
phaseCOM_total = [phaseCOM_total,phase_COMerror_OK];
phasePEAK_total = [phasePEAK_total,phase_PEAKerror_OK];
end
nbars = 36;
figure;histogram(phaseCOM_total,nbars);title('COM')
figure
[A,B] = hist(phaseCOM_total,nbars); %#ok<*HIST>
bar(B,A,'histc')
hold on
plot(B,smooth(A,round(nbars/4)),'r','LineWidth',1.5)
hold off
title('COM')

% saveas(gcf,[file_output folder_name 'COM'],'png')
figure;histogram(phasePEAK_total,nbars);title('PPP')
figure
[A,B] = hist(phasePEAK_total,nbars); %#ok<*HIST>
bar(B,A,'histc')
hold on
plot(B,smooth(A,round(nbars/4)),'r','LineWidth',1.5)
hold off
title('PPP')

% saveas(gcf,[file_output folder_name 'PPP'],'png')


%% 2. plot_thetaseq

% load(trackdata)
% load(Data_angle_ontrack.mat)
close all
phase_onset_total = [];
phase_offset_total = [];
for nl = 1:6
plot_thetaseq
phase_onset_total = [phase_onset_total,phase_onset_lim];
phase_offset_total = [phase_offset_total,phase_offset_lim];
end
figure
set(gcf,'Position',[490 310 1010 420])
subplot(1,2,1)
polarhistogram(phase_onset_total,36);title('onset')
subplot(1,2,2)
polarhistogram(phase_offset_total,36);title('offset')


%% 3.1 spike density
close all
phase_total = [];
for nl = 1:6
Untitledspikedensity;
phase_total = [phase_total phase_spkdensity];
saveas(gcf,['lap' num2str(nl) '-cut.png'])
saveas(gcf,['lap' num2str(nl) '-cut.fig'])
end
figure
nbars = 36;
[A,B] = hist(phase_total,nbars); %#ok<*HIST>  %A = count; B = center(phase)
histogram(phase_total,nbars);title('phase')
hold on
plot(B,smooth(A,round(nbars/4)),'r','LineWidth',1.5)
hold off
saveas(gcf,['all-lap-phase-cut.png'])

[~,indp]= max(smooth(A,round(nbars/4)));
max_phase = B(indp);
%% 3.2 after spike density
nseg = 1;%哪个部分？ 1=pre-running, 2=sample, 3=test, 4=post
% nl = 2;%哪一圈？
for nl = 1:6
lap_ts = scores{nseg}{nl,6};%行为学时间
lap_ang = scores{nseg}{nl,5};%位置（角度）
lap_pxn = scores{nseg}{nl,3};%解码位置（Pxn）
lap_pos = scores{nseg}{nl,5}(scores{nseg}{nl,4});%实际位置
lap_EEG = scores{nseg}{nl,11};%EEG数据
lap_EEG_ts = scores{nseg}{nl,10};%EEG时间
lap_theta = scores{nseg}{nl,18};%scores中的theta
lap_thetaphs = scores{nseg}{nl,21};%scores中的theta相位
lap_timepoint_track = size(lap_pxn,2);%路径总时间点
lap_timepoint_EEG = size(lap_EEG,2);%EEG总时间点
lap_speed = scores{nseg}{nl,8}(1,:);
%% 载入行为学数据,根据行为学筛选运动时的时间
% load Data_angle_ontrack.mat
speed = data_angle{nseg}{nl}(:,3);% cm/s
speed_ang = data_angle{nseg}{nl}(:,4);% rad/s

ind_ts_start = max(find(lap_speed'>=10 & lap_pos<=Ang_RewardLoc_ontrack(1)));
ind_ts_end = max(find(lap_speed'>=10 & lap_pos<=Ang_RewardLoc_ontrack(19)));
limit = [lap_ts(ind_ts_start), lap_ts(ind_ts_end)];
%% 用最大spike密度处的相位切割theta
% 对应的时间点

% [~,I] = min(abs(lap_thetaphs(4,:)-max_phase),[],1);
[~,I]=findpeaks(-abs(lap_thetaphs(4,:)-max_phase),'MinPeakHeight',-5);

% 对应的位置
IND = [];
for i= 1:length(I)
    [~,IND(i)] = min(abs(lap_ts - lap_EEG_ts(I(i))));
end
figure
set(gcf,'Position',[2 42 1438 952]);
subplot(2,1,1)
imagesc(lap_ts,lap_ang,lap_pxn)
hold on
plot(lap_ts,lap_pos,'--w','LineWidth',1.5) % 画运动轨迹
stem(lap_ts(IND),repmat( [lap_ang(90)],1,length(IND)),'w','LineWidth',1,'Marker','none')
hold off
colormap jet
colorbar('Position',[0.91 0.7093,0.015,0.21])
caxis([0,0.25])
axis xy
xlim(limit)

title('decoding')
subplot(2,1,2)
plot(lap_EEG_ts,lap_theta(4,:)','k')
hold on
scatter(lap_EEG_ts(I),lap_theta(4,I)')
hold off
xlim(limit)

end
