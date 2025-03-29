file_output = 'Rat150_0620-1_thetaseq_phase\';
folder_name = 'prerunning_thetaseq\';
mkdir([file_output folder_name])
%%
nseg = 1;%哪个部分？ 1=pre-running, 2=sample, 3=test, 4=post
% nl = 2;%哪一圈？

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

%% 截取相位
% ind_phase = find
lap_EEG_Fs = 2000;
n_fft = 2^nextpow2(lap_timepoint_EEG);
power2 = abs(fft(lap_theta',n_fft)/lap_timepoint_EEG);
power1 = power2(1:n_fft/2,:);
power1(2:end-1,:) = 2*power1(2:end-1,:);
[~,ind_theta] = max(sum(power1,1));
% plot(0:(lap_EEG_Fs/n_fft):(lap_EEG_Fs/2-lap_EEG_Fs/n_fft),power1);

%% 计算theta sequence的开始（onset）和结束（offset）

% Default parameters
maxJump_thr = 30;  % in position bins
timeWin = 6;     % in time bin
timeStep = 1;     % in time bin
Distance_thr = 0; % in position bins
jump_prop_thr = 0;

[onset, offset, para] = DetectSequenceEvents_pos(scores{nseg}{nl,3},...
                    maxJump_thr,timeWin,timeStep,Distance_thr,jump_prop_thr);
%% 画原始解码图和路径

figure
set(gcf,'Position',[2 42 958 952]);
subplot(3,2,1:2)
imagesc(lap_ts,lap_ang,lap_pxn)
hold on
plot(lap_ts,lap_pos,'--w','LineWidth',1.5) % 画运动轨迹
stem(lap_ts(onset),repmat([10],size(onset,1),1),'-','LineWidth',1.5,'Marker','none')
stem(lap_ts(offset),repmat([10],size(offset,1),1),'-','LineWidth',1.5,'Marker','none')
hold off
colormap jet
colorbar('Position',[0.92,0.7093,0.02,0.2157])
caxis([0,0.25])
set(gca,'FontSize',14);

% limit = [2943.5,2947.5];% 2/3
% limit = [3129,3137];% 2/4
% limit = [1805.4,1811.3];% 1/3
xlim(limit)
axis xy


%% 画theta节律 找到相位和时间
length_EEG = length(lap_EEG_ts);
ind_onset = [];
ind_offset = [];
for i = 1:size(onset,1)
[~,ind_onset(i)] = min(abs(lap_EEG_ts-lap_ts(onset(i))));
end

for i = 1:size(offset,1)
[~,ind_offset(i)] = min(abs(lap_EEG_ts-lap_ts(offset(i))));
end

subplot(3,2,3:4)
plot(lap_EEG_ts,lap_theta(4,:)','k')
axis tight
hold on
stem(lap_EEG_ts(ind_onset),repmat([350],size(onset,1),1),'-','LineWidth',1.5,'Marker','none')
stem(lap_EEG_ts(ind_offset),repmat([350],size(offset,1),1),'-','LineWidth',1.5,'Marker','none')
hold off
set(gca,'FontSize',14);
xlim(limit)


% find the phase between timestamps of 'limit' 
phase_onset = lap_thetaphs(ind_onset);
phase_onset_ts = lap_EEG_ts(ind_onset);
phase_onset_lim = phase_onset(phase_onset_ts>= limit(1) ...
   & phase_onset_ts<=limit(2));
subplot(3,2,5)
polarhistogram(phase_onset_lim,18)
title('onset')


phase_offset = lap_thetaphs(ind_offset);
phase_offset_ts = lap_EEG_ts(ind_offset);
phase_offset_lim = phase_offset(phase_offset_ts>= limit(1) ...
   & phase_offset_ts<=limit(2));
subplot(3,2,6)
polarhistogram(phase_offset_lim,18)
title('offset')
% saveas(gcf,[file_output folder_name 'lap' num2str(nl)],'png')