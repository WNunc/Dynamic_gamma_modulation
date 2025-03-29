% figure;
% plot(1:3827,abs(error));
% hold on
% scatter(idxpE,abs(error(idxpE)))
% hold off
% xlim([3400,3827])
% 画老师的scores
% 解码误差最大处分割theta序列 

file_output = 'Rat148_0808-2_error_in_theta_phase\';
folder_name = 'prerunning_PeakProminence0.3\';
mkdir([file_output folder_name])
%%
nseg = 1;%哪个部分？ 1=pre-running, 2=sample, 3=test, 4=post
nl = 2;%哪一圈？

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
%% 画原始解码图和路径

figure
yyaxis left
imagesc(lap_ts,lap_ang,lap_pxn)
hold on
plot(lap_ts,lap_pos,'--w','LineWidth',2) % 画运动轨迹
hold off
colormap jet
colorbar
caxis([0,0.25])

% limit = [2943.5,2947.5];% 2/3
% limit = [3129,3137];% 2/4
% limit = [1805.4,1811.3];% 1/3
xlim(limit)
axis xy
yyaxis right
plot(densityts,densityspk,'w')

%% 画解码概率质心的连线
lap_angs = repmat(lap_ang,[1,lap_timepoint_track]);
lap_pxnCOM = sum(lap_pxn .* lap_angs,1) ./ sum(lap_pxn,1);
error_COM0 = abs(lap_pxnCOM' - lap_pos);% 含有nan数据
ind_notnan = find(~isnan(error_COM0));
error_COM = error_COM0(ind_notnan);% 不包含nan数据
[~,ind_pk] = findpeaks(error_COM,'MinPeakProminence',th,'MinPeakDistance',5);
% [~,ind_pk] = findpeaks(error_COM,'MinPeakProminence',th);
ind_COM = ind_notnan(ind_pk);% 解码误差最大的那些点
%======找到相位======
L_ind_COM = length(ind_COM);
ind_COMerror_phs=zeros(1,L_ind_COM);
for s=1:L_ind_COM
    dis_t=abs(lap_EEG_ts-lap_ts(ind_COM(s)));
    [~,ind_COMerror_phs(1,s)]=min(dis_t);
end
phase_COMerror = lap_thetaphs(ind_COMerror_phs);
phase_COMerror_ts = lap_EEG_ts(ind_COMerror_phs);
phase_COMerror_OK = phase_COMerror(phase_COMerror_ts>=limit(1) ...
    & phase_COMerror_ts<=limit(2));

figure
set(gcf,'Position',[2 42 958 952]);
set(gca,'FontSize',14);
suptitle('COM')
subplot(3,1,1)
imagesc(lap_ts,lap_ang,lap_pxn)
hold on
plot(lap_ts,lap_pxnCOM,'w','LineWidth',2) % 画COM轨迹
plot(lap_ts,lap_pos,'--g','LineWidth',2) % 画运动轨迹
scatter(lap_ts(ind_COM),lap_pxnCOM(ind_COM))
hold off
colormap jet
% colorbar
caxis([0,0.25])
xlim(limit)
axis xy

subplot(3,1,2)
plot(lap_ts,error_COM0);
hold on
scatter(lap_ts(ind_COM),error_COM0(ind_COM))
hold off
xlim(limit)
subplot(3,1,3)
plot(lap_EEG_ts,lap_theta(4,:)','k')
hold on
scatter(lap_ts(ind_COM),lap_theta(4,ind_COMerror_phs)')
hold off
xlim(limit)
saveas(gcf,[file_output folder_name 'lap' num2str(nl)],'png')
%% 画解码概率最大值的连线
[~,ind_max] = max(lap_pxn);
lap_pxnPP = [];
for i = 1:lap_timepoint_track
    lap_pxnPP(i) = lap_pxn(ind_max(i),i);
end

lap_pxnPEAK = lap_ang(ind_max);
lap_pxnPEAK(isnan(lap_pxnPP)) = nan;

error_PEAK0 = abs(lap_pxnPEAK- lap_pos);% 含有nan数据
ind_notnan = find(~isnan(error_PEAK0));
error_PEAK = error_PEAK0(ind_notnan);% 不包含nan数据
[~,ind_pk] = findpeaks(error_PEAK,'MinPeakProminence',th,'MinPeakDistance',5);
% [~,ind_pk] = findpeaks(error_PEAK,'MinPeakProminence',th);
ind_PEAK = ind_notnan(ind_pk);% 解码误差最大的那些点

%======找到相位======
L_ind_PEAK = length(ind_PEAK);
ind_PEAKerror_phs=zeros(1,L_ind_PEAK);
for s=1:L_ind_PEAK
    dis_t=abs(lap_EEG_ts-lap_ts(ind_PEAK(s)));
    [~,ind_PEAKerror_phs(1,s)]=min(dis_t);
end
phase_PEAKerror=lap_thetaphs(ind_PEAKerror_phs);
phase_PEAKerror_ts = lap_EEG_ts(ind_PEAKerror_phs);
phase_PEAKerror_OK = phase_PEAKerror(phase_PEAKerror_ts>=limit(1) ...
    & phase_PEAKerror_ts<=limit(2));



figure
set(gcf,'Position',[962 42 958 952]);
suptitle('PPP')
subplot(3,1,1)
imagesc(lap_ts,lap_ang,lap_pxn)
hold on
plot(lap_ts,lap_pxnPEAK,'w','LineWidth',2) % 画PEAK轨迹
plot(lap_ts,lap_pos,'--g','LineWidth',2) % 画运动轨迹
scatter(lap_ts(ind_PEAK),lap_pxnPEAK(ind_PEAK))
hold off
colormap jet
% colorbar
caxis([0,0.25])
xlim(limit)
axis xy

subplot(3,1,2)
plot(lap_ts,error_PEAK0);
hold on
scatter(lap_ts(ind_PEAK),error_PEAK0(ind_PEAK))
hold off
xlim(limit)
subplot(3,1,3)
plot(lap_EEG_ts,lap_theta(4,:)','k')
hold on
scatter(lap_ts(ind_PEAK),lap_theta(4,ind_PEAKerror_phs)')
hold off
xlim(limit)
