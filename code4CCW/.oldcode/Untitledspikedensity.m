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
%% Load spike firing data

celllist = 'H:\colgin data\Rat150\2017-06-12-CT-1\TTList_dCA1_pyr.txt';
fid=fopen(celllist);
if (fid == -1)
    warning([ 'Could not open tfile ' celllist]);
else
    % read the file names from the t-file list
    TT0 = ReadFileList(celllist);
    numCells0 = length(TT0);
    if numCells0==1 && max(TT0{1}==-1)
        % no cells in ttlist
        numcells=0;
    else
        numcells=numCells0;
    end
end
path = 'H:\colgin data\Rat150\2017-06-12-CT-1\';

NaFile = cell(1,numcells); % Store non existing file names here
NaFile = NaFile(1,1:0);
S = loadSpikes(TT0,path,NaFile);
spikes=cell(numcells,1);

for nc=1:numcells
    if ~isa(S{nc},'ts') % Empty cell in this session
        ts = 1e64; % use a single ridicilous time stamp if the cell is silent
    else
        % Convert t-file data to timestamps in second
        ts = Data(S{nc}) / 10000;
    end
    
%     [tempspkang_ontrack,tempspkang,Ind] = GetSpikePos(ts,posang_ontrack,posang,t);
%     spk_vel = vel(Ind);
%     spk_vel_ang = vel_ang(Ind);
%     
%     % use all spikes occurred when the running speed > vel_threshold
%     ind = find(spk_vel > vel_threshold);  % limit the running speed
    spikes{nc,1} = ts;
%     spikes{nc,1} = ts(ind);
%     spikes{nc,2} = tempspkang_ontrack(ind);
%     spikes{nc,3} = spk_vel(ind);
%     spikes{nc,4} = spk_vel_ang(ind);
end


%% 拼接所有spike
spikesall = [];
for i = 1:numcells
    spikesall = [spikesall;spikes{i}];
end

spikesall_sort = sort(spikesall);
indspike = find(spikesall_sort>=limit(1) & spikesall_sort<=limit(2));

spikesinlap = spikesall_sort(indspike);
step = 0.005 ;%s
win = 0.020 ;%s
i=1;
numspikes = [];
densityts = [];
while limit(1)+(i-1)*step+win<=limit(2)

    numspikes(i)=length(find(spikesinlap>=limit(1)+(i-1)*step & spikesinlap<=limit(1)+(i-1)*step+win));
    densityts(i) = limit(1)+(i-1)*step+0.5*win;
    i=i+1;
end



densityspikes = numspikes./win;
densityspk = smooth(densityspikes,10);
[~,ind_pk] = findpeaks(densityspk,'MinPeakProminence',50,'MinPeakDistance',5);
% [~,ind_pk] = findpeaks(densityspk,'MinPeakProminence',50);
% 找到峰值ind_pk对应的phase
ind_pkts = [];
for i = 1:length(ind_pk)
[~,ind_pkts(i)] = min(abs(lap_EEG_ts-densityts(ind_pk(i))));
end
% 找到峰值ind_pk对应的pos
ind_pkpos = [];
for i = 1:length(ind_pk)
[~,ind_pkpos(i)] = min(abs(lap_ts-densityts(ind_pk(i))));
end
%% 画原始解码图和路径

figure
set(gcf,'Position',[2 42 1438 952]);
subplot(3,1,1)
imagesc(lap_ts,lap_ang,lap_pxn)
hold on
plot(lap_ts,lap_pos,'--w','LineWidth',1.5) % 画运动轨迹
stem(lap_ts(ind_pkpos),repmat( [lap_ang(90)],1,length(ind_pkpos)),'w','LineWidth',1,'Marker','none')
hold off
colormap jet
colorbar('Position',[0.91 0.7093,0.015,0.21])
caxis([0,0.25])
xlim(limit)
axis xy
title('decoding')
% yyaxis right
% hold on
% plot(densityts,densityspk,'w')
% scatter(densityts(ind_pk),densityspk(ind_pk))
% hold off
subplot(3,1,2)
plot(lap_EEG_ts,lap_theta(4,:)','k')
hold on
scatter(lap_EEG_ts(ind_pkts),lap_theta(4,ind_pkts)')
% stem(lap_EEG_ts(ind_pkts),repmat([350],size(onset,1),1),'-','LineWidth',1.5,'Marker','none')
hold off
xlim(limit)
title('theta')
subplot(3,1,3)
plot(densityts,densityspk,'k')
hold on
scatter(densityts(ind_pk),densityspk(ind_pk))
hold off
phase_spkdensity = lap_thetaphs(ind_pkts);
xlim(limit)
title('SPK density')
% polarhistogram(phase_spkdensity,18)
% title('SPK density peak phase')