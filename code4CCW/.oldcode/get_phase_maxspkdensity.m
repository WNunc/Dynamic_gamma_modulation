%% calulate the MAX density of all spike in TTlist
% created by WN on 2022/01/22
% output:
% 1.
% 2.



%% 单圈运行启用
nseg = 1;%哪个部分？ 1=pre-running, 2=sample, 3=test, 4=post
% nl = 2;%哪一圈？

% % Load spike firing data
% celllist = [path 'TTList_dCA1_pyr.txt'];
% fid=fopen(celllist);
% if (fid == -1)
%     warning([ 'Could not open tfile ' celllist]);
% else
%     % read the file names from the t-file list
%     TT0 = ReadFileList(celllist);
%     numCells0 = length(TT0);
%     if numCells0==1 && max(TT0{1}==-1)
%         % no cells in ttlist
%         numcells=0;
%     else
%         numcells=numCells0;
%     end
% end
% NaFile = cell(1,numcells); % Store non existing file names here
% NaFile = NaFile(1,1:0);
% S = loadSpikes(TT0,path,NaFile);
% spikes=cell(numcells,1);

%%
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
%
[TSlimit, speed, speed_ang]= TSfind_laplimit(data_angle,Ang_RewardLoc_ontrack,lap_speed,lap_pos,lap_ts,nseg,nl);

% %%
% for nc=1:numcells
%     if ~isa(S{nc},'ts') % Empty cell in this session
%         ts = 1e64; % use a single ridicilous time stamp if the cell is silent
%     else
%         % Convert t-file data to timestamps in second
%         ts = Data(S{nc}) / 10000;
%     end
%     
%     %     [tempspkang_ontrack,tempspkang,Ind] = GetSpikePos(ts,posang_ontrack,posang,t);
%     %     spk_vel = vel(Ind);
%     %     spk_vel_ang = vel_ang(Ind);
%     %
%     %     % use all spikes occurred when the running speed > vel_threshold
%     %     ind = find(spk_vel > vel_threshold);  % limit the running speed
%     spikes{nc,1} = ts;
%     %     spikes{nc,1} = ts(ind);
%     %     spikes{nc,2} = tempspkang_ontrack(ind);
%     %     spikes{nc,3} = spk_vel(ind);
%     %     spikes{nc,4} = spk_vel_ang(ind);
% end
    % find cell on track
    spikes = S.spikes;
%     peak0 = 1;
%     peak_all = max(S.Ratemap_singlelap{nl,1}); %1 是segment的编号
%             ind = find(peak_all >= peak0);
%             spikes = spikes(ind,:);
%     X = cell2struct(S.fieldProp_singlelap{nl,1}(ind),'placefield',1);%1 是segment的编号
% %     load(trackdata_ns,'Ang_RewardLoc')
%     COM = [];
%     for ncell  = 1:length(ind)
%         COM(ncell) = X(ncell).placefield(1).x_COM;
%     end
%     
%     ind_ot = find(COM>Ang_RewardLoc_ontrack(1) & COM<Ang_RewardLoc_ontrack(18));
%     spikes = spikes(ind_ot,1);
    numcells = size(spikes,1);
%% 拼接所有spike
spikesall = [];
for i = 1:numcells
    spikesall = [spikesall;spikes{i}];
end
spikesall_sort = sort(spikesall);
indspike = find(spikesall_sort>=TSlimit(1) & spikesall_sort<=TSlimit(2));
%
spikesinlap = spikesall_sort(indspike);
step = 0.005 ;%s
win = 0.020 ;%s
i=1;
numspikes = [];
densityts = [];
while TSlimit(1)+(i-1)*step+win<=TSlimit(2)
    numspikes(i)=length(find(spikesinlap>=TSlimit(1)+(i-1)*step & spikesinlap<=TSlimit(1)+(i-1)*step+win));
    densityts(i) = TSlimit(1)+(i-1)*step+0.5*win;
    i=i+1;
end

densityspikes = numspikes./win;
densityspk = smooth(densityspikes,10);
[~,ind_pk] = findpeaks(densityspk,'MinPeakProminence',40,'MinPeakDistance',5);
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

% select max theta FFT channel
% thetafft = abs(fft(lap_theta));
% thetafftsum = sum(thetafft,2);
% [~,mth] = max(thetafftsum);

disp(['theta rythm on CSC' num2str(CSClist_CA1{ns}(mth)) ' accept'])

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
xlim(TSlimit)
axis xy
title('decoding')
% yyaxis right
% hold on
% plot(densityts,densityspk,'w')
% scatter(densityts(ind_pk),densityspk(ind_pk))
% hold off
subplot(3,1,2)
plot(lap_EEG_ts,lap_theta(mth,:)','k')
hold on
scatter(lap_EEG_ts(ind_pkts),lap_theta(mth,ind_pkts)')
% stem(lap_EEG_ts(ind_pkts),repmat([350],size(onset,1),1),'-','LineWidth',1.5,'Marker','none')
hold off
xlim(TSlimit)
title('theta')
subplot(3,1,3)
plot(densityts,densityspk,'k')
hold on
scatter(densityts(ind_pk),densityspk(ind_pk))
hold off
phase_spkdensity = lap_thetaphs(mth,ind_pkts);
phase_mxspk{nl} = phase_spkdensity;
xlim(TSlimit)
title('SPK density')
% polarhistogram(phase_spkdensity,18)
% title('SPK density peak phase')







%% function
function [TSlimit, speed, speed_ang]= TSfind_laplimit(data_angle,Ang_RewardLoc_ontrack,lap_speed,lap_pos,lap_ts, nsegment,nlap)
% input:
% data_angle ---- from file named Data_angle_ontrack.mat
% Ang_RewardLoc_ontrack ---- from file named like date_CT_tracking.mat
% output:
% TSlimit ---- the start and the end timgstamp of one lap
% speed ---- raw speed from data_angle
% speed_ang ---- raw angle speed from data_angle
speed = data_angle{nsegment}{nlap}(:,3);% cm/s
speed_ang = data_angle{nsegment}{nlap}(:,4);% rad/s
ind_ts_start = max(find(lap_speed'>=10 & lap_pos<=Ang_RewardLoc_ontrack(1)));
ind_ts_end = max(find(lap_speed'>=10 & lap_pos<=Ang_RewardLoc_ontrack(18)));
TSlimit = [lap_ts(ind_ts_start), lap_ts(ind_ts_end)];
end