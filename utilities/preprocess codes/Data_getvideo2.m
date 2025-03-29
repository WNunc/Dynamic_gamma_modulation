% to load video data in the session of circular track
% edited by Guo M and Zheng C, 2021/5/16

% clear
function Data_getvideo2(pathT)
%% load video data
% track=dir('*tracking.mat');
% pathT=track.name;
load(pathT);
track_center = [x_center,y_center];
track_center_orig = [x_center,y_center]./[x_sign,y_sign];
[t,x,y] = loadPos_rescaled_cz_v3('VT1.nvt',scale_x,scale_y,vfs,track_center_orig);
posx = x_sign*x-track_center(1);
posy = y_sign*y-track_center(2);

%% plot posx-posy on 2-dimention
ffa = figure('Units','normalized','Position',[0 0.1 0.8 0.8]);
figure(ffa)
subplot(2,1,1)
hold on
plot([-100,100],[0,0],'r--')
plot([0,0],[-100,100],'r--')
plot(posx,posy,'b.')
hold off
xlim([-80,80])
ylim([-80,80])
axis square
x0 = -80:40:80;
set(gca,'FontSize',20);
set(gca, 'XTick',x0);
set(gca, 'YTick',x0);
xlabel('Position x (cm)')
ylabel('Position y (cm)')
title('2-dimention map')
set(gcf,'PaperPositionMode','auto')

%% get angle data from the position data
ind0=find(posx==0);
posx(ind0)=0.01;
posang=nan(size(posx));
ind=find(posx>0);
posang(ind)=mod(atan(-posy(ind)./posx(ind)),2*pi);
ind=find(posx<0);
posang(ind)=mod(atan(-posy(ind)./posx(ind))+pi,2*pi);
rotation=rot;
if strcmp(rotation,'Counterclockwise')
    posang=2*pi-posang;
end
% use the on-track angle, treat the start-zone is located at angle=0
posang_ontrack=mod(posang-Ang_StartZone_center,2*pi);
posang_ontrack_unwrap = unwrap(posang_ontrack)';

%% Plot angle against time
figure(ffa)
subplot(2,1,2)
plot(t, posang_ontrack,'b.')
xlim(t([1,end]))
ylim([0,2*pi])
y0 = 0:pi/2:2*pi;
set(gca,'FontSize',20);
set(gca, 'YTick',x0);
set(gca, 'YTickLabel',{'0','pi/2','pi','3pi/2','2pi'});
xlabel('Time (s)')
ylabel('Angle (rad)')
title('1-dimention map')
set(gcf,'PaperPositionMode','auto')

%% calculate the running speed
% 1. 2D velocity
timelimit=120;   % ======= set the highest speed, unit is cm/s ========
vel=speed2D(posx,posy,t); %velocity in cm/s
vel(vel>=timelimit) = 0.5*(vel(circshift((vel>=timelimit),-3)) + vel(circshift((vel>=timelimit),3)));

% 2. 1D angular velocity
diam_center = (diam_outer+diam_inner);
vellimit = (2/diam_center)*100; % rad, approximately equals to 100cm/s
velrun = (2/diam_center)*5; % rad, approximately equals to 5cm/s, assuming the rat is running above this threshold
vel_ang = (posang_ontrack_unwrap(3:end)-posang_ontrack_unwrap(1:end-2))./(t(3:end)-t(1:end-2))';
vel_ang(vel_ang>=vellimit) = 0.5*(vel_ang(circshift((vel_ang>=vellimit),-3)) + vel_ang(circshift((vel_ang>=vellimit),3)));
vel_ang = [vel_ang(1);vel_ang;vel_ang(end)];

% 3. 2D acceleration
accelerate = (vel(3:end)-vel(1:end-2))./(t(3:end)-t(1:end-2))';
accelerate = [accelerate(1);accelerate;accelerate(end)];
accelerate_ang = (vel_ang(3:end)-vel_ang(1:end-2))./(t(3:end)-t(1:end-2))';
accelerate_ang = [accelerate_ang(1);accelerate_ang;accelerate_ang(end)];
%% Save data
data_video = nan(length(t),6);
data_video(:,1)=t';% Timesample
data_video(:,2)=posx';% x on circilar track
data_video(:,3)=posy';% y on circular track
data_video(:,4)=posang_ontrack';
data_video(:,5)=vel;% velocity
data_video(:,6)=vel_ang;% angle velocity
data_video(:,7)=accelerate;% acceleration
data_video(:,8)=accelerate_ang;% angle acceleration

scale=[scale_x,scale_y];% scale of x and y
save('Data_video.mat','data_video','vfs','scale');
end