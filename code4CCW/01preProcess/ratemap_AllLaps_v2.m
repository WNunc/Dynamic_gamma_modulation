function [Ratemap,spikes,mapAxis] = ratemap_AllLaps_v2(celllist,trackdata,vel_threshold)

% edited on 03/16/2017
% ratemap_AllLaps_v1: 
% For each place cell, calculate overall rate map of all running laps
% including all pre-running laps, sample and test trials, and post-test
% trials

if nargin < 3
    vel_threshold = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read tracking data
load(trackdata);

numbins = 90; % Number of bins
h = 0.1; % Smoothing factor when calculating the ratemap

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load position information
track_center = [x_center,y_center];
track_center_orig = [x_center,y_center]./[x_sign,y_sign];
[t,x,y] = loadPos_rescaled_cz_v3('VT1.nvt',scale_x,scale_y,vfs,track_center_orig);
posx = x_sign*x-track_center(1);
posy = y_sign*y-track_center(2);

% get angle data from the position data
ind0=find(posx==0);
posx(ind0)=0.01;
posang=nan(size(posx));
ind=find(posx>0);
posang(ind)=mod(atan(-posy(ind)./posx(ind)),2*pi);
ind=find(posx<0);
posang(ind)=mod(atan(-posy(ind)./posx(ind))+pi,2*pi);
rotation=rot;
if strcmp(rotation,'Counterclockwise');
    posang=2*pi-posang;
end
% use the on-track angle, treat the start-zone is located at angle=0
posang_ontrack=mod(posang-Ang_StartZone_center,2*pi);
posang_ontrack_unwrap = unwrap(posang_ontrack)';

%% 2-D running speed, acceleration, jerk, and angular velocity
% calculate the running speed
timelimit=120;   % ======= set the highest speed, unit is cm/s ========
vel=speed2D(posx,posy,t); %velocity in cm/s
vel(vel>=timelimit) = 0.5*(vel(circshift((vel>=timelimit),-3)) + vel(circshift((vel>=timelimit),3)));

% calculate the angular velocity
vellimit = 2; % rad, approximately equals to 100cm/s
velrun = 0.1; % rad, approximately equals to 5cm/s, assuming the rat is running above this threshold
vel_ang = (posang_ontrack_unwrap(3:end)-posang_ontrack_unwrap(1:end-2))./(t(3:end)-t(1:end-2))';
vel_ang(vel_ang>=vellimit) = 0.5*(vel_ang(circshift((vel_ang>=vellimit),-3)) + vel_ang(circshift((vel_ang>=vellimit),3)));
vel_ang = [vel_ang(1);vel_ang;vel_ang(end)];

% calculate the acceleration
accelerate = (vel(3:end)-vel(1:end-2))./(t(3:end)-t(1:end-2))';
accelerate = [accelerate(1);accelerate;accelerate(end)];
accelerate_ang = (vel_ang(3:end)-vel_ang(1:end-2))./(t(3:end)-t(1:end-2))';
accelerate_ang = [accelerate_ang(1);accelerate_ang;accelerate_ang(end)];

% calculate the jerk
jerk = (accelerate(3:end)-accelerate(1:end-2))./(t(3:end)-t(1:end-2))';
jerk = [jerk(1);jerk;jerk(end)];
jerk_ang = (accelerate_ang(3:end)-accelerate_ang(1:end-2))./(t(3:end)-t(1:end-2))';
jerk_ang = [jerk_ang(1);jerk_ang;jerk_ang(end)];

% limit the posang_ontrack with running speed > vel_threshold
ind = find(vel > vel_threshold);
t_vel = t(ind);
posang_ontrack_vel = posang_ontrack(ind);

%% Load spike firing data
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
path = pwd;
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
    
    [tempspkang_ontrack,tempspkang,Ind] = GetSpikePos(ts,posang_ontrack,posang,t);
    spk_vel = vel(Ind);
    spk_vel_ang = vel_ang(Ind);
    
    % use all spikes occurred when the running speed > vel_threshold
    ind = find(spk_vel > vel_threshold);  % limit the running speed
    
    spikes{nc,1} = ts(ind);
    spikes{nc,2} = tempspkang_ontrack(ind);
    spikes{nc,3} = spk_vel(ind);
    spikes{nc,4} = spk_vel_ang(ind);
end

%% for each cell, rate map across all running laps
% Ts_start_stop{1,1} = Ts_prerunning{1,1}./ 1000000;% CW
Ts_start_stop{1,1} = Ts_prerunning{1,2}./ 1000000;% CCW
% Ts_start_stop{1,2} = Ts_sample{1,1}(:,[1,3])./ 1000000;
% Ts_start_stop{1,3} = Ts_test{1,1}(:,[1,3])./ 1000000;
% Ts_start_stop{1,4} = Ts_posttest{1,1}(:,[1,3])./ 1000000;
nsegment = 1;
bin_ang = 2*pi/numbins;
mapAxis = bin_ang/2:bin_ang:(2*pi-bin_ang/2);
Ratemap = zeros(numbins,numcells);
for nc = 1:numcells
    ts = spikes{nc,1};
    
    t_all = [];
    posang_ontrack_all = [];
    ts_all = [];
    tempspkang_ontrack_all = [];
    for nseg = 1:nsegment
        Ts_start_stop0 = Ts_start_stop{nseg};
        nlap = size(Ts_start_stop0,1);
        for nl = 1:nlap
            ind = find(t_vel>=Ts_start_stop0(nl,1) & t_vel<=Ts_start_stop0(nl,2));
            t_nl = t_vel(ind);
            posang_ontrack_nl = posang_ontrack_vel(ind);
            t_all = [t_all,t_nl];
            posang_ontrack_all = [posang_ontrack_all,posang_ontrack_nl];
            
            ind = find(ts>=Ts_start_stop0(nl,1) & ts<=Ts_start_stop0(nl,2));
            ts_nl = ts(ind);
            tempspkang_ontrack_nl = spikes{nc,2}(ind);
            ts_all = [ts_all;ts_nl];
            tempspkang_ontrack_all = [tempspkang_ontrack_all;tempspkang_ontrack_nl];
        end
    end
    
    if length(ts_all) <= 1
        % no spike firing in this lap
        map_all = ones(numbins,1)*0.01;
    else
        map_all = ratemap_ang_cz(tempspkang_ontrack_all,posang_ontrack_all,t_all,h,mapAxis,vfs);
    end
    Ratemap(:,nc) = map_all;
end