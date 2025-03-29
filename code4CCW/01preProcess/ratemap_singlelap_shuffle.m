function [Ratemap,spikes,mapAxis,SInfo] = ratemap_singlelap_shuffle(celllist,trackdata,shuffle_num,vel_threshold)

% edited on 03/14/2017
% ratemap_singlelap_v1: 
% For each place cell, calculate single-lap ratemap
% including all pre-running laps, sample and test trials, and post-test
% trials
% Almost understand on 11/08/2020

if nargin < 4
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
track_center_orig = [x_center,y_center]./[x_sign,y_sign]; % Is sign for direction? orignal center?
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
vel(vel>=timelimit) = 0.5*(vel(circshift((vel>=timelimit),-3)) + vel(circshift((vel>=timelimit),3))); % why shift 3

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
path = pwd;% current path
NaFile = cell(1,numcells); % Store non existing file names here
NaFile = NaFile(1,1:0); % what for
S = loadSpikes(TT0,path,NaFile);
spikes=cell(numcells,1);
for nc=1:numcells
    if ~isa(S{nc},'ts') % Empty cell in this session
        ts = 1e64; % use a single ridicilous time stamp if the cell is silent
    else
        % Convert t-file data to timestamps in second?
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

%% shuffled Single lap rate map
% Ts_start_stop{1,1} = Ts_prerunning{1,1}./ 1000000;%CW
Ts_start_stop{1,1} = Ts_prerunning{1,2}./ 1000000;%CCW
% Ts_start_stop{2,1} = Ts_sample{1,1}./ 1000000;
% Ts_start_stop{3,1} = Ts_test{1,1}./ 1000000;
% Ts_start_stop{4,1} = Ts_posttest{1,1}./ 1000000;
nsegment = 1;% I only need the prerunning session
bin_ang = 2*pi/numbins;
mapAxis = bin_ang/2:bin_ang:(2*pi-bin_ang/2);
Ratemap = cell(1,nsegment);
SInfo = cell(1,nsegment);
for nseg = 1:nsegment
    Ts_start_stop0 = Ts_start_stop{nseg};
    nlap = size(Ts_start_stop0,1);
    for nl = 1:nlap
        ind = find(t_vel>=Ts_start_stop0(nl,1) & t_vel<=Ts_start_stop0(nl,2));
        t_nl = t_vel(ind);
        posang_ontrack_nl = posang_ontrack_vel(ind);
        Ratemap{nl,nseg} = [];
        dura = Ts_start_stop0(nl,2) - Ts_start_stop0(nl,1);
        mTS1 = dura*0.25;
        % mTS1 = 15;
        mTS2 = mTS1*2;
        for nc = 1:numcells
            ts = spikes{nc,1};
            ind = find(ts>=Ts_start_stop0(nl,1) & ts<=Ts_start_stop0(nl,2));
            time = mTS1:(dura-mTS2)/(shuffle_num+10):dura-mTS1;
            for nshf = 1:shuffle_num
                if length(ind) <= 1
                    % no spike firing in this lap
                    map_nl = zeros(numbins,1);
                    SI1 = 0;
                    SI2 = 0;
                else
                    ts_nl = ts(ind);
                    T_randind = randperm(length(time),1);
                    T_rand = time(T_randind);
                    disp(['random time: ' num2str(T_rand) 's'])
                    ts_nl_shuffle = ts_nl + T_rand;
                    ind_over = ts_nl_shuffle > Ts_start_stop0(nl,2);
                    temp_ts = ts_nl_shuffle(ind_over);
                    ts_nl_shuffle(ind_over) = temp_ts - Ts_start_stop0(nl,2) + Ts_start_stop0(nl,1);
                    ts_nl_shuffle = sort(ts_nl_shuffle);
                    
                    [tempspkang_ontrack_nl,~,Ind] = GetSpikePos(ts_nl_shuffle,posang_ontrack,posang,t);
                    %                 tempspkang_ontrack_nl = spikes{nc,2}(ind);
                    %                 map_nl = ratemap_ang_circ_cz(tempspkang_ontrack_nl,posang_ontrack_nl,t_nl,h,mapAxis,vfs);
                    map_nl = ratemap_ang_cz(tempspkang_ontrack_nl,posang_ontrack_nl,t_nl,h,mapAxis,vfs);
                    [SI1,SI2] = SpatialInfo(map_nl,tempspkang_ontrack_nl,mapAxis);
                end
                Ratemap{nl,nseg}{nshf}(:,nc) = map_nl;
                SInfo{nl,nseg}{nshf}(:,nc) = [SI1;SI2];
            end
            
            if length(ind) <= 1
                % no spike firing in this lap
                map_nl = zeros(numbins,1);
                SI1 = 0;
                SI2 = 0;
            else
                ts_nl = ts(ind);
                tempspkang_ontrack_nl = spikes{nc,2}(ind);
                % map_nl = ratemap_ang_circ_cz(tempspkang_ontrack_nl,posang_ontrack_nl,t_nl,h,mapAxis,vfs);
                map_nl = ratemap_ang_cz(tempspkang_ontrack_nl,posang_ontrack_nl,t_nl,h,mapAxis,vfs);
                [SI1,SI2] = SpatialInfo(map_nl,tempspkang_ontrack_nl,mapAxis);
            end
            SInfo{nl,nseg}{nshf+1}(:,nc) = [SI1;SI2];
        end
    end
end

