function [scores,Ts_start_stop] = decodewindows_CT_v3(trackdata,spikes,Ratemap,mincells,minspikes,decodwin,trialid,lapid)

% Edited on 2021/1/20
% add trialid as input


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read video data
path = pwd;
load([path '/Data_video.mat']);
% Read tracking data
load(trackdata,'Ts_prerunning','Ts_sample','Ts_test','Ts_posttest','rot');

switch trialid
    case 1
        % Ts_start_stop{1,1} = Ts_prerunning{1,1}./ 1000000; % CW
        Ts_start_stop{1,1} = Ts_prerunning{1,2}./ 1000000; % CCW
    case 2
        Ts_start_stop{1,1} = Ts_sample{1,1}./ 1000000;
    case 3
        Ts_start_stop{1,1} = Ts_test{1,1}./ 1000000;
    case 4
        Ts_start_stop{1,1} = Ts_posttest{1,1}./ 1000000;
end

nsegment = trialid;

% binsize = 4; % in degree
numbins = 90; % Number of bins
bin_ang = 2*pi/numbins;
mapAxis = bin_ang/2:bin_ang:(2*pi-bin_ang/2);
% dt = .04;
% step = .01;
dt = decodwin(1);
step = decodwin(2);
h = 0.1; % Smoothing factor when calculating the ratemap

numcells=size(spikes,1);

% %% load position information
% track_center = [x_center,y_center];
% track_center_orig = [x_center,y_center]./[x_sign,y_sign];
% [t,x,y] = loadPos_rescaled_cz_v3('VT1.nvt',scale_x,scale_y,vfs,track_center_orig);
% posx = x_sign*x-track_center(1);
% posy = y_sign*y-track_center(2);
% 
% % get angle data from the position data
t = data_video(:,1);
posx = data_video(:,2);
posy = data_video(:,3);

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
% 
% % use the on-track angle, treat the start-zone is located at angle=0
% posang_ontrack=mod(posang-Ang_StartZone_center,2*pi);
% posang_ontrack_unwrap = unwrap(posang_ontrack)';
% 
% %% 2-D running speed, acceleration, jerk, and angular velocity
% % calculate the running speed
% timelimit=120;   % ======= set the highest speed, unit is cm/s ========
% vel=speed2D(posx,posy,t); %velocity in cm/s
% vel(vel>=timelimit) = 0.5*(vel(circshift((vel>=timelimit),-3)) + vel(circshift((vel>=timelimit),3)));
% 
% % calculate the angular velocity
% vellimit = 2; % rad, approximately equals to 100cm/s
% velrun = 0.1; % rad, approximately equals to 5cm/s, assuming the rat is running above this threshold
% vel_ang = (posang_ontrack_unwrap(3:end)-posang_ontrack_unwrap(1:end-2))./(t(3:end)-t(1:end-2))';
% vel_ang(vel_ang>=vellimit) = 0.5*(vel_ang(circshift((vel_ang>=vellimit),-3)) + vel_ang(circshift((vel_ang>=vellimit),3)));
% vel_ang = [vel_ang(1);vel_ang;vel_ang(end)];
% 
% % calculate the acceleration
% accelerate = (vel(3:end)-vel(1:end-2))./(t(3:end)-t(1:end-2))';
% accelerate = [accelerate(1);accelerate;accelerate(end)];
% accelerate_ang = (vel_ang(3:end)-vel_ang(1:end-2))./(t(3:end)-t(1:end-2))';
% accelerate_ang = [accelerate_ang(1);accelerate_ang;accelerate_ang(end)];
% 
% % calculate the jerk
% jerk = (accelerate(3:end)-accelerate(1:end-2))./(t(3:end)-t(1:end-2))';
% jerk = [jerk(1);jerk;jerk(end)];
% jerk_ang = (accelerate_ang(3:end)-accelerate_ang(1:end-2))./(t(3:end)-t(1:end-2))';
% jerk_ang = [jerk_ang(1);jerk_ang;jerk_ang(end)];
% 
% % limit the posang_ontrack with running speed > vel_threshold
% vel_threshold = 5; % cm/s  %IMPORTANT: change to 0 if do not want to limit
% ind = find(vel > vel_threshold);
% t_vel = t(ind);
% posang_ontrack_vel = posang_ontrack(ind);
vel = data_video(:,5);
vel_ang = data_video(:,6);
accelerate = data_video(:,7);
accelerate_ang = data_video(:,8);

posang_ontrack = data_video(:,4);
%% get p(x) (first iteration of inference)
% Pxx = probX(posang_ontrack_vel,bin_ang);
% method 1: to simply use the distibution ofposx (as in jensen,lisman)
%   Px = Pxx;
% method 2: to contol for directoin ,and use uniform prior
Px = ones(numbins,1)./numbins;

%% get p(x|n)
scores = {};
for nseg = nsegment
    Ts_start_stop0 = Ts_start_stop{1};
    nlap = size(Ts_start_stop0,1);
    if nlap<lapid
     continue;
    end % wn 
    scores{1,nseg} = {};
    for nl = lapid
        tstart = Ts_start_stop0(nl,1);
        tstop = Ts_start_stop0(nl,2);
        % decoding
        tstartstep = tstart;
        tstopstep = tstart + dt;
        nstep = 0;
        Pxn = nan(numbins,1);
        xbin = nan(0);
        xpos = nan(0);
        tbin = nan(0);
        vel_nl = nan(0);
        acc_nl = nan(0);
        jerk_nl = nan(0);
        vel_ang_nl = nan(0);
        acc_ang_nl = nan(0);
        jerk_ang_nl = nan(0);
        while(tstopstep <= tstop)
            nstep = nstep+1;
            % spikes for the small time window (dt time step)
            n = zeros(1,numcells);
            for i = 1:numcells
                n(1,i) = length(find(spikes{i,1}>tstartstep & spikes{i,1}<=tstopstep));
            end
            
            %get probabilty the animals was in each bin
            if sum(n) >= minspikes && length(find(n>0)) >= mincells
                for i=1:numbins
                    for j = 1:numcells
                        Pnx(j) = ((dt*Ratemap(i,j))^n(1,j))/factorial(n(1,j));
                        Pnx(j) = Pnx(j) * exp(-dt*Ratemap(i,j));
                    end
                    Pxn(i,nstep) =  Px(i) * prod(Pnx);
                    clear Pnx;
                end
            else
                Pxn(:,nstep) =  nan(numbins,1);
            end
            
            %get "measured position" in terms of bin for the time window.
            t0 = (tstartstep+tstopstep)/2;
            [ang0,~,~] = GetSpikePos(t0,posang_ontrack,posang,t);
            [~,ind] = min(abs(mapAxis - ang0));
            xbin(1,nstep) = ind; %the measured bin
            xpos(1,nstep) = ang0;
            tbin(1,nstep) = tstartstep;
            ind_dt = find(t >= tstartstep & t <= tstopstep);
            if isempty(ind_dt)
                [~,ind_dt] = min(abs(t-t0));
            end
            vel_nl(1,nstep) = mean(vel(ind_dt));
            acc_nl(1,nstep) = mean(accelerate(ind_dt));
%             jerk_nl(1,nstep) = mean(jerk(ind_dt));
            vel_ang_nl(1,nstep) = mean(vel_ang(ind_dt));
            acc_ang_nl(1,nstep) = mean(accelerate_ang(ind_dt));
%             jerk_ang_nl(1,nstep) = mean(jerk_ang(ind_dt));
            
            tstartstep = tstartstep + step;
            tstopstep = tstopstep + step;
        end
        
        % get normalized version of Pxn
        mx = nansum(Pxn);
        for nstep = 1:size(Pxn,2)
            Pxn(:,nstep) = Pxn(:,nstep)/mx(nstep);
        end
        
        % save in scores
        scores{nseg}{nl,1} = [nl]; % trial id, sample trial id
        scores{nseg}{nl,2} = [tstart,tstop];
        scores{nseg}{nl,3} = Pxn;
        scores{nseg}{nl,4} = xbin;% actual position bin
        scores{nseg}{nl,5} = mapAxis';% position of each bin
        scores{nseg}{nl,6} = tbin;% time point of each bin
        scores{nseg}{nl,8} = [vel_nl;acc_nl;jerk_nl;vel_ang_nl;acc_ang_nl;jerk_ang_nl]; % [mean vel, mean_accelerate, mean_jerk]
        scores{nseg}{nl,9} = xpos;  % actual position on the track
    end
end
