% edited on 4/24/2017

% read video and transform data to angle on the track，把视频的行为学信息转化成角度信息
% save this data for later use
% column 1:4
% time, angle on the track, vel, vel_ang

% updated on 4/27/2017: 
% add 'data_angle_all': angle and speed data for the whole session

clear

directories_allData_v2

file_output = 'Tseq\Data_angle_ontrack.mat';

for ns = 1:isession
   
    path_ns = path{ns};
    cd(path_ns);
    disp(path_ns);
    mkdir('Tseq')
    trackdata_ns = trackdata{ns};
    load(trackdata_ns)
    
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
    
    
    % divide angle data into each trials
    Ts_start_stop{1,1} = Ts_prerunning{1,1}./ 1000000;
    Ts_start_stop{1,2} = Ts_sample{1,1}(:,[1,2])./ 1000000;
    Ts_start_stop{1,3} = Ts_test{1,1}(:,[1,2])./ 1000000;
%     Ts_start_stop{1,4} = Ts_posttest{1,1}(:,[1,3])./ 1000000;
%     nsegment = 4;
    nsegment = 3;
    
    data_angle_all = [t',posang_ontrack',vel,vel_ang];
    
    data_angle = {};
    for nseg = 1:nsegment
        Ts_start_stop0 = Ts_start_stop{nseg};
        nlap = size(Ts_start_stop0,1);
        for nl = 1:nlap
            tstart = Ts_start_stop0(nl,1);
            tstop = Ts_start_stop0(nl,2);
            ind = find(t >= tstart & t <= tstop);
            data_nl = [t(ind)',posang_ontrack(ind)',vel(ind),vel_ang(ind)];
            data_angle{1,nseg}{nl,1} = data_nl;
        end
    end
    
    save(file_output,'data_angle','data_angle_all');
    cd ../
end