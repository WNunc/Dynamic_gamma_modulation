% edited on 4/24/2017

% read video and transform data to angle on the track，把视频的行为学信息转化成角度信息
% save this data for later use
% column 1:4
% time, angle on the track, vel, vel_ang

% updated on 4/27/2017:
% add 'data_angle_all': angle and speed data for the whole session
%2021/03/09 wxl
%以动物第一次出现在第一个点，和，最后一次出现在最后一个点，为开始和结束时间
%pre有10圈
clear
for icon=1:2
    if icon==1
        cd('E:\matlab\bin\code for CT\directories');
        directories_allData_v2
    else
        cd('E:\matlab\bin\code for CT\directories');
        directories_allData_v3 % AD
    end
    file_output = 'Data_angle_ontrack.mat';
    
    for ns = 1:isession
%          for ns = 26:isession
        
        path_ns = path{ns};
        cd(path_ns);
        disp(path_ns);
        
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
        %     vel(vel>=timelimit) = 0.5*(vel(circshift((vel>=timelimit),-3)) + vel(circshift((vel>=timelimit),3)));
        vel(vel>=timelimit) = 0.5*(vel(circshift((vel>=timelimit),-10)) + vel(circshift((vel>=timelimit),10)));%wxl修改，位移3位，仍有大于200的速度
        
        % calculate the angular velocity
        vellimit = 2; % rad, approximately equals to 100cm/s
        velrun = 0.1; % rad, approximately equals to 5cm/s, assuming the rat is running above this threshold
        vel_ang = (posang_ontrack_unwrap(3:end)-posang_ontrack_unwrap(1:end-2))./(t(3:end)-t(1:end-2))';
        %     vel_ang(vel_ang>=vellimit) = 0.5*(vel_ang(circshift((vel_ang>=vellimit),-3)) + vel_ang(circshift((vel_ang>=vellimit),3)));
        vel_ang(vel_ang>=vellimit) = 0.5*(vel_ang(circshift((vel_ang>=vellimit),-10)) + vel_ang(circshift((vel_ang>=vellimit),10)));
        vel_ang = [vel_ang(1);vel_ang;vel_ang(end)];
        
        % calculate the accelerate speed
        acclimit = 120;
        acc = vel2acc(vel,t);
        acc(acc>=acclimit) = 0.5*(acc(circshift((acc>=timelimit),-3)) + vel(circshift((acc>=timelimit),3)));
        
        % divide angle data into each trials
        Ts_start_stop{1,1} = Ts_prerunning{1,1}./ 1000000;
        for i =n_prerunning+1:n_prerunning*2
            Ts_start_stop{1,1}(i,:) = Ts_prerunning{1,2}(i-n_prerunning,:)./ 1000000;
        end
        Ts_start_stop{1,2} = Ts_sample{1,1}(:,[1,3])./ 1000000;
        Ts_start_stop{1,3} = Ts_test{1,1}(:,[1,3])./ 1000000;
        %     Ts_start_stop{1,4} = Ts_posttest{1,1}(:,[1,3])./ 1000000;
        %     nsegment = 4;
        nsegment = 3;
        
        data_angle_all = [t',posang_ontrack',vel,vel_ang,acc];
        data_angle = {};
        data_angle_point = {};
        for nseg = 1:nsegment
            Ts_start_stop0 = Ts_start_stop{nseg};
            nlap = size(Ts_start_stop0,1);
            for nl = 1:nlap
                tstart = Ts_start_stop0(nl,1);
                tstop = Ts_start_stop0(nl,2);
                ind = find(t >= tstart & t <= tstop);
                data_nl = [t(ind)',posang_ontrack(ind)',vel(ind),vel_ang(ind),acc(ind)];
                data_angle{1,nseg}{nl,1} = data_nl;
            end
        end
        
        %以动物第一次出现在第一个点，和，最后一次出现在最后一个点，为开始和结束时间
        for nl = 1:n_prerunning
            ind = min(find(abs(data_angle{1,1}{nl,1}(:,2)-Ang_RewardLoc_ontrack(1))<0.05));
            ind2 = max(find(abs(data_angle{1,1}{nl,1}(:,2)-Ang_RewardLoc_ontrack(18))<0.05));
            ts_start_stop_point{1,1}(nl,1) = data_angle{1,1}{nl,1}(ind,1);%从1开始
            ts_start_stop_point{1,1}(nl,2) = data_angle{1,1}{nl,1}(ind2,1);
            if ts_start_stop_point{1,1}(nl,1)> ts_start_stop_point{1,1}(nl,2)
                fprintf('pre cw day:%8.5f\n',ns) ;
                fprintf('lap:%8.5f\n',nl) ;
                a = ts_start_stop_point{1,1}(nl,1);
               ts_start_stop_point{1,1}(nl,1) = ts_start_stop_point{1,1}(nl,2);
                ts_start_stop_point{1,1}(nl,2) = a;
            end
        end
        
        for nl = n_prerunning+1:n_prerunning*2
            ind = max(find(abs(data_angle{1,1}{nl,1}(:,2)-Ang_RewardLoc_ontrack(1))<0.05));%路过第一个点，结束
            ind2 = min(find(abs(data_angle{1,1}{nl,1}(:,2)-Ang_RewardLoc_ontrack(18))<0.05));%路过第18个点,开始
            if ~isempty(ind2)
                ts_start_stop_point{1,1}(nl,2) = data_angle{1,1}{nl,1}(ind,1);
                ts_start_stop_point{1,1}(nl,1) = data_angle{1,1}{nl,1}(ind2,1);%从18开始
            else
                ts_start_stop_point{1,1}(nl,2) = Ts_start_stop{1,1}(nl,2);
                ts_start_stop_point{1,1}(nl,1) =Ts_start_stop{1,1}(nl,1);
            end
            
            if ts_start_stop_point{1,1}(nl,1)> ts_start_stop_point{1,1}(nl,2)
                fprintf('pre ccw day:%8.5f\n',ns) ;
                fprintf('lap:%8.5f\n',nl) ;
            end
            
        end
        
        n_test= length(sign_correct_test)-numel(find(isnan(ts_test_end)));%一共跑了几次test
        for nl = 1:n_test
            %sample
            ind = min(find(abs(data_angle{1,2}{nl,1}(:,2)-Ang_RewardLoc_ontrack(1))<0.1));
            ind2 = max(find(abs(data_angle{1,2}{nl,1}(:,2)-Ang_RewardLoc_ontrack(18))<0.05));
            ts_start_stop_point{1,2}(nl,1) = data_angle{1,2}{nl,1}(ind,1);
            if ~isempty(ind2)
                ts_start_stop_point{1,2}(nl,2) = data_angle{1,2}{nl,1}(ind2,1);
            else
                ts_start_stop_point{1,2}(nl,2) = Ts_start_stop{1,2}(nl,2);
            end
            if ts_start_stop_point{1,2}(nl,1)> ts_start_stop_point{1,2}(nl,2)
                fprintf('sample day:%8.5f\n',ns) ;
                fprintf('lap:%8.5f\n',nl) ;
                ts_start_stop_point{1,i}(nl,1)= Ts_start_stop{nseg}(nl,1);
                ts_start_stop_point{1,i}(nl,2)= Ts_start_stop{nseg}(nl,2);
            end
            
            %test
            ind = max(find(abs(data_angle{1,3}{nl,1}(:,2)-Ang_RewardLoc_ontrack(1))<0.05));
            ind2 = min(find(abs(data_angle{1,3}{nl,1}(:,2)-Ang_RewardLoc_ontrack(18))<0.05));
            if ~isempty(ind)
                ts_start_stop_point{1,3}(nl,2) = data_angle{1,3}{nl,1}(ind,1);
            else
                ts_start_stop_point{1,3}(nl,2) = Ts_start_stop{1,3}(nl,2);
            end
            ts_start_stop_point{1,3}(nl,1) = data_angle{1,3}{nl,1}(ind2,1);
            if ts_start_stop_point{1,3}(nl,1)> ts_start_stop_point{1,3}(nl,2)
                fprintf('test day:%8.5f\n',ns) ;
                fprintf('lap:%8.5f\n',nl) ;
            end
            
        end
        
        
        %划分每个圈
        for nseg = 1:nsegment
            ts_start_stop0 = ts_start_stop_point{nseg};
            nlap = size(ts_start_stop0,1);
            for nl = 1:nlap
                tstart = ts_start_stop0(nl,1);
                tstop = ts_start_stop0(nl,2);
                ind = find(t >= tstart & t <= tstop);
                data_nl = [t(ind)',posang_ontrack(ind)',vel(ind),vel_ang(ind)];
                data_angle_point{1,nseg}{nl,1} = data_nl;
            end
        end
        
        save(file_output,'data_angle','data_angle_all','data_angle_point','ts_start_stop_point','n_test');
        clear ts_start_stop_point
        cd ../
    end
end