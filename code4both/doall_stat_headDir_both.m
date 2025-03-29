% 统计 头部方向与运动方向的夹角
%%
clear
close all
directories_allData_v0_allgood
resuletFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
Dx = {'-cw','-ccw'};
case1 = {'-ontrack','-ontrack_exfg','-ontrack_dsfg'};
ncells = 3;nspkikes = 5;
case2 = sprintf('_%ucell_%uspk',ncells, nspkikes);
midmod = '_cycmid';

nsn = 0;
peak0 = 1; % firing rate threshold to remove place cells
hd_ot = cell(2,5);hdm_ot = cell(2,5);
hd_seq = cell(2,5);hdm_seq = cell(2,5);
for ns = [1,2,4:6,8,10:13,15,16]%1:isession%
    nsn = nsn+1;
    path_ns = path{ns};
    disp(path_ns);
    cd(path_ns);
    outFolder = [resuletFolder path_ns(13:end)];
    file_input1 = 'Data_HeadDirection.mat';
    load(file_input1);
    file_input2 = 'Tseq\Data_angle_ontrack.mat';
    load(file_input2);
    ts = data_angle_all(:,1);
    head_dir = deg2rad(v_dir);

%     xx = cos(data_angle_all(:,2));
%     yy = sin(data_angle_all(:,2));
%     u = cos(rad_head);%x 分量
%     v = sin(rad_head);%y 分量
%     figure
%     quiver(xx,yy,u,v)
    goodphase_ns = Seq_cutPhase{ns,1};%***CW和CCW分割相位一致，用1和2都一样
    case3 = num2str(goodphase_ns);
    % all theta cycle in each lap both Direction
    file_input6=[outFolder,'data_theta_seq_info','_AllLap',case1{1},case2,case3,midmod,'_v5-both.mat'];
    TGa = load(file_input6,...
        'theta_INFO');
    TGa = TGa.theta_INFO;
    nt = 1;
    trackdata_ns = trackdata{ns};load(trackdata_ns,'Ang_RewardLoc_ontrack')
    

    for D =1:2
        disp(['======== process Direction' Dx{D} ' ========'])
        subfolder1 = Directionfolder{D};
        subfolder2 = Phasecutfolder{D};
        
        % get TS_running
        file_input4 = strcat(path_ns,subfolder2,'data_AllLap_thetaphasecut_',case3,'.mat');
        load(file_input4,'TS_running')
        
        if D==1
            m_dir = data_angle_all(:,2)+.5*pi; % 顺时针的运动方向
%             m_dir = mod(m_dir,2*pi);
            m_dir = unwrap(m_dir);
        elseif D==2
            m_dir = data_angle_all(:,2)-.5*pi; % 逆时针的运动方向
%             m_dir = mod(m_dir,2*pi);
            m_dir = unwrap(m_dir);
        end
        delta_rad = unwrap(-head_dir)-m_dir;
        delta_rad = mod(delta_rad,2*pi);
        
        for nl = 1:5
            ind = find (ts>=TS_running(nl,1) & ts<=TS_running(nl,2));
%             delta_rad = unwrap(-head_dir(ind))-m_dir(ind);
%             delta_rad = mod(delta_rad,2*pi);
            hd_ot{D,nl} = [hd_ot{D,nl};delta_rad(ind)];
            hdm_ot{D,nl} = [hdm_ot{D,nl};circ_mean(delta_rad(ind))];
            
            Da = [];Dg = [];
            disp(['======== process lap ' num2str(nl) ' ========'])
            % load good theta cycle in each lap both Direction
            file_input5 =[outFolder,'data_theta_seq_info',case1{1},'_lap',num2str(nl),case2,case3,midmod,'_v5-both.mat'];
            TG = load(file_input5,'ThetaGood');
            TG = TG.ThetaGood;
            indgood = [TG{:,1}];
            Lg = length(indgood);
            Sg = find(diff(indgood)<0);% theta good 中CW和CCW的分界
            
            if ~isempty(Sg)
                Dg(1:Sg) = 1;
                Dg(Sg+1:Lg) = 2;
            else
                cyc_ts = cell2mat(TG(:,3));
                Sg = find(diff(cyc_ts(:,1))>60);
                Dg(1:Sg) = 1;
                Dg(Sg+1:Lg) = 2;
            end
            
            indall = [TGa{nl}{:,1}];
            La = length(indall);
            Sa = find(diff(indall)<0);% theta info 中CW和CCW的分界
            Da(1:Sa) = 1;
            Da(Sa+1:La) = 2;
            
            Num_TGa = size(TGa{nl},1);
            Num_TG = size(TGa{nl},1);
            temp = [];
            for t = find(Dg==D)
                ind = [];
                ind = find (ts>=TG{t,3}(1) & ts<=TG{t,3}(2));
                hd_seq{D,nl} = [hd_seq{D,nl};delta_rad(ind)];
                temp = [temp;delta_rad(ind)];
            end
            hdm_seq{D,nl} = [hdm_seq{D,nl};circ_mean(temp)];
        end
    end
end

figure('Position',[79 315 1756 420])
for nl = 1:5
    subplot(1,5,nl)
    polarhistogram(cell2mat(hd_seq(:,nl)),36)
    title(sprintf('Lap %u',nl))
end
suptitle('difference between head direction and moving direction (in sequence)')

figure('Position',[79 315 1756 420])
for nl = 1:5
    subplot(1,5,nl)
    polarhistogram(cell2mat(hd_ot(:,nl)),36)
    title(sprintf('Lap %u',nl))
end
suptitle('difference between head direction and moving direction (on running)')


figure('Position',[79 315 1756 420])
for nl = 1:5
    subplot(1,5,nl)
    polarhistogram(cell2mat(hdm_seq(:,nl)),12)
    title(sprintf('Lap %u',nl))
end
suptitle('difference between head direction and moving direction (in sequence)')

figure('Position',[79 315 1756 420])
for nl = 1:5
    subplot(1,5,nl)
    polarhistogram(cell2mat(hdm_ot(:,nl)),12)
    title(sprintf('Lap %u',nl))
end
suptitle('difference between head direction and moving direction (on running)')
