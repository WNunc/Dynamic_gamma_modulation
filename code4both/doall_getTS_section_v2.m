% 获取ontrack时,fastgamma domain 的每个loc的时间点
% 根据这些时间点，将整个奔跑路径分为3段
% 前期中期和后期，即6个点一个阶段
% 数据格式为 3个阶段*2个时间戳（阶段开始，阶段结束）


clear
close all
directories_allData_v0
resuletFolder = 'H:\neuralynx\gamma in sequence result';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
case1 = '-ontrack';
ncells = 3;nspkikes = 5;
case2 = sprintf('_%ucell_%uspk',ncells, nspkikes);
midmod = '_cycmid';
TPsweep = 7;
SPsweep = 10;
peak0 = 1;

for ns = 13%1:isession
    path_ns = path{ns};
    cd(path_ns)
    disp(path_ns);cd(path_ns);
    trackdata_ns = trackdata{ns};
    load(trackdata_ns)
    nseg = 1;
    
    outFolder = [resuletFolder path_ns(13:end)];
    if ~exist(outFolder,'dir')
        mkdir([resuletFolder path_ns(13:end)])
    end
    
    RunSpeed = cell(5,1);
    TS_section = cell(3,2);
    TS_section_temp = [];
    for nl = 1:5
        for D = 1:2
            goodphase_ns = Seq_cutPhase{ns,D};%***
            subfolder1 = Directionfolder{D};
            subfolder2 = Phasecutfolder{D};
            case3 = num2str(goodphase_ns);
            load([path{ns} subfolder1 'Data_angle_ontrack.mat'])
            
            file_input1 = strcat(path_ns,subfolder1,'scores',num2str(nseg),'-',num2str(nl),case1,'_v5.mat');
            load(file_input1)
            file_input2 = strcat(path_ns,subfolder2, 'data_AllLap_thetaphasecut_',case3,'.mat');
            load(file_input2,'TS_running')
            fileinput3 = [subfolder1 'Cells_allsegment_v1_vel_0.mat'];
            S1 = load(fileinput3);
            mapAxis = S1.mapAxis;
            lap_ts = scores{nl,6};%行为学时间
            lap_pos = scores{nl,5}(scores{nl,4});%实际位置
            lap_speed = scores{nl,8}(1,:);
            lap_ts = lap_ts(1:end-20);
            lap_pos = lap_pos(1:end-20);
            lap_speed = lap_speed(1:end-20);
            
            if D == 1
                [TSlimit, speed, speed_ang]= TSfind_laplimit_CW(data_angle,mapAxis,lap_speed,lap_pos,lap_ts,nseg,nl);
            else
                [TSlimit, speed, speed_ang]= TSfind_laplimit_CCW(data_angle,mapAxis,lap_speed,lap_pos,lap_ts,nseg,nl);
            end
            
            TS_section_temp(nl,1,D) = TS_running(nl,1);
            TS_section_temp(nl,4,D) = TS_running(nl,2);
            TS_section_temp(nl,[2,3],D) = TSlimit;
            TS_section{1,D} = TS_section_temp(:,1:2,D);
            TS_section{2,D} = TS_section_temp(:,2:3,D);
            TS_section{3,D} = TS_section_temp(:,3:4,D);
        end
    end
%     save([outFolder,'TS_section_v2.mat'],'TS_section')
end


function [TSlimit, speed, speed_ang]= TSfind_laplimit_CW(data_angle,Ang_RewardLoc_ontrack,lap_speed,lap_pos,lap_ts, nsegment,nlap)
% input:
% data_angle ---- from file named Data_angle_ontrack.mat
% Ang_RewardLoc_ontrack ---- from file named like date_CT_tracking.mat
% output:
% TSlimit ---- the start and the end timgstamp of one lap
% speed ---- raw speed from data_angle
% speed_ang ---- raw angle speed from data_angle
speed = data_angle{nsegment}{nlap}(:,3);% cm/s
speed_ang = data_angle{nsegment}{nlap}(:,4);% rad/s
ind_ts_1 = max(find(lap_speed'>=10 & lap_pos<=Ang_RewardLoc_ontrack(17)));
ind_ts_2 = max(find(lap_speed'>=10 & lap_pos<=Ang_RewardLoc_ontrack(52)));
TSlimit = [lap_ts(ind_ts_1), lap_ts(ind_ts_2)];
end

function [TSlimit, speed, speed_ang]= TSfind_laplimit_CCW(data_angle,Ang_RewardLoc_ontrack,lap_speed,lap_pos,lap_ts, nsegment,nlap)
% input:
% data_angle ---- from file named Data_angle_ontrack.mat
% Ang_RewardLoc_ontrack ---- from file named like date_CT_tracking.mat
% output:
% TSlimit ---- the start and the end timgstamp of one lap
% speed ---- raw speed from data_angle
% speed_ang ---- raw angle speed from data_angle
speed = data_angle{nsegment}{nlap}(:,3);% cm/s
speed_ang = data_angle{nsegment}{nlap}(:,4);% rad/s
ind_ts_1 = max(find(lap_speed'>=10 & lap_pos>=Ang_RewardLoc_ontrack(52)));
ind_ts_2 = max(find(lap_speed'>=10 & lap_pos>=Ang_RewardLoc_ontrack(17)));
TSlimit = [lap_ts(ind_ts_1), lap_ts(ind_ts_2)];
end
