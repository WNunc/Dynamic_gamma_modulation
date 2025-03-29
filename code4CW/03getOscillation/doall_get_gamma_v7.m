% edited on 2021
% add gamma measures to the "scores"
% edited from v6:
% Only get gamma information on tetrodes from csclist but not TTList
% for experiments only with CW trials

clear
directories_allData_v2

TTList0 = 'TTList_dCA1_pyr.txt';
subfix = '-ontrack';% 文件后缀
%%
for ns = 1:isession
    path_ns = path{ns};
    cd(path_ns);
    subfolder = 'Tseq\';
    file_input_cell = 'Cells_singleLap_v2_vel_5.mat';
    file_input_speed = 'Data_angle_ontrack.mat';  % input the speed information
    
    trackdata_ns = trackdata{ns};
    csclist_ns = CSClist_CA1{ns};
    load(trackdata_ns,'diam_inner');
    speedlimit = [5,5/(diam_inner/2)]; %[vel_limit,ang_vel_limit]
    load([subfolder file_input_cell],'spikes');
    load([subfolder file_input_speed],'data_angle_all');
    for nt = 1:1 % trial number
        switch nt
            case 1
                laps=5;
            case 2
                laps=8;
            case 3
                laps=8;
        end
        Scores = {};
        for nl = 1:laps
            fprintf(['%' num2str(nt) '$s'  '  lap-%4$d \n'],'pre', 'samp', 'test', nl);
            file_input = strcat(path{ns},'Tseq\scores',num2str(nt)','-',num2str(nl),subfix,'_v5.mat');
            
            ind_tet = 1;
            load (file_input)
            Scores(nl,:) = scores(nl,:);
            clear scores
            %             switch nt % 'Ts_start_stop','cind_ot'是解码程序生成的有用的变量
            %                 case 1
            %                     save(file_output,'scores','Ts_start_stop','cind_ot');
            %                 case 2
            %                     save(file_output,'scores','Ts_start_stop','Ts_reward','cind_ot');
            %                 case 3
            %                     save(file_output,'scores','Ts_start_stop','Ts_reward','cind_ot');
            %             end
        end
        disp('Getting gamma info');
        [Scores] = bayesGamma_CT_v6(Scores,csclist_ns,TTList0, spikes, ind_tet,data_angle_all,speedlimit,nl);
        for nl = 1:laps
            fprintf(['saving %' num2str(nt) '$s'  '  lap-%4$d \n'],'pre', 'samp', 'test', nl);
            file_output = strcat(path{ns},'Tseq\scores',num2str(nt)','-',num2str(nl),subfix,'_v5.mat');
            load (file_output)
            scores(nl,10:21) = Scores(nl,10:21);
            save(file_output,'scores','Ts_start_stop','cind_ot');
        end
    end
end