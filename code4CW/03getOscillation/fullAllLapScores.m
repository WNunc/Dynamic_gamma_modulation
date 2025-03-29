

clear
directories_allData_v2
subfix = '-ontrack';
%% 把单圈解码得到的scores中的节律信息填到全部圈解码的scores中
for ns = 1:isession
    path_ns = path{ns};
    cd(path_ns);
    file_output = 'Tseq\BayesData_CircMap_wn_dt40ms_5ms_1cells_1.mat';
    load (file_output)
    for nt = 1:1 % trial number
        switch nt
            case 1
                laps=5;
            case 2
                laps=8;
            case 3
                laps=8;
        end
        for nl = 1:laps
            fprintf(['%' num2str(nt) '$s'  '  lap-%4$d \n'],'pre', 'samp', 'test', nl);
            dir_input = strcat(path{ns},'Tseq\scores',num2str(nt)','-',num2str(nl),subfix,'_v5.mat');
            file_input = dir_input;
            load (file_input,'scores')
            Scores{nt}(nl,10:21) = scores(nl,10:21); 
        end
    end
    save(file_output,'Scores','Ts_start_stop')
    cd ..\
end

%% 当不小心算错了scores，或者是新的条件下得到了scores，而又不希望浪费时间来运行滤波时使用
% for ns = [4 5 6 7 8]%:isession
%     path_ns = path{ns};
%     cd(path_ns);
%     
%     
%     for nt = 1:1 % trial number
%         switch nt
%             case 1
%                 laps=5;
%             case 2
%                 laps=8;
%             case 3
%                 laps=8;
%         end
%         for nl = 1:laps
%             file_output = strcat(path{ns},'Tseq\scores',num2str(nt),'-',num2str(nl),'-ontrack-wn.mat');
%             load (file_output)
%             Scores = scores;
%             fprintf(['%' num2str(nt) '$s'  '  lap-%4$d \n'],'pre', 'samp', 'test', nl);
%             dir_input = strcat(path{ns},'Tseq\scores',num2str(nt),'-',num2str(nl),'-wn.mat');
%             file_input = dir_input;
%             load (file_input,'scores')
%             Scores{nt}(nl,10:21) = scores(nl,10:21); 
%             scores = Scores{1};
%             save(file_output,'scores','Ts_start_stop')
%         end
%     end
%     
%     cd ..\
% end