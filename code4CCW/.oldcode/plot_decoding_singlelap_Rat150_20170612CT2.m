% edited on 03/24/2017
% updated on 04/26/2017: add group_GammaPow_sampletest_v2

% Plot the Baysian decoding for each trial, including:
% pre-running, sample trials, test trials, and post-test trials
% modified on plot_decoding_eachtrial_batch.m for old paradigm
% this code is for new paradigm from Rat139

clear
% (1) use the overall place field as d decoder
% file_input = 'BayesData_CircMap_v1_dt40ms_10ms_1cells_1spk.mat';
% fig_folder ='Figure_decoding_singlelap_v1\';

% (2) use the overall place field excluding pre-ruuning trials as a decoder
file_input = 'BayesData_CircMap_v1_dt40ms_10ms_1cells_1spk_v2.mat';
file_input_speed = 'Data_angle_ontrack.mat';  % input the speed information
fig_folder ='Figure_decoding_singlelap_40ms_10ms_1cells_1spk_v2\';

% file_input = 'BayesData_CircMap_v1_dt20ms_10ms_1cells_1spk_v2.mat';
% file_input_speed = 'Data_angle_ontrack.mat';  % input the speed information
% fig_folder ='Figure_decoding_singlelap_20ms_10ms_1cells_1spk_v2\';


directories_allData_v1

file_input_thetawin = 'thetawindows_minspk_v2.mat';

sign_plot = 1;   % 1 = plot; 0 = not plot
numcell_seq = 2;  % use the theta sequence with at least 3 cells firing
fg_color = [238 59 59]./255;
sg_color = [0 0 205]./255;
dt = 0.04; % 40 ms window
Nlaps = 0;
ncell_threshold = 30; % Only use data with more than 30 cells

fig_dir = 'G:\MATLAB\MainFunctions\Code & Data Circular Track\';
TTList0 = 'TTList_dCA1_pyr.txt';

cd(fig_dir)
if sign_plot == 1
    mkdir(fig_dir,fig_folder);
    delete(strcat(fig_folder,'*.*'));
    ffa = figure('Units','normalized','Position',[0 0 1 1]);
end

for ns = 21
    path_ns = path{ns};
    cd(path_ns);
    
    Ncell = getnumberofcells_cz_v1(TTList0);
    
    if exist(file_input,'file')>0 & Ncell >= ncell_threshold
        disp(pwd)
        load (file_input);
        load (file_input_thetawin);
        load(file_input_speed,'data_angle');
        
        trackdata_ns = trackdata{ns};
        load(trackdata_ns);
        
        % Only plot sample-test laps for now
        scores_sample = scores{1,2};
        scores_test = scores{1,3};
        nlaps = size(scores_sample,1);
        
        % Limit the time period between RewardLoc(1) and stop location
        ts_start_stop_sample = Ts_sample{1}(:,[1,2])./ 1000000;
        ts_start_stop_test = Ts_test{1}(:,[1,2])./ 1000000;
        
        ang_vel_limit = 5/(Diam_inner/2); % 5cm/s
        
        for nl = 6
            Nlaps = Nlaps+1;
            ts_start_stop_sample_nl = ts_start_stop_sample(nl,:);
            ts_start_stop_test_nl = ts_start_stop_test(nl,:);
            
            %% start plotting when the rat passed reward location #1
            ang_reward1 = Ang_RewardLoc_ontrack(1);
            [~,reward_angle_bin] = min(abs(ang_reward1-scores_sample{nl,5}));
            
            pos_angle_bin = scores_sample{nl,4};
            pos_angle_bin_2 = [pos_angle_bin(1:end-1)',pos_angle_bin(2:end)'];
            ind = find(pos_angle_bin_2(:,1) <= reward_angle_bin & pos_angle_bin_2(:,2) > reward_angle_bin & abs(diff(pos_angle_bin_2'))'<5);
            ind = ind(end);
            ts_start_stop_sample_nl(1) = scores_sample{nl,6}(ind);
            
            pos_angle_bin = scores_test{nl,4};
            pos_angle_bin_2 = [pos_angle_bin(1:end-1)',pos_angle_bin(2:end)'];
            ind = find(pos_angle_bin_2(:,1) <= reward_angle_bin & pos_angle_bin_2(:,2) > reward_angle_bin & abs(diff(pos_angle_bin_2'))'<5);
            ind = ind(end);
            ts_start_stop_test_nl(1) = scores_test{nl,6}(ind);
            
            %% stop plotting when the rat passed correct/incorrect reward location
            if ~isnan(Sign_correct_sample{1}(nl))
                ang_reward_sample_border = Ang_RewardLoc_ontrack_zone_narrow(Ind_rewardloc_sample{1}(nl),1);
            else
                ang_reward_sample_border = Ang_RewardLoc_ontrack_zone_narrow(Ind_rewardloc,1);
            end
            [~,reward_angle_bin] = min(abs(ang_reward_sample_border-scores_sample{nl,5}));
            
            pos_angle_bin = scores_sample{nl,4};
            pos_angle_bin_2 = [pos_angle_bin(1:end-1)',pos_angle_bin(2:end)'];
            ind = find(pos_angle_bin_2(:,1) <= reward_angle_bin & pos_angle_bin_2(:,2) > reward_angle_bin & abs(diff(pos_angle_bin_2'))'<5);
            ind = ind(1);
            ts_start_stop_sample_nl(2) = min(ts_start_stop_sample_nl(2),scores_sample{nl,6}(ind));
            
            if ~isnan(Sign_correct_test{1}(nl))
                ang_reward_test = Ang_RewardLoc_ontrack_zone_narrow(Ind_rewardloc_test{1}(nl),1);
            else
                ang_reward_test = Ang_RewardLoc_ontrack_zone_narrow(Ind_rewardloc,1);
            end
            [~,reward_angle_bin] = min(abs(ang_reward_test-scores_test{nl,5}));
            
            pos_angle_bin = scores_test{nl,4};
            pos_angle_bin_2 = [pos_angle_bin(1:end-1)',pos_angle_bin(2:end)'];
            ind = find(pos_angle_bin_2(:,1) <= reward_angle_bin & pos_angle_bin_2(:,2) > reward_angle_bin & abs(diff(pos_angle_bin_2'))'<5);
            ind = ind(1);
            ts_start_stop_test_nl(2) = min(ts_start_stop_test_nl(2),scores_test{nl,6}(ind));
            
            %% plot in each sample-test lap
            if sign_plot == 1
                if nl == 1
                    ts_start_stop_sample_nl(1) = ts_start_stop_sample_nl(1)+1.5;
                elseif nl == 5
                    ts_start_stop_sample_nl(1) = ts_start_stop_sample_nl(1)+2;
                end
                
                figure(ffa)
                plot_decoding_sampletest_v1;
                
                ind = find(path_ns == '\');
                day = path_ns(ind(4)+1:ind(5)-1);
                bmpImage = sprintf('%s%s%s%s%s%s%s',fig_dir,fig_folder,day,'-trial',num2str(nl),'.bmp');
                f = getframe(gcf);
                [pic, cmap] = frame2im(f);
                imwrite(pic,bmpImage,'bmp');
                % clf
            end
            
        end
    end
    
    cd ../
end
close all
