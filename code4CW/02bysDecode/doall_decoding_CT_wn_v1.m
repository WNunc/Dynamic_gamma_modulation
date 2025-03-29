% do Bayesian decoding on all running laps in 4 types:
% pre-running, sample trials, test trials, and post-test trials
% modified from doall_decoding_CT_v4: 
% all decoding

clear
directories_allData_v2
%%
minspikes = 1;
mincells = 1;
dt = .04;% 40ms time bin
step = .005;% 5ms step length
peak0 = 1; % firing rate threshold to remove place cells

TTList0 = 'TTList_dCA1_pyr.txt';
subfolder = 'Tseq\';
file_spike_input1 = [subfolder 'Cells_ALLLaps_v2_vel_5.mat'];  % use to get all spikes
file_spike_input2 = [subfolder 'Cells_ALLLaps_v2_vel_5.mat'];  % use itself as the decoder
file_output = strcat('BayesData_CircMap_wn_dt40ms_5ms_',num2str(mincells),'cells_',num2str(minspikes),'.mat');

for ns = 1:isession
    path_ns = path{ns};
    cd(path_ns);
    Ncell = getnumberofcells_cz_v1(TTList0);
    trackdata_ns = trackdata{ns};
    if Ncell < 1
        continue
    end
    disp(pwd)
    for nt = 1:3 % not decoding test segment
        switch nt
            case 1
                laps=6;
            case 2
                laps=8;
            case 3
                laps=8;
            case 4
                laps=6;
        end
        for nl = 1:laps
            % Load Spikes data
            S1 = load(file_spike_input1);  % used to get all spikes
            spikes = S1.spikes;
            Ratemap_Alllaps = S1.Ratemap_AllLaps; % y reflect the trial
            S2 = load(file_spike_input2);  % used to make the decoder
            spikes_vel = S2.spikes;
            Ratemap_Alllaps_vel = S2.Ratemap_AllLaps;
            % Remove the place cells whose peak firing rate<1Hz in all laps
            % excluding prerunning laps
            peak_all = max(S2.Ratemap_AllLaps); % use S2 for detecting
            ind = find(peak_all >= peak0);
            spikes = spikes(ind,:);
            Ratemap_Alllaps = Ratemap_Alllaps(:,ind);
            spikes_vel = spikes_vel(ind,:);
            Ratemap_Alllaps_vel = Ratemap_Alllaps_vel(:,ind);
            % use Ratemap_SingleLaps_vel as decoder
            [scores,Ts_start_stop] = decodewindows_CT_v3(trackdata_ns,spikes,Ratemap_Alllaps_vel,mincells,minspikes,[dt,step],nt,nl);
            
            % save scores
            if ~isempty(scores)
            Scores{nt}(nl,:) = scores{nt}(nl,:);
            end
%             Ts_start_stop=Ts_start_stop{1,1};
%             if nt ~= 1
%                 A(:,1)=Ts_start_stop(:,1);
%                 A(:,2)=Ts_start_stop(:,3);
%                 Ts_reward=Ts_start_stop(:,2);
%                 Ts_start_stop=A;
% %                 file_output=strcat('scores',num2str(nt)','-',num2str(nl),'-wn.mat');
%                 save([subfolder file_output],'scores','Ts_start_stop','Ts_reward');
%             else
%                 file_output=strcat('scores',num2str(nt)','-',num2str(nl),'-wn.mat');
%                 save([subfolder file_output],'scores','Ts_start_stop');
%             end
        end
    end
    save([subfolder file_output],'Scores','Ts_start_stop');
end