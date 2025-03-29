% do Bayesian decoding on all running laps in 4 types:
% pre-running, sample trials, test trials, and post-test trials
% modified from v3:

clear
disp('++++++start singlelap decoding++++++')
directories_allData_v2
%%
minspikes = 1;
mincells = 1;
dt = .04;% 40ms time bin
step = .005;% 10ms step length
peak0 = 1; % firing rate threshold to remove place cells

TTList0 = 'TTList_dCA1_pyr.txt';
subfolder = 'Tseq\';
file_spike_input1 = [subfolder 'Cells_singleLap_v2_vel_5.mat'];  % use to get all spikes
file_spike_input2 = [subfolder 'Cells_singleLap_v2_vel_5.mat'];  % use itself as the decoder
subfix = '-ontrack';%
for ns = 2:isession
    path_ns = path{ns};
    cd(path_ns);
    Ncell = getnumberofcells_cz_v1(TTList0);
    trackdata_ns = trackdata{ns};
    if Ncell < 1
        continue
    end
    disp(pwd)
    
    
    for nt = 1
        switch nt
            case 1
                laps=5;
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
            Ratemap_SingleLaps = S1.Ratemap_singlelap{nl,nt}; % y reflect the trial
            S2 = load(file_spike_input2);  % used to make the decoder
            spikes_vel = S2.spikes;
            Ratemap_SingleLaps_vel = S2.Ratemap_singlelap{nl,nt};
            % Remove the place cells whose peak firing rate<1Hz in all laps
            % excluding prerunning laps
            peak_all = max(S2.Ratemap_singlelap{nl,nt}); % use S2 for detecting
            ind = find(peak_all >= peak0);
            spikes = spikes(ind,:);
            Ratemap_SingleLaps = Ratemap_SingleLaps(:,ind);
            spikes_vel = spikes_vel(ind,:);
            Ratemap_SingleLaps_vel = Ratemap_SingleLaps_vel(:,ind);
              
            % find cell on track
            X = cell2struct(S2.fieldProp_singlelap{nl,nt}(ind),'placefield',1);
            load(trackdata_ns,'Ang_RewardLoc_ontrack')
            COM = [];
            for ncell  = 1:length(ind)
                COM(ncell) = X(ncell).placefield(1).x_COM;
            end
            
            cind_ot = find(COM>Ang_RewardLoc_ontrack(1) & COM<Ang_RewardLoc_ontrack(18));
            cellnumontrack = length(cind_ot)
            % use Ratemap_SingleLaps_vel as decoder
            [scores1,Ts_start_stop] = decodewindows_CT_v3(trackdata_ns,spikes(cind_ot,:),Ratemap_SingleLaps_vel(:,cind_ot),mincells,minspikes,[dt,step],nt,nl);%
            
            % save scores
%             scores = scores{1,nt};
            if ~isempty(scores1)
                file_output=strcat('scores',num2str(nt)','-',num2str(nl),subfix,'-wn.mat');
                load ([subfolder file_output],'scores')
                scores(:,1:9) = scores1{1,nt};
            end
            
            Ts_start_stop=Ts_start_stop{1,1};
            if nt ~= 1
                A(:,1)=Ts_start_stop(:,1);
                A(:,2)=Ts_start_stop(:,3);
                Ts_reward=Ts_start_stop(:,2);
                Ts_start_stop=A;
                file_output=strcat('scores',num2str(nt)','-',num2str(nl),subfix,'-wn.mat');
                save([subfolder file_output],'scores','Ts_start_stop','Ts_reward');
            else
                file_output=strcat('scores',num2str(nt)','-',num2str(nl),subfix,'-wn.mat');
                save([subfolder file_output],'scores','Ts_start_stop');
            end
        end
    end
end