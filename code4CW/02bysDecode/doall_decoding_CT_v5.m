% do Bayesian decoding on all running laps in 4 types:
% pre-running, sample trials, test trials, and post-test trials
% modified from v3:

clear
disp('++++++start singlelap decoding++++++')
% directories_allData_v2
directories_allData_v0
%%
minspikes = 1;
mincells = 1;
dt = .04;% 40ms time bin
step = .005;% 10ms step length
peak0 = 1; % firing rate threshold to remove place cells

TTList0 = 'TTList_dCA1_pyr.txt';
subfolder = 'Tseq\';
file_spike_input1 = [subfolder 'Cells_allsegment_v1_vel_0.mat'];  % use to get all spikes
file_spike_input2 = [subfolder 'Cells_allsegment_v1_vel_5.mat'];  % use to find cell on the track
file_spike_input3 = [subfolder 'Cells_singleLap_v2_vel_5.mat'];  % use single lap as the decoder
subfix = '-ontrack';%
for ns = 1:isession
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
        % Load Spikes data
        S1 = load(file_spike_input1);  % used to get all spikes
        spikes = S1.spikes;
        S2 = load(file_spike_input2);  % used to get all cell on the track
        Ratemap_seg = S2.Ratemap_seg{nt}; 
        % Remove the place cells whose peak firing rate<1Hz in all laps
        peak_all = max(Ratemap_seg); 
        ind = find(peak_all >= peak0);
        spikes = spikes(ind,:);
        % find cell on track as detecting
        X = cell2struct(S2.fieldProp_seg{nt}(ind),'placefield',1);
        load(trackdata_ns,'Ang_RewardLoc_ontrack')
        COM = [];
        for ncell  = 1:length(ind)
            COM(ncell) = X(ncell).placefield(1).x_COM;
        end
        cind_ot = find(COM>Ang_RewardLoc_ontrack(1) & COM<Ang_RewardLoc_ontrack(18));
        cellnumontrack = length(cind_ot);
        
        S3 = load(file_spike_input3);  % used to make the decoder
        
        for nl = 1:laps
            Ratemap_SingleLaps_vel = S3.Ratemap_singlelap{nl,nt};
            Ratemap_SingleLaps_vel = Ratemap_SingleLaps_vel(:,ind);
            % use Ratemap_Ratemap_seg_vel as decoder
            [scores1,Ts_start_stop] = decodewindows_CT_v3(trackdata_ns,spikes(cind_ot,:),Ratemap_SingleLaps_vel(:,cind_ot),mincells,minspikes,[dt,step],nt,nl);%
            
            % save scores
            if ~isempty(scores1)
%                 scores = scores1{1,nt};
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                file_old=strcat('scores',num2str(nt)','-',num2str(nl),subfix,'_v5.mat');
                load ([subfolder file_old],'scores')
                scores(:,1:9) = scores1{1,nt};
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
            Ts_start_stop=Ts_start_stop{1,1};
            if nt ~= 1
                A(:,1)=Ts_start_stop(:,1);
                A(:,2)=Ts_start_stop(:,3);
                Ts_reward=Ts_start_stop(:,2);
                Ts_start_stop=A;
                file_output=strcat('scores',num2str(nt)','-',num2str(nl),subfix,'_v5.mat');
                save([subfolder file_output],'scores','Ts_start_stop','Ts_reward','cind_ot');
            else
                file_output=strcat('scores',num2str(nt)','-',num2str(nl),subfix,'_v5.mat');
                save([subfolder file_output],'scores','Ts_start_stop','cind_ot');
            end
        end
    end
end