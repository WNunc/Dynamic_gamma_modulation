% do Bayesian decoding on all running laps in 4 types:
% pre-running, sample trials, test trials, and post-test trials
% modified from v4:
% this script should be run after doall_decoding_CT_v4.m having run

clear
disp('++++++start singlelap decoding++++++')
disp('downsample same number of slow gamma phase locking cell')
directories_allData_v2
%%
minspikes = 1;
mincells = 1;
dt = .04;% 40ms time bin
step = .005;% 5ms step length
peak0 = 1; % firing rate threshold to remove place cells

TTList0 = 'TTList_dCA1_pyr.txt';
subfolder = 'Tseq\';
file_spike_input1 = [subfolder 'Cells_allsegment_v1_vel_5.mat'];  % use to get all spikes
file_spike_input2 = [subfolder 'Cells_singleLap_v2_vel_5.mat'];  % use itself as the decoder
file_input3 = [subfolder 'data_phaselocking_spkmin.mat'];  % phase locking feature
% file_input3 = [subfolder 'data_phaselocking_TSlap.mat'];  % ** phase locking feature trial TS line 109
fprintf(1,'all spike:\t%s\ndecoder:\t%s\nphase lock:\t%s\n',file_spike_input1,file_spike_input2,file_input3)
subfix = '-ontrack_dssg';%
for ns = 1:isession
    path_ns = path{ns};
    cd(path_ns);
    Ncell = getnumberofcells_cz_v1(TTList0);
    trackdata_ns = trackdata{ns};
    load(file_input3,'plFeature_singlap','plFeature_alllap','plFeature')
    
    if Ncell < 1
        continue
    end
    disp(pwd)
    for nt = 1% this file for pre_running only
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
        Ratemap_seg = S1.Ratemap_seg{nt}; % y reflect the trial
        % Remove the place cells whose peak firing rate<1Hz in all laps
        peak_all = max(S1.Ratemap_seg{nt}); % use S2 for detecting
        ind = find(peak_all >= peak0);
        spikes = spikes(ind,:);
        % find cell on track as detecting
        X = cell2struct(S1.fieldProp_seg{nt}(ind),'placefield',1);
        load(trackdata_ns,'Ang_RewardLoc_ontrack')
        COM = [];
        for ncell  = 1:length(ind)
            COM(ncell) = X(ncell).placefield(1).x_COM;
        end
        
        cind_ot = find(COM>Ang_RewardLoc_ontrack(1) & COM<Ang_RewardLoc_ontrack(18));
        cellnumontrack = length(cind_ot);
        sgamma_rayleighP = [plFeature_alllap{1,2}.rayleighP];
        sgamma_rayleighP_ontrack = sgamma_rayleighP(ind);
        sgamma_rayleighP_ontrack = sgamma_rayleighP_ontrack(cind_ot);
        % find(sgamma_rayleighP_ontrack>0&sgamma_rayleighP_ontrack<0.05)
        cind_ot_exsg = find(sgamma_rayleighP_ontrack>=0.05|sgamma_rayleighP_ontrack==0);
        cellnumontrack_exsg = length(cind_ot_exsg);
        Ncell_sg = cellnumontrack-cellnumontrack_exsg;
        fprintf(1,'pre-running :\t%u cells ontrack and %u cell phase locking in sgamma\n',...
            cellnumontrack,Ncell_sg)
        S2 = load(file_spike_input2);  % used to make the decoder
        cind_ot_exsg = cind_ot(cind_ot_exsg);
        cind_ot_sg = cind_ot(~ismember(cind_ot,cind_ot_exsg))      
        if Ncell_sg == 0
            fprintf(1,'on slow gamma phase in this session (%s)\n', path_ns);
            continue
        else
            if cellnumontrack_exsg<Ncell_sg
                cind_dsamp_ex = cind_ot_exsg
                cind_dsamp = cind_ot_sg
                disp('causion: slow gamma phase locking cell > non phase locking cell')
            else
            rand('seed',ns)
            tempind_dsamp = sort(randperm(cellnumontrack_exsg,Ncell_sg),'ascend');%随机生成和相锁神经元数目相同的序列
            cind_dsamp_ex = cind_ot_exsg(tempind_dsamp)%被排除的非相锁神经元
            cind_dsamp = cind_ot(~ismember(cind_ot,cind_dsamp_ex));%排除这个神经元后剩下的用于detecting的
        end
        
        for nl = 1:laps
            Ratemap_SingleLaps_vel = S2.Ratemap_singlelap{nl,nt};
            Ratemap_SingleLaps_vel = Ratemap_SingleLaps_vel(:,ind);
            % use Ratemap_Ratemap_seg_vel as decoder
            [scores1,Ts_start_stop] = decodewindows_CT_v3(trackdata_ns,spikes(cind_dsamp,:),Ratemap_SingleLaps_vel(:,cind_dsamp),mincells,minspikes,[dt,step],nt,nl);%
            % save scores
            % scores = scores{1,nt};
            if ~isempty(scores1)
                file_old=strcat('scores',num2str(nt)','-',num2str(nl),'-ontrack_v5.mat');
                load ([subfolder file_old],'scores')
                scores(:,1:9) = scores1{1,nt};
            end
            
            Ts_start_stop=Ts_start_stop{1,1};
            if nt ~= 1
                A(:,1)=Ts_start_stop(:,1);
                A(:,2)=Ts_start_stop(:,3);
                Ts_reward=Ts_start_stop(:,2);
                Ts_start_stop=A;
                file_output=strcat('scores',num2str(nt)','-',num2str(nl),subfix,'.mat');
%                 file_output=strcat('scores',num2str(nt)','-',num2str(nl),subfix,'_TSlap.mat');% **
                save([subfolder file_output],'scores','Ts_start_stop','Ts_reward','cind_ot','cind_dsamp');
            else
                file_output=strcat('scores',num2str(nt)','-',num2str(nl),subfix,'.mat');
%                 file_output=strcat('scores',num2str(nt)','-',num2str(nl),subfix,'_TSlap.mat');% **
                save([subfolder file_output],'scores','Ts_start_stop','cind_ot','cind_dsamp');
            end
        end
    end
end