% 统计 spk数量
% 在奔跑时的spk = spk_ot
% 在sequence中的spk = spk_seq

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

file_input5 = [];% good theta cycle in each lap
file_input6 = [];% all theta cycle in each lap

lockat = {'firstlap','alllap','f2lap'};
L = 2;
SPK_ot_fg = [];SPK_ot_nfg = [];
SPK_seq_fg  = [];SPK_seq_nfg  = [];
nsn = 0;
peak0 = 1; % firing rate threshold to remove place cells
for ns = [1,2,4:6,8,10:13,15,16]%1:isession%
    nsn = nsn+1;
    path_ns = path{ns};
    disp(path_ns);
    cd(path_ns);
    outFolder = [resuletFolder path_ns(13:end)];
    
    pause(0.1)
    goodphase_ns = Seq_cutPhase{ns,1};%***CW和CCW分割相位一致，用1和2都一样
    case3 = num2str(goodphase_ns);
    % all theta cycle in each lap both Direction
    file_input6=[outFolder,'data_theta_seq_info','_AllLap',case1{1},case2,case3,midmod,'_v5-both.mat'];
    TGa = load(file_input6,...
        'theta_INFO');
    TGa = TGa.theta_INFO;
    nt = 1;
    trackdata_ns = trackdata{ns};load(trackdata_ns,'Ang_RewardLoc_ontrack')
    
    spk_ot_fg = zeros(2,5); spk_ot_nfg = zeros(2,5);
    spk_seq_fg = zeros(2,5); spk_seq_nfg = zeros(2,5);
    for D =1:2
        disp(['======== process Direction' Dx{D} ' ========'])
        subfolder1 = Directionfolder{D};
        subfolder2 = Phasecutfolder{D};
        % use to get all spikes
        file_input1 = strcat(path_ns,subfolder1, 'Cells_allsegment_v1_vel_5.mat');
        % Load Spikes data
        S1 = load(file_input1);  % used to get all spikes
        spikes = S1.spikes;
        Ratemap_seg = S1.Ratemap_seg{nt};
        % Remove the place cells whose peak firing rate<1Hz in all laps
        peak_all = max(Ratemap_seg); % use S2 for detecting
        ind = find(peak_all >= peak0);
        spikes = spikes(ind,:);
        
        % get cellind
        file_input2 = [subfolder1,'scores1-1-ontrack_dsfg_TSlap_vel0',lockat{L},'v3.mat'];
        load(file_input2,'cind_ot_nonfg');
        file_input3 = [subfolder1,'scores1-1-ontrack_exfg_TSlap_vel0',lockat{L},'v2.mat'];
        load(file_input3,'cind_ot_fg');
        
        % get TS_running
        file_input4 = strcat(path_ns,subfolder2,'data_AllLap_thetaphasecut_',case3,'.mat');
        load(file_input4,'TS_running')
        
        for nl = 1:5
            for ic = 1:length(cind_ot_fg)
                spk = spikes{cind_ot_fg(ic),1};
                ind = find (spk>=TS_running(nl,1) & spk<=TS_running(nl,2));
                spk_ot_fg(D,nl) = spk_ot_fg(D,nl)+length(ind);
            end
            for ic = 1:length(cind_ot_nonfg)
                spk = spikes{cind_ot_nonfg(ic),1};
                ind = find (spk>=TS_running(nl,1) & spk<=TS_running(nl,2));
                spk_ot_nfg(D,nl) = spk_ot_nfg(D,nl)+length(ind);
            end
            
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
            
            for t = find(Dg==D)
                
                for ic = 1:length(cind_ot_fg)
                    spk = spikes{cind_ot_fg(ic),1};
                    ind = find (spk>=TG{t,3}(1) & spk<=TG{t,3}(2));
                    spk_seq_fg(D,nl) = spk_seq_fg(D,nl)+length(ind);
                end
                for ic = 1:length(cind_ot_nonfg)
                    spk = spikes{cind_ot_nonfg(ic),1};
                    ind = find (spk>=TG{t,3}(1) & spk<=TG{t,3}(2));
                    spk_seq_nfg(D,nl) = spk_seq_nfg(D,nl)+length(ind);
                end
            end
        end
    end
    SPK_ot_fg(:,:,nsn) = spk_ot_fg;
    SPK_ot_nfg(:,:,nsn) = spk_ot_nfg;
    SPK_seq_fg(:,:,nsn) = spk_seq_fg;
    SPK_seq_nfg(:,:,nsn) = spk_seq_nfg;
end
%%
nSPK_fg = squeeze(sum(SPK_ot_fg,[1,2]));
nSPK_nfg = squeeze(sum(SPK_ot_nfg,[1,2]));

spkN = [nSPK_fg,nSPK_nfg];
plot_spkN(spkN)
ylim([0,2000])
title('on running')
saveas(gcf,'C:\Users\tju\OneDrive\gamma整理\投稿2\comment\cell characteristic\spkN on running.png')
saveas(gcf,'C:\Users\tju\OneDrive\gamma整理\投稿2\comment\cell characteristic\spkN on running.eps','epsc')
nSPK_fg = squeeze(sum(SPK_seq_fg,[1,2]));
nSPK_nfg = squeeze(sum(SPK_seq_nfg,[1,2]));
spkN = [nSPK_fg,nSPK_nfg];
plot_spkN(spkN)
ylim([0,1000])
title('in sequence')
saveas(gcf,'C:\Users\tju\OneDrive\gamma整理\投稿2\comment\cell characteristic\spkN in sequence.png')
saveas(gcf,'C:\Users\tju\OneDrive\gamma整理\投稿2\comment\cell characteristic\spkN in sequence.eps','epsc')
function plot_spkN(spkN)
color = {'#F24444','#F2CA50'};
figure('Position', [1179 316 283 278])
for x = 1:2
    [bar,er] = barwitherror(x,spkN(:,x));
    bar.FaceColor = color{x};
    hold on
end
hold off
xticks([1,2])
xticklabels({'fg cell','nfg cell'})
set(gca,'FontSize',14)
ax = gca;
% ax.YAxis.Exponent = 3;
ylabel('Spike number')
% saveas(gcf,'spkN_insession','epsc')
% saveas(gcf,'spkN_insession.png')
end
    