% 只为了画example的图
% 把每个cell在fast gamma或者slow gamma episode时的放电画出来
% 根据doall_phaselocking_wn_v2修改

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% origin file is phaselocking_feature_catch.m
%
% this file can read those result *.mat files
% extract and filt the LFP signal from *.ncs
% calculate the phase of LFP the spike occured
% and PLV
%
% output each folder PLV result and LFP phase when spikes occured
% spikePhase  cell numbers × frequency bands
%                       every cell array is 1 × the number of spikes
% spikePhase_alllap  cell numbers × frequency bands × segment numbers
%                       every cell array is 1 × the number of spikes
% spikePhase_singlap  cell numbers × frequency bands × lap numbers
%                       every cell array is 1 × the number of spikes
% plFeature  1 × frequency bands
% plFeature_alllap  segment numbers × frequency bands
% plFeature_singlap  lap numbers × frequency bands
%
% created by WN on 20220412
%
% ==========use timestamp at trial start and end==========

%%
clear
close all
% directories_allData_v2
directories_allData_v0
subfolder = 'Tseq-CCW\';
% bandpass
f1 = 65;
f2 = 100;
s1 = 25;
s2 = 45;
th1 = 4;
th2 = 12;
% frequency window length in TFR
delta_f=5;
win_len=0.2;
Fs = 2000;

for ns = 11:isession
    
    path_ns = path{ns};
    cd(path_ns);
    disp(path_ns)
    outfolder = 'Tseq-CCW\gammaWin\';
    mkdir(outfolder)
    
    % step 1 Load spike firing data
    file_input1 = [subfolder 'Cells_singleLap_v2_vel_0.mat'];  % use to get all spikes
    load(file_input1)
    
    file_input2 = trackdata{ns};
    load(file_input2,'Ts_prerunning','Ts_sample','Ts_test','Ts_posttest')
    % read timestamps
    % Ts_start_stop{1,1} = Ts_prerunning{1,1}./ 1000000;%CW
    Ts_start_stop{1,1} = Ts_prerunning{1,2}./ 1000000;%CCW
    % Ts_start_stop{1,2} = Ts_sample{1,1}(:,[1,3])./ 1000000;
    % Ts_start_stop{1,3} = Ts_test{1,1}(:,[1,3])./ 1000000;
    % Ts_start_stop{1,4} = Ts_posttest{1,1}(:,[1,3])./ 1000000;
    
    % step 2 load trackdata
    % load(trackdata{ns})
    
    
    % open TTlist
    celllist = TTList0;
    fid=fopen(celllist);
    if (fid == -1)
        warning([ 'Could not open tfile ' celllist]);
    else
        % read the file names from the t-file list
        TT0 = ReadFileList(celllist);
        numCells0 = length(TT0);
        if numCells0==1 && max(TT0{1}==-1)
            % no cells in ttlist
            numcells=0;
        else
            numcells=numCells0;
        end
    end
    
    % create the CSC file list
    CSC0 = cell(Ncell,1);
    for ii = 1:Ncell
        TTind = strfind(  TT0{ii}, '_' );
        CSC0{ii,1} = ['CSC' TT0{ii}(3:TTind-1) '.ncs'] ;
    end
    nseg = 1;
    
    spikePhase = cell(Ncell,1);
    spikePhase_alllap = cell(Ncell,1,1);
    spikePhase_singlap = cell(Ncell,1,1);
    plFeature = cell(1,3);
    plFeature_alllap = cell(1,3);
    plFeature_singlap = cell(1,3);
    
    Current = [];
    
    fprintf(['=====session %.3g start=====' '\n'], ns)
    for nc = 1:Ncell
        
        if nc == 1 || (nc>=2 && (~strcmp(Current, CSC0{nc,1})))
            Current = CSC0{nc,1};
            disp(Current)
            [sample,tt,tt_raw] = loadCSC_new_WN(Current); % sample unit:uv, tt and tt_raw unit: s
            
            [sigTheta,sigTheta_filtwts]= eegfilt(sample', 2000, th1, th2,length(sample'),[],0,'fir1',0); % pop-up window mode
            [sigGammaS,sigGammaS_filtwts]= eegfilt(sample', 2000, s1, s2,length(sample'),[],0,'fir1',0);
            [sigGammaF,sigGammaF_filtwts]= eegfilt(sample', 2000, f1, f2,length(sample'),[],0,'fir1',0);
            
            % detect slow and fast gamma
            [slow_gamma_windows_EEG,slow_gamma_windows_ts,slow_gamma_windows_TFR_z,...
                slow_gamma_windows_TFR_z0,slow_gamma_windows_bp,start2_slow,stop2_slow,...
                fast_gamma_windows_EEG,fast_gamma_windows_ts,fast_gamma_windows_TFR_z,...
                fast_gamma_windows_TFR_z0,fast_gamma_windows_bp,start2_fast,stop2_fast]...
                = detect_gamma_window_bypeak_sg_fg_v2(tt,sample, sigGammaS,sigGammaF,s1,s2,f1, f2, Fs,delta_f,win_len);
            
            fprintf('\n')
            % 找到pre时候的Windows
            for nseg=1%%%only pre-running
                Ts_start_stop0 = Ts_start_stop{nseg};
                nlap = size(Ts_start_stop0,1);
                IND_sg = [];IND_fg = [];
                for nl = 1:nlap
                    ind_sg = find(slow_gamma_windows_ts(:,1)>=Ts_start_stop0(nl,1) & ...
                        slow_gamma_windows_ts(:,end)<=Ts_start_stop0(nl,2));
                    ind_fg = find(fast_gamma_windows_ts(:,1)>=Ts_start_stop0(nl,1) & ...
                        fast_gamma_windows_ts(:,end)<=Ts_start_stop0(nl,2));
                    
                    IND_sg = [IND_sg; ind_sg];
                    IND_fg = [IND_fg; ind_fg];
                end
            end
            EEG_sgwin = slow_gamma_windows_EEG(IND_sg,:);
            SG_sgwin = slow_gamma_windows_bp(IND_sg,:);
            Win_sg = slow_gamma_windows_ts(IND_sg,:);
            
            EEG_fgwin = fast_gamma_windows_EEG(IND_fg,:);
            FG_fgwin = fast_gamma_windows_bp(IND_fg,:);
            Win_fg = fast_gamma_windows_ts(IND_fg,:);
            %开始画图
            for iw = 1:length(IND_sg)%slow gamma
                ispk = find(spikes{nc,1}>Win_sg(iw,1) & spikes{nc,1}<Win_sg(iw,end));
                if isempty(ispk) || length(ispk)<3
                    continue
                end
                [spikesEeg] = SpikeTStoEEGind(spikes{nc,1}(ispk), Win_sg(iw,:));
                plot_spk_in_wins(spikesEeg,EEG_sgwin(iw,:),SG_sgwin(iw,:),Win_sg(iw,:))
                subplot(2,1,1)
                title(['ccw-' TT0{nc} '-SGwin' num2str(iw)])
                saveas(gcf,[outfolder 'ccw-' TT0{nc} '-SGwin' num2str(iw) '.png'])
                saveas(gcf,[outfolder 'ccw-' TT0{nc} '-SGwin' num2str(iw) '.eps'],'epsc')
            end
            for iw = 1:length(IND_fg)%fast gamma
                ispk = find(spikes{nc,1}>Win_fg(iw,1) & spikes{nc,1}<Win_fg(iw,end));
                if isempty(ispk) || length(ispk)<3
                    continue
                end
                [spikesEeg] = SpikeTStoEEGind(spikes{nc,1}(ispk), Win_fg(iw,:));
                plot_spk_in_wins(spikesEeg,EEG_fgwin(iw,:),FG_fgwin(iw,:),Win_fg(iw,:))
                subplot(2,1,1)
                title(['ccw-' TT0{nc} '-FGwin' num2str(iw)])
                saveas(gcf,[outfolder 'ccw-' TT0{nc} '-FGwin' num2str(iw) '.png'])
                saveas(gcf,[outfolder 'ccw-' TT0{nc} '-FGwin' num2str(iw) '.eps'],'epsc')
            end
            close all
        else
            if nc>=2 && strcmp(Current, CSC0{nc,1})%同一个电极上的cell直接画图
                disp(Current)
                % 直接画图
                for iw = 1:length(IND_sg)%slow gamma
                    ispk = find(spikes{nc,1}>Win_sg(iw,1) & spikes{nc,1}<Win_sg(iw,end));
                    if isempty(ispk) || length(ispk)<3
                        continue
                    end
                    [spikesEeg] = SpikeTStoEEGind(spikes{nc,1}(ispk), Win_sg(iw,:));
                    plot_spk_in_wins(spikesEeg,EEG_sgwin(iw,:),SG_sgwin(iw,:),Win_sg(iw,:))
                    subplot(2,1,1)
                    title(['ccw-' TT0{nc} '-SGwin' num2str(iw)])
                    saveas(gcf,[outfolder 'ccw-' TT0{nc} '-SGwin' num2str(iw) '.png'])
                    saveas(gcf,[outfolder 'ccw-' TT0{nc} '-SGwin' num2str(iw) '.eps'],'epsc')
                end
                for iw = 1:length(IND_fg)%fast gamma
                    ispk = find(spikes{nc,1}>Win_fg(iw,1) & spikes{nc,1}<Win_fg(iw,end));
                    if isempty(ispk) || length(ispk)<3
                        continue
                    end
                    [spikesEeg] = SpikeTStoEEGind(spikes{nc,1}(ispk), Win_fg(iw,:));
                    plot_spk_in_wins(spikesEeg,EEG_fgwin(iw,:),FG_fgwin(iw,:),Win_fg(iw,:))
                    subplot(2,1,1)
                    title(['ccw-' TT0{nc} '-FGwin' num2str(iw)])
                    saveas(gcf,[outfolder 'ccw-' TT0{nc} '-FGwin' num2str(iw) '.png'])
                    saveas(gcf,[outfolder 'ccw-' TT0{nc} '-FGwin' num2str(iw) '.eps'],'epsc')
                end
                close all
            end
        end
        %         nb = length(sprintf([repmat('>' , 1, round((nc-1)/Ncell*10)) '%.2f'],(nc-1)/Ncell*100));
        %         per = sprintf([repmat('>' , 1, round(nc/Ncell*10)) '%.2f'],nc/Ncell*100);
        %         fprintf(1,[repmat('\b',1,nb+1) '%1$s%2$s'],    per, '%');
    end
    fprintf('\ndone!!\n');
    fclose all;
end


function plot_spk_in_wins(spikesEeg,EEG,EEG_bp,win)
figure('Position',[480 180 350 610],'Visible',0)
subplot(2,1,1)
plot(win,EEG)
ylabel('Amp(μV)')
axis tight
subplot(2,1,2)
plot(win,EEG_bp)
ylabel('Amp(μV)')
xlabel('Time(s)')
axis tight
hold on
scatter(win(spikesEeg),EEG_bp(spikesEeg),200,...
    'LineWidth',1.5,'Marker','|','MarkerEdgeColor','r')
hold off
end

%% 找到离spike发生最近的eeg时间点
function [spikesEeg] = SpikeTStoEEGind(spikeTS, eegTS)
spikesEeg = [];
%eliminate spikes that occur earlier than the first time stamp of the EEG
for k = 1:length(spikeTS)
    % Find closest eeg timestamp to the current spike timestamp
    tDiff = (eegTS-spikeTS(k)).^2;
    [~,eegTS_ind] = min(tDiff);
    spikesEeg = [spikesEeg,eegTS_ind];
end
spikesEeg(spikesEeg<1) = 1;
end