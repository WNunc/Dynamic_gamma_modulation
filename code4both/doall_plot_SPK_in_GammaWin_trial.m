% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 还是画example的图！
% 真实受够了！
% 一个lap画一个图
% 每圈好几根TT，都画出来了

%%
clear
close all
% directories_allData_v2
% directories_allData_v0
directories_allData_v0_allgood
% subfolder = 'Tseq-CCW\';
Directionfolder = {'Tseq\','Tseq-CCW\'};
Phasecutfolder = {'CuttingPhaseCW\','CuttingPhaseCCW\'};
DX = {'CW','CCW'};
case1 = '-ontrack';
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

for ns = [1,2,4:6,8,10:13,15,16]
    
    path_ns = path{ns};
    cd(path_ns);
    disp(path_ns)
    
    for D = 1:2
        goodphase_ns = Seq_cutPhase{ns,D};%***
        subfolder1 = Directionfolder{D};
        subfolder2 = Phasecutfolder{D};
        case3 = num2str(goodphase_ns);
        outfolder = [Directionfolder{D},'gammaWin\'];
        mkdir(outfolder)
        
        % step 1 Load spike firing data
        file_input1 = [subfolder1 'Cells_singleLap_v2_vel_0.mat'];  % use to get all spikes
        load(file_input1)
        file_input2 = strcat(path_ns,subfolder2,'data_AllLap_thetaphasecut_',case3,'.mat');
        load(file_input2,'TS_running')
        
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
        for nl =1:5
            file_input3 = strcat(path_ns,subfolder1,'scores',num2str(nseg),'-',num2str(nl),case1,'_v5.mat');
            load(file_input3)
            Current = [];
            
            fprintf(['=====session %.3g start=====' '\n'], ns)
            CSCn = 0;
            for nc = 1:Ncell
                
                if nc == 1 || (nc>=2 && (~strcmp(Current, CSC0{nc,1})))
                    CSCn = CSCn+1;
                    Current = CSC0{nc,1};
                    disp(Current)
                    % [sample,tt,tt_raw] = loadCSC_new_WN(Current); % sample unit:uv, tt and tt_raw unit: s
                    
                    %                 [sigTheta,sigTheta_filtwts]= eegfilt(sample', 2000, th1, th2,length(sample'),[],0,'fir1',0); % pop-up window mode
                    %                 [sigGammaS,sigGammaS_filtwts]= eegfilt(sample', 2000, s1, s2,length(sample'),[],0,'fir1',0);
                    %                 [sigGammaF,sigGammaF_filtwts]= eegfilt(sample', 2000, f1, f2,length(sample'),[],0,'fir1',0);
                    
                    %             % detect slow and fast gamma
                    %             [slow_gamma_windows_EEG,slow_gamma_windows_ts,slow_gamma_windows_TFR_z,...
                    %                 slow_gamma_windows_TFR_z0,slow_gamma_windows_bp,start2_slow,stop2_slow,...
                    %                 fast_gamma_windows_EEG,fast_gamma_windows_ts,fast_gamma_windows_TFR_z,...
                    %                 fast_gamma_windows_TFR_z0,fast_gamma_windows_bp,start2_fast,stop2_fast]...
                    %                 = detect_gamma_window_bypeak_sg_fg_v2(tt,sample, sigGammaS,sigGammaF,s1,s2,f1, f2, Fs,delta_f,win_len);
                    %
                    fprintf('\n')
                    % detect slow and fast gamma in lap
                    tt = scores{nl,10}';
                    sample = scores{nl,11}(CSCn,:)';
                    sigGammaS = scores{nl,17}(CSCn,:)';
                    sigGammaF = scores{nl,16}(CSCn,:)';
                    EEGlap = find(tt>=TS_running(nl,1) & tt<=TS_running(nl,2));
                    [slow_gamma_windows_EEG,slow_gamma_windows_ts,slow_gamma_windows_TFR_z,...
                        slow_gamma_windows_TFR_z0,slow_gamma_windows_bp,start2_slow,stop2_slow,...
                        fast_gamma_windows_EEG,fast_gamma_windows_ts,fast_gamma_windows_TFR_z,...
                        fast_gamma_windows_TFR_z0,fast_gamma_windows_bp,start2_fast,stop2_fast]...
                        = detect_gamma_window_bypeak_sg_fg_v2(tt,sample, sigGammaS,sigGammaF,s1,s2,f1, f2, Fs,delta_f,win_len);
                    ind_sg = [];ind_fg = [];
                    if ~isempty(slow_gamma_windows_ts)
                    ind_sg = find(slow_gamma_windows_ts(:,1)>=TS_running(nl,1) & ...
                        slow_gamma_windows_ts(:,end)<=TS_running(nl,2));
                    end
                    if ~isempty(fast_gamma_windows_ts)
                    ind_fg = find(fast_gamma_windows_ts(:,1)>=TS_running(nl,1) & ...
                        fast_gamma_windows_ts(:,end)<=TS_running(nl,2));
                    end
                    %开始画图
                    F1 = figure('Position',[0 180 3020 750],'Renderer','painters','Visible',0);
                    ttst = tt(EEGlap)-tt(EEGlap(1));
                    axes('Position', [0.050 0.800 0.9 0.11])
                    plot(ttst,sample(EEGlap));xticklabels({})
                    axis tight;ylim([-500,500])
                    axes('Position', [0.050 0.682 0.9 0.11])
                    plot(ttst,sigGammaS(EEGlap));xticklabels({})
                    axis tight;ylim([-150,150])
                    axes('Position', [0.050 0.564 0.9 0.11])
                    plot(ttst,sigGammaF(EEGlap));xticklabels({})
                    axis tight;ylim([-100,100])
                    xmax = max(xlim);
                    axes('Position', [0.050 0.1 0.9 0.45])
                    for c = 1:Ncell
                        ispk = find(spikes{c,1}>TS_running(nl,1) & spikes{c,1}<TS_running(nl,2));
                        if isempty(ispk) || length(ispk)<3
                            continue
                        end
                        
                        [spikesEeg] = SpikeTStoEEGind(spikes{c,1}(ispk), tt(EEGlap));
                        %             subplot(5,1,[4,5])
                        scatter(ttst(spikesEeg),ones(1,length(spikesEeg))*c,'Marker','|','MarkerEdgeColor','k')
                        hold on
                    end
                    xlim([0,xmax]);
                    ymax = max(ylim);
                    plot_fill_gammawin(slow_gamma_windows_ts(ind_sg,:),fast_gamma_windows_ts(ind_fg,:),ymax,tt(EEGlap))
                    SSS = sprintf('%s lap-%s %s',DX{D},num2str(nl),Current);
                    annotation('textbox',[0.32,0.92,0.5,0.05],'LineStyle','none','String',SSS,...
                        'FontSize',20)
                    saveas(F1,[outfolder,SSS(1:end-4),'.png']);
                    saveas(F1,[outfolder,SSS(1:end-4)],'epsc');
                    clf
                end
                
                %         nb = length(sprintf([repmat('>' , 1, round((nc-1)/Ncell*10)) '%.2f'],(nc-1)/Ncell*100));
                %         per = sprintf([repmat('>' , 1, round(nc/Ncell*10)) '%.2f'],nc/Ncell*100);
                %         fprintf(1,[repmat('\b',1,nb+1) '%1$s%2$s'],    per, '%');
            end
        end
    end
    fprintf('\ndone!!\n');
    fclose all;
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

%% 画出gamma window
function plot_fill_gammawin(slow_gamma_windows_ts,fast_gamma_windows_ts,ymax,tt)
hold on
% sgamma win
n = size(slow_gamma_windows_ts, 1);
% 循环遍历每个window
for i = 1:n
    % 获取当前window的数据
    window = slow_gamma_windows_ts(i, :);
    window = window-tt(1);
    fx = [window(1),window(end),window(end),window(1)];
    fy = [0,0,ymax,ymax];
    % 绘制矩形并填充颜色
    f = fill(fx,fy, 'k');
    f.FaceAlpha = 0.2;
    f.EdgeAlpha = 0;
end
% fgamma win
n = size(fast_gamma_windows_ts, 1);
for i = 1:n
    % 获取当前window的数据
    window = fast_gamma_windows_ts(i, :);
    window = window-tt(1);
    fx = [window(1),window(end),window(end),window(1)];
    fy = [0,0,ymax,ymax];
    % 绘制矩形并填充颜色
    f = fill(fx,fy, 'r');
    f.FaceAlpha = 0.1;
    f.EdgeAlpha = 0;
    
end

end
