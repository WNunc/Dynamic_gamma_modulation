% 批处理
clear

close all
directories_allData_v2
subfolder = 'Tseq\';
%

%%
for ns = 1:isession
    %%
    path_ns = path{ns};
    cd(path_ns);
    disp(path_ns)
    % step 1.1 Load spike firing data
    file_input = [subfolder 'Cells_singleLap_v2_vel_5.mat'];  % use to get all spikes
    celllist = [path{ns} 'TTList_dCA1_pyr.txt'];
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
    
    cellonCSC = [];
    for ncell = 1:numCells0
        cn = find(TT0{ncell}=='_');
        cellonCSC(ncell,1) = str2num(TT0{ncell}(3:cn-1));
    end
    mth =find(CSClist_CA1{ns} == mode(cellonCSC));
    
    
    %     NaFile = cell(1,numcells); % Store non existing file names here
    %     NaFile = NaFile(1,1:0);
    %
    %     S = loadSpikes(TT0,path{ns}(1:end-1),NaFile);
    S = load(file_input);
    
    
    
    % step 1.2 Load trackingdata
    load(trackdata{ns})
    load([path{ns} subfolder 'Data_angle_ontrack.mat'])
    % step 1.3 Load scores
    Scores_input = [subfolder 'BayesData_CircMap_wn_dt40ms_5ms_1cells_1.mat'];
    load(Scores_input)
    scores = Scores;
    
    
    
    
    %% step 2 calculate spike density
    close all
    phase_total = [];%所有圈的spike density最大相位
    TS_running = [];
    phase_mxspk = {};%每一圈的spike density最大相位
    for nl = 1:5
        get_phase_maxspkdensity;
        phase_total = [phase_total phase_spkdensity];
        TS_running(nl,1:2) = TSlimit;
        saveas(gcf,[subfolder 'decoding_maxspk_phase_cutting_lap' num2str(nl) '.png'])
        saveas(gcf,[subfolder 'decoding_maxspk_phase_cutting_lap' num2str(nl) '.fig'])
    end
    %     for b = [12 15 18 20 24 30 36]%测试用
    figure
    bars = 12;
    %         bars = b;%测试用
    nbars = 3*bars;
    barwidth = 1080/nbars;
    edges = 0:barwidth:1080;
    phase_total3 = [phase_total phase_total+360 phase_total+720];
    h = histogram(phase_total3,'BinEdges',edges);
    value = h.Values;
    bincentre = h.BinEdges(1:nbars)+0.5.*diff(h.BinEdges);
%     phasesmooth = smooth(value,round(0.4*bars),'lowess');
    phasesmooth = smooth(value,7,'lowess');
    hold on
    plot(bincentre,phasesmooth,'r','LineWidth',1.5)
    [~,indp]= max(phasesmooth(bars+1:2*bars));
    max_phase = mod(bincentre(indp),360);
    title(['max spike density phase: ' num2str(max_phase)])
    scatter(bincentre(indp+bars),phasesmooth(indp+bars),'Marker','v','LineWidth',10)
    hold off
    saveas(gcf,[subfolder 'hist_AllLap_thetaphase_',num2str(nbars),'bars.png'])
    save([path{ns},subfolder, 'data_AllLap_thetaphase_',num2str(nbars),'bars.mat'],...
        'max_phase','phase_total','TS_running','phase_mxspk');
    theta_phase_cutting
    %     end %测试用
end