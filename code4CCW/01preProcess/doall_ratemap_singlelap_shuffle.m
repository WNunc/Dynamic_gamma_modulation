% edited on 3/22/2017
% For each place cell, calculate single-lap ratemap
% use spikes at time points with running speed > 0cm/s
% including all pre-running laps, sample and test trials, and post-test
% Save the ratemaps for all later analyses
% Locations are binned in angle, and angle 0 corresponds to the center of
% animal box

clear

% directories_allData_v2
directories_allData_v0
TTList0 = 'TTList_dCA1_pyr.txt';
file_output = 'Tseq-CCW\Cells_singleLap_shuffle_vel_5_v3.mat';

% limit the posang_ontrack with running speed > vel_threshold
vel_threshold = 5; % cm/s  %IMPORTANT: change to 0 if do not want to limit

shuffle_num = 100;

for ns = 1:isession
    path_ns = path{ns};
    cd(path_ns);
    
    % 位置域参数
    p.minNumBins = 3;
    % Bins with rate at p.fieldTreshold * peak rate and higher will be
    % considered as part of a place field
    % Actually, these parameters are not used- see peakAll variable at approximately line 429, which is
    % set to 1 Hz
    p.fieldTreshold = 0.1;
    % Lowest field rate in Hz.
    p.lowestFieldRate = 1;
    peak_all = 0.5;
    fclose all;
    Ncell = getnumberofcells_cz_v1(TTList0);
    trackdata_ns = trackdata{ns};
    
    % shuffled rate maps for single laps
    tic
    [Ratemap_singlelap_sh,spikes,mapAxis,Spatialinfo_sh] = ratemap_singlelap_shuffle(TTList0,trackdata_ns,shuffle_num,vel_threshold);
    %% calculate the place fields
    
    %     % place field for single laps
    %     nFields_singlelap = cell(size(Ratemap_singlelap));
    %     fieldProp_singlelap = cell(size(Ratemap_singlelap));
    %     for nseg = 1:size(Ratemap_singlelap,2)% row we need prerunning so nseg=1
    %         for nl = 1:size(Ratemap_singlelap,1)% column
    %             ratemap = Ratemap_singlelap{nl,nseg};
    %             if ~isempty(ratemap)
    %                 for nc = 1:Ncell
    %                     ratemap_nc = ratemap(:,nc);
    %                     [nFields,fieldProp] = placefield_circ_CZ(ratemap_nc,p,mapAxis,peak_all);
    %                     nFields_singlelap{nl,nseg}(1,nc) = nFields;
    %                     fieldProp_singlelap{nl,nseg}{1,nc} = fieldProp;
    %                 end
    %             end
    %         end
    %     end
    ratemap_shuffled = Ratemap_singlelap_sh;
    %     nFields_shuffled{Nshuf} = nFields_singlelap;
    %     fieldProp_shuffle{Nshuf} = fieldProp_singlelap;
    Spatialinfo_shuffled = Spatialinfo_sh;
    
    save(file_output,'TTList0','path_ns','Ncell','trackdata_ns','mapAxis',...
        'p','peak_all','ratemap_shuffled','Spatialinfo_shuffled');
    toc
    
    %% Return
    %     cd ../
    cd E:\code\theta_precession_gamma\code4CCW
end
