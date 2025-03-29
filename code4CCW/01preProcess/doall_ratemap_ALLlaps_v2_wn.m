% edited on 3/22/2017
% For each place cell, calculate ALLlaps ratemap
% use spikes at time points with running speed > 5cm/s
% including all pre-running laps, sample and test trials, and post-test
% Save the ratemaps for all later analyses
% Locations are binned in angle, and angle 0 corresponds to the center of
% animal box
% edited on 3/14/2022 by WN
% file from v2

clear

% directories_allData_v2
directories_allData_v0
TTList0 = 'TTList_dCA1_pyr.txt';
file_output = 'Tseq-CCW\Cells_ALLLaps_v2_vel_0.mat';

% limit the posang_ontrack with running speed > vel_threshold
vel_threshold = 0; % cm/s  %IMPORTANT: change to 0 if do not want to limit

for ns = 1:isession
    path_ns = path{ns};
    cd(path_ns);
    
    Ncell = getnumberofcells_cz_v1(TTList0);
    trackdata_ns = trackdata{ns};
    
    % rate maps for overall session
    [Ratemap_AllLaps,spikes,mapAxis] = ratemap_AllLaps_v2(TTList0,trackdata_ns,vel_threshold);
    
    save(file_output,'TTList0','path_ns','Ncell','trackdata_ns',...
        'Ratemap_AllLaps','spikes','mapAxis');
    
    %% calculate the place fields
    p.minNumBins = 3;
    % Bins with rate at p.fieldTreshold * peak rate and higher will be
    % considered as part of a place field
    % Actually, these parameters are not used- see peakAll variable at approximately line 429, which is
    % set to 1 Hz
    p.fieldTreshold = 0.1;
    % Lowest field rate in Hz.
    p.lowestFieldRate = 1;
    peak_all = 0.5;
    
    % place field for overall ratemap across all running laps
    nFields_AllLaps = nan(1,Ncell);
    fieldProp_AllLaps = cell(1,Ncell);
    for nc = 1:Ncell
        ratemap_nc = Ratemap_AllLaps(:,nc);
        [nFields,fieldProp] = placefield_circ_CZ(ratemap_nc,p,mapAxis,peak_all);
        nFields_AllLaps(1,nc) = nFields;
        fieldProp_AllLaps{1,nc} = fieldProp;
    end
    
    save(file_output,'p','peak_all',...
        'nFields_AllLaps','fieldProp_AllLaps','-append');
    
    %% Return
    cd ../
end

cd E:\code\theta_precession_gamma\code4CCW