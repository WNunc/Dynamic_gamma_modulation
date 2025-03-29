% edited on 3/22/2017
% For each place cell, calculate single-lap ratemap
% use spikes at time points with running speed > 0cm/s
% including all pre-running laps, sample and test trials, and post-test
% Save the ratemaps for all later analyses
% Locations are binned in angle, and angle 0 corresponds to the center of
% animal box
% edited on 2021年1月7日 by 王宁
% edited on 2021年1月15日 by 王宁
% edited on 2022年4月17日 by 王宁
clear

% directories_allData_v2
directories_allData_v0
TTList0 = 'TTList_dCA1_pyr.txt';
file_output = 'Tseq-CCW\Cells_allsegment_v1_vel_5.mat';

% limit the posang_ontrack with running speed > vel_threshold
vel_threshold = 5; % cm/s  %IMPORTANT: change to 0 if do not want to limit

for ns = 1:isession
    path_ns = path{ns};
    cd(path_ns);
    
    Ncell = getnumberofcells_cz_v1(TTList0);
    trackdata_ns = trackdata{ns};
    
    % rate maps for single laps
    [Ratemap_seg,spikes,mapAxis] = ratemap_Allsegment_v1(TTList0,trackdata_ns,vel_threshold);
    
    save(file_output,'TTList0','path_ns','Ncell','trackdata_ns',...
        'Ratemap_seg','spikes','mapAxis');
    
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
    
    % place field for segment
    nFields_seg = cell(size(Ratemap_seg));
    fieldProp_seg = cell(size(Ratemap_seg));
    for nseg = 1:size(Ratemap_seg,2)
        ratemap = Ratemap_seg{nseg};
        if ~isempty(ratemap)
            for nc = 1:Ncell
                ratemap_nc = ratemap(:,nc);
                [nFields,fieldProp] = placefield_circ_CZ(ratemap_nc,p,mapAxis,peak_all);
                nFields_seg{nseg}(1,nc) = nFields;
                fieldProp_seg{nseg}{1,nc} = fieldProp;
            end
        end
    end
    
    save(file_output,'p','peak_all',...
        'nFields_seg','fieldProp_seg','-append');
    
    %% Return
    cd E:\code\theta_precession_gamma\code4CCW
end