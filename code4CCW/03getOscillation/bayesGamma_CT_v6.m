function [scores] = bayesGamma_CT_v6(scores,csclist,TTList,spikes,ind_tet,data_video,speedlimit,num)

% Edit on 09/30/2020
% This code is for new behavior paradigm from Rat 139
% Modified on bayesGamma_CT_v5
% Only get gamma information on tetrodes from csclist but not TTList
% Same as v5, zscored gamma power excluding time points with vel < 5 cm/s
% scores: pre-running, sample trials, test trials, and post-test trials
% the type of scores has changed
% num is the number of laps

spikes_v2 = spikes;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = 65;
f2 = 100;
s1 = 25;
s2 = 45;
th1 = 4;
th2 = 12;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen(TTList);
if (fid == -1)
    warning([ 'Could not open tfile ' TTList]);
else
    % read the file names from the t-file list
    TT0 = ReadFileList(TTList);
    numCells0 = length(TT0);
    if numCells0==1 && max(TT0{1}==-1)
        % no cells in ttlist
        numcells=0;
    else
        numcells=numCells0;
    end
end
fclose(fid);
if size(spikes_v2,1) ~= numcells
    warning(['cell number not matched']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%go through each cell and get the gamma info for that tetrode
% read the EEG file of whole recording
tetOld = 0;
neeg = 0;
alleeg = [];
for nc = 1:size(TT0,1)
    if TT0{nc,1}(4) ~= '_';
         tetNew = str2double(TT0{nc,1}(3:4));
    elseif TT0{nc,1}(4) == '_';
         tetNew = str2double(TT0{nc,1}(3));
    else
         disp('problem with something');
    end
    
    if ismember(tetNew,csclist)
        if tetNew ~= tetOld
            neeg = neeg +1;
            fprintf(num2str(neeg));
            %load new eeg stuffs
            [e, eTSi, eTS, Fs_raw] = loadCSC_new_cz_allrecording(strcat(pwd,'\CSC',num2str(tetNew),'.ncs'));
            alleeg(neeg,:) = e';
            tetOld = tetNew;
        end
        spikes_v2{nc,5} = TT0{nc};
        spikes_v2{nc,6} = neeg;
    else
        spikes_v2{nc,5} = TT0{nc};
        spikes_v2{nc,6} = nan;
    end
end

% get theta and gamma powers
powf = [];pows = [];
phsf = [];phss = [];phsth = [];
bpf = [];bps = [];bpth = [];
for neeg = 1:size(alleeg,1)
    e = alleeg(neeg,:)';
%     powf = [powf;  TFR_frequency_band_cz(e,2000,5,f1,f2)];
%     pows = [pows;  TFR_frequency_band_cz(e,2000,5,s1,s2)];
    
    bpf0 = eegfilt(e', 2000, f1, f2,length(e'),[],0,'fir1',0);%fast
    bps0 = eegfilt(e', 2000, s1, s2,length(e'),[],0,'fir1',0);%slow
    bpth0 = eegfilt(e', 2000, th1, th2,length(e'),[],0,'fir1',0);%theta
    bpf = [bpf; bpf0];
    bps = [bps; bps0];
    bpth = [bpth; bpth0];
    
    phsf = [phsf;   angle(hilbert(bpf0))*180/pi+180];%fast trough=0, peaks = 180
    phss = [phss;   angle(hilbert(bps0))*180/pi+180];%slow
    phsth = [phsth; angle(hilbert(-bpth0))*180/pi+180];%theta peaks = 0 and 360
end

% cut the time series into segments and laps
powf_seg = [];pows_seg = [];
phsf_seg = [];phss_seg = [];phsth_seg = [];
bpf_seg = [];bps_seg = [];bpth_seg = [];
t_seg = []; alleeg_seg = [];

nsegment = 1;
n_eeg = zeros(1,nsegment);
for nseg = 1:nsegment
    scores_nseg = scores;%
    for nl = 1:num %
        if ~isempty(scores_nseg{nl,3})
            t_nl = scores_nseg{nl,2};
            ind = find(eTSi >= t_nl(1) & eTSi <= t_nl(2));
            n_eeg(nl,nseg) = length(ind);
            
            t_seg = [t_seg,eTSi(ind)'];
            alleeg_seg = [alleeg_seg,alleeg(:,ind)];
%             powf_seg = [powf_seg,powf(:,ind)];
%             pows_seg = [pows_seg,pows(:,ind)];
            bpf_seg = [bpf_seg, bpf(:,ind)];
            bps_seg = [bps_seg, bps(:,ind)];
            bpth_seg = [bpth_seg, bpth(:,ind)];
            phsf_seg =[phsf_seg, phsf(:,ind)];
            phss_seg =[phss_seg, phss(:,ind)];
            phsth_seg =[phsth_seg, phsth(:,ind)];
        end
    end
end
% get zscored gamma power
powf_seg_z = [];
pows_seg_z = [];
% powf_seg_z = zscore(powf_seg, 0 , 2);
% pows_seg_z = zscore(pows_seg, 0 , 2);
%% save EEG information in each sequences
for nseg = 1:nsegment
    scores_nseg = scores;%
    for nl = 1:num %
        if ~isempty(scores_nseg{nl,3})
            t_nl = scores_nseg{nl,2};
            
            % ==== save the EEG information into scores
            ind = find(t_seg >= t_nl(1) & t_seg <= t_nl(2));
            
            %save EEG time points for the sequence
            scores{nl,10} = t_seg(:,ind);
            
            %save raw traces for the sequence
            scores{nl,11} = alleeg_seg(:,ind);
%             
%             %save TFRs for the sequence
%             scores{nl,12} = powf_seg(:,ind);
%             scores{nl,13} = pows_seg(:,ind);
%             
%             %save zscore TFRs for the sequence
%             scores{nl,14} = powf_seg_z(:,ind);
%             scores{nl,15} = pows_seg_z(:,ind);
            
            scores{nl,12} = [];
            scores{nl,13} = [];
            scores{nl,14} = [];
            scores{nl,15} = [];
            
            % bandpass filter of gamma and theta
            scores{nl,16} = bpf_seg(:,ind);
            scores{nl,17} = bps_seg(:,ind);
            scores{nl,18} = bpth_seg(:,ind);
            
            % save phases
            scores{nl,19} = phsf_seg(:,ind);
            scores{nl,20} = phss_seg(:,ind);
            scores{nl,21} = phsth_seg(:,ind); 
%            
        end
    end
end

fprintf(' done.\n');
