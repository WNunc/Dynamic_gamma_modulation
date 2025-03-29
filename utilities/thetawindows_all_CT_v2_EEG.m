function [windows,ind_tet,name_tet] = thetawindows_all_CT_v2_EEG(csclist)

% cut the theta phase at theta peaks
% get all theta windows (successive theta cycles)  for the whole recording session

% edited on 07/17/2017
% modified from thetawindows_all_CT_v2
% use loadCSC_new_cz_allrecording.m istead of loadCSC_new_cz_breaktime.m
% use EEG from CSClist, but not TTList

path = pwd;


%get all tetrodes whose EEG cound be used
tets = csclist;

%returns the tetrode with highest overall theta power ratio cutoff
ratiocut = 3;
%smooth ratio by
smo = 1000; %in index points 
%theta
f1 = 6;
f2 = 10;
%delta
f1d = 2;
f2d = 4;
%sampling of eeg
fs = 2000;


%this will get the eeg with the highest overall theta power
TFRt = 0;
ind_tet = 0;
for i = 1:size(tets,1)
    eegfile = strcat(path,'\','CSC',int2str(tets(i)),'.ncs');
    disp(strcat('Loading Tetrode:',eegfile));
    
    [neeg, eTSi, eTS, Fs_raw] = loadCSC_new_cz_allrecording(eegfile);
    % calculate the TFR for each segment in EEG
    ind = find(diff(eTS) > 1)';
    ind_segment = [[1;ind+1],[ind;length(eTS)]];
    t_segment = eTS(ind_segment);
    nTFRt = [];
    for nseg = 1:size(t_segment,1)
        [~,i_start] = min(abs(eTSi-t_segment(nseg,1)));
        [~,i_stop] = min(abs(eTSi-t_segment(nseg,2)));
        nTFRt0 = TFR_frequency_band(neeg(i_start:i_stop),fs,5,f1,f2); %theta
        nTFRt = [nTFRt,nTFRt0];
    end

    if mean(nTFRt) > mean(TFRt)
        TFRt = nTFRt;
        eeg = neeg;
        ind_tet = i;
        name_tet = strcat('CSC',int2str(tets(i)),'.ncs');
        disp(strcat('Using: ',eegfile));
    end
end
clear nTFRt neeg

%get periods of time when theta is happening:
eegfile = strcat(path,'\','CSC',int2str(tets(ind_tet)),'.ncs');
[eeg, eTSi, eTS, Fs_raw] = loadCSC_new_cz(eegfile);
TFRt = TFR_frequency_band(eeg,fs,5,f1,f2); %theta
TFRd = TFR_frequency_band(eeg,fs,5,f1d,f2d);%delta
thetadelta = smooth(TFRt./TFRd,smo);

% cut all of the theta cycles out, using the theta peaks
bp_theta = fftbandpass(eeg,2000,f1-1,f1+1,f2-1,f2+1);%theta
phs_theta = angle(hilbert(-bp_theta))*180/pi+180;  %0 = peak; 180 = trough; 360 = peak
% phs_theta = thetaPhase_from_Karel(bp_theta); %0 = peak; 180 = trough; 360 = peak

ind_theta_cycles = get_theta_cycles(phs_theta);
windows(:,1) = eTSi(ind_theta_cycles(:,1));
windows(:,2) = eTSi(ind_theta_cycles(:,2));
windows(:,3) = windows(:,2)-windows(:,1);

% for each theta cycle, get the mean theta/delta ratio
for nw = 1:size(windows,1)
    windows(nw,4) = mean(thetadelta(ind_theta_cycles(nw,1):ind_theta_cycles(nw,2)));
end


% ==================================================================
function tets = gettetnumbers(file)

tets = nan(numel(textread(file,'%1c%*[^\n]')),1); 

fid = fopen(file,'r');
if fid == -1
    msgbox('Could not open the input file','ERROR');
end

for i = 1:size(tets,1)
  tline = fgetl(fid);
  if ~ischar(tline) 
      break 
  else
      if tline(4)=='_'
        tets(i) = str2double(tline(3));
      elseif tline(5) == '_'
        tets(i) = str2double(tline(3:4));
      end
  end          
end

fclose(fid);

tets(isnan(tets)) = [];
%tets(tets>6) = [];
% tets = mode(tets);
tets = unique(tets);
% tets = num2cell(tets);
