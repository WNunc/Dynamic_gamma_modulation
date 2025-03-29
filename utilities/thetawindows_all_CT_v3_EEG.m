function [windows,ind_tet,name_tet] = thetawindows_all_CT_v3_EEG(csclist,TTList,celldir)
% cut the theta phase at the minimum firing phase
% get all theta windows (successive theta cycles)  for the whole recording session

% edited on 07/17/2017
% modified from thetawindows_all_CT_v3
% use loadCSC_new_cz_allrecording.m istead of loadCSC_new_cz_breaktime.m
% use EEG from CSClist, but not TTList

path = pwd;

%get all tetrodes with cells for this day
tets = csclist;

%returns the tetrode with highest overall theta power
%ratio cutoff
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

minwindow = 1/f2;%in seconds
maxwindow = 1/f1;

%this will get the eeg with the highest overall theta power
TFRt = 0;
ind_tet = 0;
for i = 1:size(tets,1)
    eegfile = strcat(celldir,'\','CSC',int2str(tets(i)),'.ncs');
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

% w = getthetawindows(thetadelta,ratiocut,eTSi);%for thetafiltering spikes


bp = fftbandpass(eeg,2000,f1-2,f1,f2,f2+2);%theta
bpd = fftbandpass(eeg,2000,f1d-1,f1d,f2d,f2d+1);%delta

%get theta phase (nned to figure out which theta signal to use...or averaged)
%note the following methods are offset as noted:
tphase = thetaPhase_from_Karel(bp); %0 = peak; 180 = trough; 360 = peak
% tphase = 180/pi*(angle(hilbert(bp))); %-180 = trough; 0 = peak; 180 = trough

%get all spikes - needed to find minimum firing phase of theta for cutting
cd(celldir)
fid=fopen(TTList);
if (fid == -1)
    warning([ 'Could not open tfile ' TTList]);
else
    % read the file names from the t-file list
    TT0 = ReadFileList(TTList);
    numCells0 = length(TT0);
    if numCells0==1 && max(TT0{1}==-1)
        % no cells in ttlist
        numCells=0;
    else
        numCells=numCells0;
    end
end
NaFile = cell(1,numCells); % Store non existing file names here
NaFile = NaFile(1,1:0);
S = loadSpikes(TT0,celldir,NaFile);
spikes=cell(numCells,1);
for nc=1:numCells
    if ~isa(S{nc},'ts') % Empty cell in this session
        ts = 1e64; % use a single ridicilous time stamp if the cell is silent
    else
        % Convert t-file data to timestamps in second
        ts = Data(S{nc}) / 10000;
    end
    tempspikes=ts;
%     spikes{nc,1} = get_spike_times_within_specific_eeg_windows(w(:,1),w(:,2),tempspikes); %to theta filter the spikes
    spikes{nc,1} = tempspikes; %to NOT theta filter the spikes
end
cd ../

%convert spikes into theta phase
spkphases = [];
for i = 1:size(spikes,1)
    spkphases = [spkphases, tphase(SpikeTStoEegIndex(spikes{i,1},eTS,fs))];
end

%find minimum firing point along theta
binsize = 30;
% xph = binsize/2+floor(nanmin(tphase)):  binsize  :ceil(nanmax(tphase))-binsize/2; %this is set up for the hist which centers the bars
xph = binsize/2+0:  binsize  :360-binsize/2; %this is set up for the hist which centers the bars
[~,ind] = min(hist(spkphases,xph));
cutph = xph(ind);

%make a cut point at everyone of these phases, (if the other criteria fit)
%ind = find(tphase == cutph); this doesnt work
temp = tphase - cutph;
ind = [];
for i = 1:length(tphase)-1
    if ~isnan(temp(i)) && ~isnan(temp(i+1))
        if temp(i) == 0
            ind = [ind;i];
        elseif sign(temp(i)) == -1 && sign(temp(i+1)) == 1 && abs(temp(i+1)-temp(i))<binsize-1
            ind = [ind;i];
        end
    end
end

thetacyclewindows = [];
for i = 1:length(ind)-1
   %cycle through and create a start and stop for the window
   starts = ind(i);
   stops = ind(i+1)-1;
   thetacyclewindows = [thetacyclewindows ; [starts,stops]];
end
thetacyclewindows(:,3) = thetacyclewindows(:,2) - thetacyclewindows(:,1); %still in index points

windows = [];
for i = 1:size(thetacyclewindows,1)
    %    %make sure it is during theta occurrence%make sure it is between acceptable length
    %    if mean(thetadelta(thetacyclewindows(i,1):thetacyclewindows(i,2)) ) >= ratiocut...
    %            && (thetacyclewindows(i,3) >= minwindow*fs && thetacyclewindows(i,3) <= maxwindow*fs)
    %        windows = [windows; [eTSi(thetacyclewindows(i,1)),eTSi(thetacyclewindows(i,2))]];
    %    end
    windows0 = [eTSi(thetacyclewindows(i,1)),eTSi(thetacyclewindows(i,2))];
    windows0(1,3) = windows0(1,2)-windows0(1,1);
    windows0(1,4) = mean(thetadelta((thetacyclewindows(i,1)):(thetacyclewindows(i,2))));
    windows = [windows;windows0];
end



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

function wins = getthetawindows(thetadelta,cutoff,eTSi)
higher = 0;
startt = [];
stopp = [];
data = thetadelta;
for i = 1:length(eTSi)
     if data(i) >= cutoff && ~higher
       startt = [startt;eTSi(i)];
       higher = 1;
     end
   if data(i) < cutoff && higher
       stopp = [stopp;eTSi(i)];
       higher = 0;
   end
end
if length(startt) > length(stopp)
    startt(end) = [];
end
wins = [startt,stopp];

