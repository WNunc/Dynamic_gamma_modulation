function [windows,ind_tet,name_tet] = thetawindows_all_CT_v5_EEG(csclist,TTList,celldir,cutph,binsize)
% cut the theta phase at a global theta phase: "cutph"
% get all theta windows (successive theta cycles)  for the whole recording session

% edited on 12/13/2017
% modified from thetawindows_all_CT_v4_EEG
% use the tetrode with max cell number but not max theta/delta ratio
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


%this will get the eeg with the most number of cells
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
fclose(fid);

tet_ncell = [];
tetOld = 0;
ncell0 = 0;
ntet = 0;
for nc = 1:numCells
    ind = find(TT0{nc,1} == '_');
    tetNew = str2num(TT0{nc,1}(3:ind-1));
    if tetNew == tetOld
        ncell0 = ncell0+1;
    else
        tetOld = tetNew;
        ntet = ntet+1;
        ncell0 = 1;
    end
    tet_ncell(ntet,:) = [tetNew,ncell0];
end

ind = find(ismember(tet_ncell(:,1),tets) == 1);
tet_ncell = tet_ncell(ind,:);
[~,ind] = max(tet_ncell(:,2)); ind = ind(1);
ind_tet = find(tets == tet_ncell(ind,1));
name_tet = strcat('CSC',int2str(tets(ind_tet)),'.ncs');
eegfile = strcat(celldir,'\','CSC',int2str(tets(ind_tet)),'.ncs');
disp(strcat('Using: ',eegfile));
[eeg, eTSi, eTS, Fs_raw] = loadCSC_new_cz_allrecording(eegfile);

bp = fftbandpass(eeg,2000,f1-2,f1,f2,f2+2);%theta
bpd = fftbandpass(eeg,2000,f1d-1,f1d,f2d,f2d+1);%delta

%get theta phase (nned to figure out which theta signal to use...or averaged)
%note the following methods are offset as noted:
tphase = thetaPhase_from_Karel(bp); %0 = peak; 180 = trough; 360 = peak
% tphase = 180/pi*(angle(hilbert(bp))); %-180 = trough; 0 = peak; 180 = trough


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
    % windows0(1,4) = mean(thetadelta((thetacyclewindows(i,1)):(thetacyclewindows(i,2))));
    windows0(1,4) = nan;
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

