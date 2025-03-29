% to load spike data
% edited by Guo M and Zheng C, 2021/5/17
function Data_getspikes
% clear
TTList= 'TTList_dCA1_pyr.txt';
data_video = 'Data_video.mat';

%% Load spikes data
% Load spike firing data
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
pathS = pwd;% current path
NaFile = cell(1,numcells); % Store non existing file names here
NaFile = NaFile(1,1:0); % what for
S = loadSpikes(TT0,pathS,NaFile);
spikes=cell(numcells,1);
spikes(:,1) = TT0;
for nc=1:numcells
    if ~isa(S{nc},'ts') % Empty cell in this session
        ts = 1e64; % use a single ridicilous time stamp if the cell is silent
    else
        % Convert t-file data to timestamps in second?
        ts = Data(S{nc}) / 10000;
    end
    spikes{nc,2} = ts;
end

%% find rat's current positions for spikes
load(data_video);
for nc=1:numcells
    ts = spikes{nc,2};
    Ind = GetSpikePosind(ts,data_video(:,1));
    spikes{nc,3} = data_video(Ind,2:6);  % pos_x, pos_y, pos_ang, vel(+), angle_vel(+/-)
end

%% Save data
save('Data_spikes.mat','spikes');
end
