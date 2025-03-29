function [onset, offset,para] = DetectSequenceEvents_pos(pxn,varargin)
%% Extract sequential activity from continuous Bayesian decoding
% is a part of DetectSequenceEvents_cz.m
% for positive slopes only


% Default parameters
maxJump_thr = 10;  % in position bins
timeWin = 10;     % in time bin
timeStep = 1;     % in time bin
Distance_thr = 15; % in position bins
jump_prop_thr = 0.5;

if nargin > 1
    maxJump_thr = varargin{1};
end
if nargin > 2
    timeWin = varargin{2};
end
if nargin > 3
    timeStep = varargin{3};
end
if nargin > 4
    Distance_thr = varargin{4};
end
if nargin > 5
    jump_prop_thr = varargin{5};
end

onset_pos = [];
offset_pos = [];
onset_neg = [];
offset_neg = [];
for ii = 1:timeStep:size(pxn,2)
    range = ii:ii+timeWin-1;
    if max(range)>size(pxn,2)
        break
    end
    emptyBin = isnan(pxn(1,range));
    if any(emptyBin)        
        continue
    end
    [~, decodedPos] = max(pxn(:,range));
    jump = diff(decodedPos);
    jump_prop = sum(jump>0)/length(jump);
    down_prop = sum(jump<0)/length(jump);
    maxJump = nanmax(abs(jump));
    Distance = decodedPos(end)-decodedPos(1);
    
    % positive slopes
    if maxJump <= maxJump_thr && Distance >= Distance_thr && jump_prop >= jump_prop_thr &&...
            down_prop == 0
       onset_pos = cat(1,onset_pos,ii);
       offset_pos = cat(1,offset_pos,range(end));
    end
end

% Combine overlapping events
if length(onset_pos) > 1
    isi = onset_pos(2:end)-offset_pos(1:end-1);
    merge = find(isi <= 0);
    onset_pos(merge+1) = [];
    offset_pos(merge) = [];
end
if length(onset_neg) > 1
    isi = onset_neg(2:end)-offset_neg(1:end-1);
    merge = find(isi <= 0);
    onset_neg(merge+1) = [];
    offset_neg(merge) = [];
end

onset = [onset_pos;onset_neg];
offset = [offset_pos;offset_neg];

para = [];
for ii = 1:length(onset)
    timeWin = offset(ii)-onset(ii)+1;
    [~, decodedPos] = max(pxn(:,onset(ii):offset(ii)));
    jump = abs(diff(decodedPos));
    jump_prop = sum(jump>0)/length(jump);
    maxJump = nanmax(jump);
    Distance = abs(decodedPos(end)-decodedPos(1));
    para(ii,:) = [maxJump,timeWin,timeStep,Distance,jump_prop];
end
