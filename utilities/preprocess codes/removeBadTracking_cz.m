% Removes position "jumps", i.e position samples that imply that the rat is
% moving quicker than physical possible. Samples in the "jump" parts are
% set to NaN
%
% (c) Raymond Skjerpeng, 2008.
function [x,y] = removeBadTracking_cz(x,y,threshold)

% Indexes to position samples that are to be removed
remInd = [];

diffX = diff(x);
diffY = diff(y);
diffR = sqrt(diffX.^2 + diffY.^2);
ind = find(diffR > threshold);
if isempty(ind)
    return
end

if ind(end) == length(x)
    offset = 2;
else
    offset = 1;
end

for ii = 1:length(ind)-offset
    if ind(ii+1) == ind(ii)+1
        % A single sample position jump, tracker jumps out one sample and
        % then jumps back to path on the next sample. Remove bad sample.
        remInd = [remInd; ind(ii)+1];
        ii = ii+1;
        continue
    else
        % Not a single jump. 2 possibilities:
        % 1. Tracker jumps out, and stay out at the same place for several
        % samples and then jumps back.
        % 2. Tracker just has a small jump before path continues as normal,
        % unknown reason for this. In latter case the samples are left
        % untouched.
        
        idx = find(x(ind(ii)+1:ind(ii+1))==x(ind(ii)+1));
        if length(idx) == length(x(ind(ii)+1:ind(ii+1)));
            remInd = [remInd; (ind(ii)+1:ind(ii+1))'];
        end
    end
end
% Remove the samples
x(remInd) = nan;
y(remInd) = nan;