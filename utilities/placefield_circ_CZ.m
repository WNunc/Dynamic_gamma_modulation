% placefield identifies the placefields in the firing map. It returns the
% number of placefields and the location of the peak within each
% placefield.
%
% map           Rate map
% pTreshold     Field treshold
% pBins         Minimum number of bins in a field
% mapAxis       The map axis
function [nFields,fieldProp] = placefield_circ_CZ(map,p,mapAxis,peak_all)

binWidth = mapAxis(2) - mapAxis(1);
nbins = length(mapAxis);

% Counter for the number of fields
nFields = 0;
% Field properties will be stored in this struct array
fieldProp = [];

% Allocate memory to the arrays
N = length(map);
% Array that contain the bins of the map this algorithm has visited
visited = zeros(N,1);
nanInd = isnan(map); % where is NaN in metrix is 1 else 0 logical
visited(nanInd) = 1;
visited2 = visited;

% Go as long as there are unvisited parts of the map left
while ~prod(visited) % when all of the map is NaN = 1
    
    % Find the current maximum
    [peak, peakBin] = nanmax(map);
      
    % Array that will contain the bin positions to the current placefield
    fieldBins = peakBin;

    
    % Check if peak rate is high enough
    if peak < p.lowestFieldRate
        break;
    end
    
    %visited2(map<p.fieldTreshold*peak_all) = 1; %changed by LLC to be the peak for slow and fast gamma together
    %visited2(map<p.fieldTreshold*peak) = 1;
    
    if p.fieldTreshold*peak < peak_all
        peak_all = p.fieldTreshold*peak;
    end
    
    visited2(map < peak_all) = 1;
    
    
    % Find the bins that construct the peak field
    [fieldBins,visited2] = recursiveBins_circ(map, visited2, fieldBins, peakBin, N);
    fieldBins = sort(fieldBins); % sort as min to max

    if length(fieldBins) >= p.minNumBins % Minimum size of a placefield
        nFields = nFields + 1;
        
        % Find centre of mass (com)
        ind = find(diff(fieldBins)>1);
        comX = 0;
        R = 0; % Total rate
        if length(ind) == 0
            % this place field is not across 0 degree border
            for ii = 1:length(fieldBins)
                R = R + map(fieldBins(ii));
                comX = comX + map(fieldBins(ii)) * fieldBins(ii);
            end
            startFieldBin = min(fieldBins);
            stopFieldBin = max(fieldBins);
        elseif  length(ind) == 1
            % this place field is across 0 degree border
            fieldBins_shift1 = [fieldBins(ind+1:end);fieldBins(1:ind)];
            fieldBins_shift2 = [fieldBins(ind+1:end);fieldBins(1:ind)+N];
            for ii = 1:length(fieldBins_shift1)
                R = R + map(fieldBins_shift1(ii));
                comX = comX + map(fieldBins_shift1(ii)) * fieldBins_shift2(ii);
            end
            startFieldBin = fieldBins(ind+1);
            stopFieldBin = fieldBins(ind);
        end
        ind_Com = comX/R;
        ind_Com = mod(ind_Com,nbins);
        if ind_Com == 0
            ind_Com =nbins;
        end
        Com = mapAxis(ceil(ind_Com))-(ceil(ind_Com)-ind_Com)*binWidth;
        
        % Average rate in field
        avgRate = nanmean(map(fieldBins));
        % Peak rate in field
        peakRate = nanmax(map(fieldBins));
        % Size of field
        fieldSize = length(fieldBins) * binWidth;
        % Put the field properties in the struct array
        
        fieldProp = [fieldProp; struct('x_COM',Com,'x_peak',mapAxis(peakBin),'avgRate',avgRate,'peakRate',peakRate,'size',fieldSize,'startBin',startFieldBin,'stopBin',stopFieldBin)];
    end
    
    visited(fieldBins) = 1;
    map(visited == 1) = 0;
end


function [binsX,visited] = recursiveBins_circ(map,visited,binsX,ii,N)
% edited by CZ on 03/16/2017
% for circular track

% If outside boundaries of map -> return.
if ii<1 || ii>N
    ii = mod(ii,N);
    if ii == 0
        ii = N;
    end
end
% If all bins are visited -> return.
if prod(visited)
    return;
end
if visited(ii) % This bin has been visited before
    return;
else
    binsX = [binsX;ii];
    visited(ii) = 1;
    % Call this function again in each of the 2 neighbour bins
    [binsX,visited] = recursiveBins_circ(map,visited,binsX,ii-1,N);
    [binsX,visited] = recursiveBins_circ(map,visited,binsX,ii+1,N);
end